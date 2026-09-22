"""
Closed-loop linear temperature ramp for the Ferrovac PBN heating stage.

Reads a type-K thermocouple through an HP/Agilent 34401A over GPIB and drives
a BK Precision 1685B DC supply over USB/serial to follow a linear temperature
ramp, then soak, then a controlled cooldown.

Control law is feedforward + PI.  The feedforward curve came out of the
28 Aug 2026 heating run: steady state is radiation-dominated, so the voltage
needed to hold a temperature goes roughly as (T_K^4 - T_amb_K^4)^(1/2).  The
PI only trims the residual, which kept that fit inside about +/- 1.4 V.
No derivative term -- dead time is 40-80 s and D would just amplify noise.

Hard limits: 40.0 V, 3.00 A, 120 W.  Those are the element's nameplate
ratings.  Power is computed from the supply's own GETD readback every loop,
never inferred from the commanded voltage, because the heater's resistance
changes a lot between cold and hot.

Requires: pip install pyvisa matplotlib

Run from a terminal:  python ramp_controller.py
"""

import atexit
import csv
import math
import os
import queue
import threading
import time
import tkinter as tk
from datetime import datetime
from tkinter import messagebox, ttk

import matplotlib.dates as mdates
import pyvisa
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

# ================================================================== config

DMM_ADDRESS = "GPIB0::6::INSTR"
PSU_ADDRESS = "ASRL7::INSTR"

READ_INTERVAL_S = 1.0        # thermocouple sampling period
CONTROL_EVERY_N = 3          # run the control law every N readings
DEFAULT_COLD_JUNCTION_C = 22.0

SAVE_DIR = os.path.join(os.path.expanduser("~"), "tc_logs")

# ---- hard limits (element nameplate: 40 V / 3 A / 120 W) ----
V_MAX = 40.0                 # volts, absolute ceiling on the command
I_LIMIT_A = 3.00             # programmed into the supply as CURR300
P_MAX_W = 120.0              # watts, enforced against the GETD readback
I_TRIP_A = 3.10              # readback above this is a fault
P_TRIP_W = 126.0             # readback above this is a fault

# ---- controller ----
V_SLEW_UP_V_PER_MIN = 4.0    # how fast the command may rise
V_SLEW_DOWN_V_PER_MIN = 20.0 # dropping power is always safe -- let it fall fast
R_SANITY_FACTOR = 2.0        # live resistance outside model/2..model*2 is distrusted
COOLDOWN_OFF_AFTER_S = 60.0  # at the voltage floor this long in cooldown -> output off
# Gains are SCHEDULED from the plant model, not fixed.  The stage's time
# constant runs from ~900 s at 50 C to ~46 s at 800 C because radiative losses
# climb so steeply, so no single Kp/Ki is right across the range.  The GUI
# fields are multipliers on the scheduled values: 1.0 means "use the model".
KP_DEFAULT = 1.0             # multiplier on the scheduled Kp
KI_DEFAULT = 1.0             # multiplier on the scheduled Ki
DEAD_TIME_S = 60.0           # measured 40-80 s on 28 Aug
LOOP_TC_S = 60.0             # closed-loop time constant target (SIMC tau_c)

# Look-ahead feedforward.  Power applied now reaches the sample roughly one
# dead time later, and the ramp is planned, so feedforward aims at where the
# setpoint WILL be LOOKAHEAD_S from now.  Near the target that drops the power
# to holding level a dead time early -- cutting off exactly the in-flight heat
# that caused 5-14 C of overshoot -- without slowing the setpoint down.
# 75 s is deliberately above the 60 s nominal: step responses on 28 Aug gave
# 38, 71 and 82 s, and erring long costs a little time while erring short
# costs overshoot.
LOOKAHEAD_S = 75.0

# Short final taper, constant deceleration (v^2 = 2 a d) so it arrives in
# finite time.  Starts TAPER_DEAD_TIMES dead-times' worth of distance out.
# Its extra time, TAPER_DEAD_TIMES * DEAD_TIME_S, is folded back into the ramp
# rate, so the ramp time you enter is the time you get.
# 22 Sep: at 1.0 the taper plus the 75 s look-ahead cut feedforward power over
# the last ~65 C of a 150 C / 5 min ramp -- more than half of it -- which is
# what bent the curve.  0.5 keeps the arrival clean with half the rounding.
TAPER_DEAD_TIMES = 0.5
TAPER_FLOOR = 0.05
I_TERM_CLAMP_V = 5.0         # integrator authority, volts

# Heat capacity of the stage + holder, J/K.  Estimated from the 28 Aug step
# responses (C = tau * dP/dT).  This sets the ramp-rate feedforward term:
# climbing needs P_loss(T) + C * dT/dt, and below ~300 C the second term is
# almost all of it.  Raise it if the ramp lags, lower it if it overshoots.
# 22 Sep hardware run (62 -> 150 C): energy balance gives C = 25-27 J/K for
# dead times of 60-90 s.  At 20 the ramp ran ~9 C behind, pinned on the leash.
HEAT_CAPACITY_J_PER_K = 26.0
LEASH_C = 10.0               # stage this far AHEAD of the plan -> move the plan up
# Stage this far BEHIND the plan -> the setpoint waits for it.  Was 10 C (shared
# with LEASH_C); 22 Sep the ramp sat 9 C behind without the leash ever firing.
LEASH_BEHIND_C = 3.0

# The start is released on the stage's measured response, not a fixed timer.
# 22 Sep (108 -> 175 C) the stage was still cooling at Start; the fixed 60 s
# delay released the setpoint while the stage was barely moving and a lag
# opened immediately.  The setpoint now stays parked until the measured slope
# reaches START_SLOPE_FRACTION of the planned rate (least-squares over
# SLOPE_WINDOW_S), with START_MIN_WAIT_S / START_MAX_WAIT_S as bounds.  On
# release the setpoint is rebased onto the measured temperature.  Extra
# waiting lengthens the ramp rather than steepening it.
START_SLOPE_FRACTION = 0.6
SLOPE_WINDOW_S = 21.0
START_MIN_WAIT_S = 30.0
START_MAX_WAIT_S = 240.0

# ---- ramp envelope ----
# 40 V holds about 828 degC at steady state, so a target near that leaves the
# PI no headroom.  Cap well below it.
T_TARGET_MAX_C = 780.0
T_OVERTEMP_MARGIN_C = 30.0   # trip if we exceed target by this much
MAX_RUN_HOURS = 8.0

COOLDOWN_RATE_C_PER_MIN = 40.0
COOLDOWN_END_C = 100.0       # below this, cut power and stop

# ---- fault tolerances ----
MAX_CONSECUTIVE_READ_FAILS = 3
MAX_CONSECUTIVE_PSU_FAILS = 3
V_READBACK_TOLERANCE = 1.0   # volts of disagreement between command and GETD

# =========================================================== ITS-90 type K

_C_POS = [
    -0.176004136860e-01, 0.389212049750e-01, 0.185587700320e-04,
    -0.994575928740e-07, 0.318409457190e-09, -0.560728448890e-12,
    0.560750590590e-15, -0.320207200030e-18, 0.971511471520e-22,
    -0.121047212750e-25,
]
_A0, _A1, _A2 = 0.118597600000e00, -0.118343200000e-03, 0.126968600000e03

_INV_NEG = [
    0.0, 2.5173462e01, -1.1662878e00, -1.0833638e00, -8.9773540e-01,
    -3.7342377e-01, -8.6632643e-02, -1.0450598e-02, -5.1920577e-04,
]
_INV_MID = [
    0.0, 2.508355e01, 7.860106e-02, -2.503131e-01, 8.315270e-02,
    -1.228034e-02, 9.804036e-04, -4.413030e-05, 1.057734e-06, -1.052755e-08,
]
_INV_HIGH = [
    -1.318058e02, 4.830222e01, -1.646031e00, 5.464731e-02,
    -9.650715e-04, 8.802193e-06, -3.110810e-08,
]

EMF_MIN_MV, EMF_MAX_MV = -5.891, 54.886


def _poly(coeffs, x):
    return sum(c * x**i for i, c in enumerate(coeffs))


def celsius_to_mv(t_c):
    """Type-K EMF in mV for a junction at t_c degrees C."""
    return _poly(_C_POS, t_c) + _A0 * math.exp(_A1 * (t_c - _A2) ** 2)


def mv_to_celsius(e_mv):
    """Type-K temperature in C for a total EMF of e_mv millivolts."""
    if e_mv < EMF_MIN_MV:
        raise ValueError(f"EMF {e_mv:.4f} mV below type-K range")
    if e_mv < 0.0:
        return _poly(_INV_NEG, e_mv)
    if e_mv < 20.644:
        return _poly(_INV_MID, e_mv)
    if e_mv <= EMF_MAX_MV:
        return _poly(_INV_HIGH, e_mv)
    raise ValueError(f"EMF {e_mv:.4f} mV above type-K range")


def measured_to_celsius(measured_mv, cold_junction_c):
    """Convert the millivolts the DMM reads into a hot-junction temperature."""
    return mv_to_celsius(measured_mv + celsius_to_mv(cold_junction_c))


# ============================================ feedforward: temperature -> V

# Anchors below 550 C come from the radiative model sqrt(c*(T_K^4 - Ta_K^4));
# anchors at and above 600 C are settled points measured on 28 Aug 2026.
# Only monotonicity and rough accuracy matter -- the PI covers the rest.
_FF_T = [21.0, 100.0, 200.0, 300.0, 400.0, 500.0,
         608.7, 733.4, 757.3, 788.3, 799.0, 828.3]
_FF_V = [0.0, 3.59, 6.79, 10.43, 14.64, 19.46,
         24.00, 32.00, 34.00, 36.50, 37.50, 40.00]


# Heater resistance against temperature, fitted to the 16 Sep ramp plus the
# cold bench measurement.  Pyrolytic graphite has a strong negative TCR: the
# element is 28.7 ohm cold and settles to 15.4 ohm above about 250 C.  This
# matters because the ramp-power term has to be converted to volts, and using
# the cold value up there would overstate the voltage by nearly 40%.
_R_INF, _R_AMP, _R_TAU = 15.40, 18.59, 74.6


def resistance_ohm(t_c):
    return _R_INF + _R_AMP * math.exp(-max(t_c, 0.0) / _R_TAU)


R_NOMINAL_OHM = resistance_ohm(22.0)

# Steady-state loss power against temperature, in watts.  Below 300 C these
# come from the 16 Sep ramp fit; at and above 600 C they are the settled
# points of the 28 Aug run converted through V^2 / R(T).  Between 300 and 600
# there is no data, so it interpolates -- the PI covers that stretch.
_P_T = [22.0, 100.0, 200.0, 300.0, 608.7, 733.4, 757.3, 788.3, 799.0, 828.3]
# 22 Sep: points at 300 C and below scaled x1.2.  The 16 Sep fit ran ~20% low:
# 150 C held at 3.70 W (table 3.07), 175 C still short at 4 W / 8.9 V (table
# 3.73), and a full-run energy balance gave a loss scale of 1.17.
_P_W = [0.0, 2.08, 5.28, 9.49, 37.30, 66.50, 75.10, 86.50, 91.30, 103.90]


def loss_power_w(t_c):
    """Power needed to hold t_c at steady state."""
    if t_c <= _P_T[0]:
        return 0.0
    if t_c >= _P_T[-1]:
        return _P_W[-1]
    for i in range(1, len(_P_T)):
        if t_c <= _P_T[i]:
            frac = (t_c - _P_T[i - 1]) / (_P_T[i] - _P_T[i - 1])
            return _P_W[i - 1] + frac * (_P_W[i] - _P_W[i - 1])
    return _P_W[-1]


def feedforward_volts(t_c, rate_c_per_s=0.0, resistance=None,
                      heat_capacity=HEAT_CAPACITY_J_PER_K):
    """
    Voltage needed to sit at t_c AND climb at rate_c_per_s.

        P = P_loss(T) + C * dT/dt        V = sqrt(P * R(T))

    Feedforward is computed in POWER and converted to volts at the end, which
    is the only way to get it right when the heater's resistance moves by a
    factor of two across the working range.  Pass a live resistance from the
    GETD readback when one is available; otherwise the fitted curve is used.
    """
    p = loss_power_w(t_c) + heat_capacity * max(rate_c_per_s, 0.0)
    r = resistance if (resistance and resistance > 1.0) else resistance_ohm(t_c)
    return min(V_MAX, math.sqrt(max(0.0, p * r)))


def scheduled_gains(t_c, v_op, heat_capacity=HEAT_CAPACITY_J_PER_K):
    """
    SIMC PI gains at temperature t_c, linearised about operating voltage v_op.

    Plant: C dT/dt = V^2/R - P_loss(T), plus dead time L.  With G = dP_loss/dT,
    the time constant is tau = C/G and the gain in C/V is (2V/R)/G.  G cancels
    out of the SIMC proportional gain, which leaves

        Kp = C R / (2 V (tau_c + L))        Ki = Kp / min(tau, 4 (tau_c + L))
    """
    v = max(v_op, 3.0)
    tl = LOOP_TC_S + DEAD_TIME_S
    kp = heat_capacity * resistance_ohm(t_c) / (2.0 * v * tl)
    kp = min(max(kp, 0.02), 0.35)
    g = max((loss_power_w(t_c + 5.0) - loss_power_w(t_c - 5.0)) / 10.0, 0.01)
    ti = min(heat_capacity / g, 4.0 * tl)
    return kp, kp / ti


def max_rate_c_per_min(t_c):
    """Roughly what the stage managed on 28 Aug, by temperature band."""
    if t_c < 500.0:
        return 45.0
    if t_c < 700.0:
        return 28.0
    if t_c < 780.0:
        return 15.0
    return 8.0


# ================================================================ PSU layer


class PsuError(Exception):
    pass


_LIVE_PSU = None   # most recent supply object, for emergency_off()


class Bk1685b:
    """
    BK 1685B over serial.  Protocol quirks worth remembering:

      VOLTxxx  sets voltage, xxx = volts * 10   (0.1 V steps)
      CURRxxx  sets current limit, xxx = amps * 100
      SOUT0    output ON      SOUT1    output OFF   (yes, inverted)
      GETD     returns VVVVIIIIS  -> volts*100, amps*100, 0=CV 1=CC
      GETS     returns VVVCCC     -> the programmed setpoints

    The supply's minimum output is 1 V.  VOLT000 is out of range and the unit
    answers it with silence rather than OK, so 1.0 V is treated as zero here
    and SOUT1 is what actually cuts power.

    Set commands are fire-and-verify: we do not trust the OK, we read the
    state back with GETS / GETD and check it.  A missing OK on a write is
    logged, not fatal.
    """

    V_FLOOR = 1.0   # supply cannot go below this

    def __init__(self, address=PSU_ADDRESS, log=print):
        self.log = log
        self.rm = pyvisa.ResourceManager()
        self.dev = self.rm.open_resource(address)
        self.dev.baud_rate = 9600
        self.dev.write_termination = "\r"
        self.dev.read_termination = "\r"
        self.dev.timeout = 2000
        self.output_on = False
        global _LIVE_PSU
        _LIVE_PSU = self

    def _drain(self, deadline_s=2.5):
        """Collect response lines until OK or timeout. Never raises."""
        lines, saw_ok = [], False
        deadline = time.time() + deadline_s
        while time.time() < deadline:
            try:
                line = self.dev.read().strip()
            except Exception:
                break
            if line == "OK":
                saw_ok = True
                break
            if line:
                lines.append(line)
        return lines, saw_ok

    def _write(self, text):
        """Send a set command. Tolerates a missing OK."""
        self.dev.write(text)
        _, saw_ok = self._drain()
        if not saw_ok:
            self.log(f"note: {text} returned no OK (continuing)")
        return saw_ok

    def _query(self, text):
        """Send a read command. Raises if no payload comes back."""
        self.dev.write(text)
        lines, _ = self._drain()
        if not lines:
            raise PsuError(f"{text}: no response")
        return lines[0]

    # ------------------------------------------------------------ state

    def get_setpoints(self):
        """Returns (volts, amps) as programmed. GETS payload is VVVCCC."""
        raw = self._query("GETS")
        if len(raw) < 6:
            raise PsuError(f"GETS payload malformed: {raw!r}")
        return int(raw[0:3]) / 10.0, int(raw[3:6]) / 100.0

    def read(self):
        """Returns (volts, amps, in_constant_current). GETD is VVVVIIIIS."""
        raw = self._query("GETD")
        if len(raw) < 9:
            raise PsuError(f"GETD payload malformed: {raw!r}")
        return int(raw[0:4]) / 100.0, int(raw[4:8]) / 100.0, raw[8] == "1"

    # ----------------------------------------------------------- control

    def configure(self):
        """Output off, current limit set, voltage at the floor. Verified."""
        self._write("SOUT1")
        self.output_on = False
        self._write(f"CURR{int(round(I_LIMIT_A * 100)):03d}")
        self._write(f"VOLT{int(round(self.V_FLOOR * 10)):03d}")

        v_set, i_set = self.get_setpoints()
        self.log(f"PSU setpoints read back: {v_set:.1f} V, {i_set:.2f} A")
        if abs(i_set - I_LIMIT_A) > 0.02:
            raise PsuError(
                f"Current limit did not take: asked {I_LIMIT_A:.2f} A, "
                f"supply reports {i_set:.2f} A"
            )
        v_out, i_out, _ = self.read()
        if v_out * i_out > 1.0:
            raise PsuError(
                f"Supply still delivering {v_out * i_out:.1f} W after SOUT1"
            )

    def set_voltage(self, volts):
        """Command a voltage. Returns what was actually commanded."""
        volts = max(0.0, min(V_MAX, volts))
        if volts < self.V_FLOOR:
            volts = self.V_FLOOR
        self._write(f"VOLT{int(round(volts * 10)):03d}")
        return round(volts, 1)

    def output(self, on):
        self._write("SOUT0" if on else "SOUT1")
        self.output_on = on

    def shutdown(self):
        """Best effort, never raises -- called from fault paths and atexit."""
        for text in (f"VOLT{int(round(self.V_FLOOR * 10)):03d}", "SOUT1"):
            try:
                self.dev.write(text)
                self._drain(deadline_s=1.0)
            except Exception:
                pass
        self.output_on = False

    def close(self):
        self.shutdown()
        try:
            self.dev.close()
        except Exception:
            pass


# ============================================================ control loop


class RampWorker(threading.Thread):
    """
    Owns both instruments.  Reads the thermocouple every READ_INTERVAL_S and
    updates the heater command every CONTROL_EVERY_N readings.  Everything
    reaches the GUI through a queue; the GUI never touches an instrument.
    """

    def __init__(self, out_q, cmd_q, params):
        super().__init__(daemon=True)
        self.out_q = out_q
        self.cmd_q = cmd_q
        self.p = params

        self.dmm = None
        self.psu = None

        self.state = "CONNECTING"
        self.t_setpoint = None
        self.t_start = None
        self.rate_c_per_s = 0.0
        self.i_term = 0.0
        self.v_cmd = 0.0
        self.last_control_t = None
        self.tick = 0
        self.read_fails = 0
        self.psu_fails = 0
        self.r_live = None
        self.leashed = False
        self.rate_eff = 0.0
        self.kp_now = self.ki_now = 0.0
        self.rebased = False
        self.r_warned = False
        self.floor_since = None
        self.start_delay_left = 0.0
        self.parked = False
        self.park_elapsed = 0.0
        self.hist = []
        self.stop_flag = threading.Event()

    # ---------------------------------------------------------- plumbing

    def say(self, text):
        self.out_q.put(("status", text))

    def fault(self, reason):
        self.state = "FAULT"
        if self.psu:
            self.psu.shutdown()
        self.out_q.put(("fault", reason))

    def drain_commands(self):
        while True:
            try:
                msg = self.cmd_q.get_nowait()
            except queue.Empty:
                return
            kind = msg[0]
            if kind == "start":
                self.begin_ramp(msg[1], msg[2])
            elif kind == "retarget":
                self.retarget(msg[1], msg[2])
            elif kind == "hold" and self.state == "RAMP":
                self.state = "HOLD"
                self.say("Ramp held. Setpoint frozen.")
            elif kind == "resume" and self.state == "HOLD":
                self.state = "RAMP"
                self.say("Ramp resumed.")
            elif kind == "cooldown":
                if self.state in ("RAMP", "HOLD", "SOAK"):
                    self.state = "COOLDOWN"
                    self.say("Controlled cooldown started.")
            elif kind == "estop":
                self.fault("Emergency stop pressed.")
            elif kind == "quit":
                self.stop_flag.set()

    # ---------------------------------------------------------- lifecycle

    def begin_ramp(self, target_c, minutes):
        if self.state not in ("IDLE",):
            self.say("Not idle; ignoring start.")
            return
        if self.last_temp is None or math.isnan(self.last_temp):
            self.say("No valid temperature yet; cannot start.")
            return

        self.t_start = self.last_temp
        self.t_setpoint = self.last_temp
        self.target_c = target_c
        self.start_delay_left = DEAD_TIME_S
        self._park()
        self.rate_c_per_s = self._honored_rate(
            target_c - self.t_start, minutes, self.start_delay_left
        )
        self.i_term = 0.0
        self.leashed = False
        self.last_control_t = time.time()
        self.run_deadline = time.time() + MAX_RUN_HOURS * 3600.0

        # Jump straight to the feedforward voltage instead of slewing up from
        # zero.  This is the open-loop voltage needed to sit at T_start and
        # climb at the requested rate, so it is a known-safe starting point --
        # and creeping up to it at the slew rate wastes minutes of a short
        # ramp.  The slew limit still governs everything after this.
        self.rate_eff = self._tapered_rate()
        self.state = "RAMP"
        sp_f, r_f = self._future_setpoint(LOOKAHEAD_S)
        v0 = feedforward_volts(sp_f, r_f, None, self.p["heat_capacity"])
        try:
            self.v_cmd = self.psu.set_voltage(v0)
            self.psu.output(True)
        except Exception as exc:
            self.fault(f"Could not enable output: {exc}")
            return

        self.state = "RAMP"
        self.say(
            f"Ramping {self.t_start:.1f} -> {target_c:.0f} C, arriving in "
            f"{minutes:.0f} min. Starting at {self.v_cmd:.1f} V."
        )

    def retarget(self, target_c, minutes):
        """
        Change target and/or remaining time without restarting.

        The new ramp starts from the CURRENT setpoint, not the original T0,
        so the command does not jump -- the integrator keeps its charge and
        the feedforward moves smoothly.  A target below the current setpoint
        gives a negative rate, which is a controlled descent: the rate term
        drops out of the feedforward and the loop simply backs the voltage
        off as the setpoint walks down.
        """
        if self.state not in ("RAMP", "HOLD", "SOAK"):
            self.say("Not ramping -- use Start Ramp.")
            return

        base = self.t_setpoint if self.t_setpoint is not None else self.last_temp
        self.target_c = target_c
        self.start_delay_left = DEAD_TIME_S if self.state == "SOAK" else 0.0
        if self.state == "SOAK":
            self._park()
        self.rate_c_per_s = self._honored_rate(
            target_c - base, minutes, self.start_delay_left
        )
        self.rate_eff = self._tapered_rate()
        self.leashed = False
        self.run_deadline = max(self.run_deadline,
                                time.time() + minutes * 60.0 + 1800.0)

        if self.state == "SOAK":
            self.state = "RAMP"
        held = " (still held)" if self.state == "HOLD" else ""
        self.say(
            f"Retargeted: {base:.1f} -> {target_c:.0f} C, arriving in "
            f"{minutes:.0f} min{held}."
        )

    # ------------------------------------------------------ the math bit

    def advance_setpoint(self, dt):
        self.leashed = False
        self.rebased = False
        now = time.time()
        self.hist.append((now, self.last_temp))
        while self.hist and now - self.hist[0][0] > SLOPE_WINDOW_S:
            self.hist.pop(0)

        if self.state == "RAMP" and self.parked:
            # Heat is on but has not reached the sample yet.  Keep the setpoint
            # parked until the stage is visibly climbing at a good fraction of
            # the planned rate, so the setpoint starts moving on the measured
            # line instead of opening a lag in the first minute.
            self.park_elapsed += dt
            self.start_delay_left = max(0.0, DEAD_TIME_S - self.park_elapsed)
            self.rate_eff = self._tapered_rate()
            slope = self._measured_slope()
            sgn = 1.0 if self.rate_c_per_s >= 0 else -1.0
            responding = (slope is not None and
                          sgn * slope >= START_SLOPE_FRACTION * abs(self.rate_c_per_s))
            if ((self.park_elapsed >= START_MIN_WAIT_S and responding)
                    or self.park_elapsed >= START_MAX_WAIT_S):
                self.parked = False
                self.start_delay_left = 0.0
                between = (self.last_temp < self.target_c if sgn > 0
                           else self.last_temp > self.target_c)
                if between:
                    self.t_setpoint = self.last_temp
                self.say(f"Stage responding after {self.park_elapsed:.0f} s -- "
                         f"setpoint released at {self.t_setpoint:.1f} C.")
            self.out_q.put(("leash", False))
            return
        if self.state == "RAMP":
            self.rate_eff = self._tapered_rate()
            proposed = self.t_setpoint + self.rate_eff * dt
            # +1 climbing, -1 descending; "ahead" means further along the ramp
            sgn = 1.0 if self.rate_c_per_s >= 0 else -1.0
            ahead = sgn * (self.last_temp - proposed)

            # Stage AHEAD of the plan by more than the leash.  That is a
            # disturbance, not something to brake against -- pulling the stage
            # back down to the old schedule wastes heat and winds the
            # integrator negative (18 Sep: 36 C ahead, then a 17-minute stall).
            # Move the plan up to meet the stage instead.
            if ahead > LEASH_C:
                proposed = self.last_temp - sgn * LEASH_C
                self.rebased = True

            reached = (proposed >= self.target_c - 0.05 if sgn > 0
                       else proposed <= self.target_c + 0.05)
            if reached:
                self.t_setpoint = self.target_c
                self.state = "SOAK"
                self.say(f"Target reached. Soaking at {self.target_c:.0f} C.")
                return

            # Stage BEHIND the plan by more than the leash: hold the setpoint.
            if -ahead > LEASH_BEHIND_C:
                self.leashed = True
                self.out_q.put(("leash", True))
                return
            self.out_q.put(("leash", False))
            self.t_setpoint = proposed

        elif self.state == "COOLDOWN":
            step = COOLDOWN_RATE_C_PER_MIN / 60.0 * dt
            # never let "cooldown" raise the setpoint: if we are already below
            # COOLDOWN_END_C, the floor is where we are, not the constant
            floor = min(COOLDOWN_END_C, self.t_setpoint)
            self.t_setpoint = max(floor, self.t_setpoint - step)

    def _park(self):
        self.parked = True
        self.park_elapsed = 0.0
        self.hist = []

    def _measured_slope(self):
        """Least-squares slope of the recent temperature history, C/s."""
        if len(self.hist) < 5 or self.hist[-1][0] - self.hist[0][0] < 0.7 * SLOPE_WINDOW_S:
            return None
        n = len(self.hist)
        tm = sum(t for t, _ in self.hist) / n
        ym = sum(y for _, y in self.hist) / n
        sxx = sum((t - tm) ** 2 for t, _ in self.hist)
        sxy = sum((t - tm) * (y - ym) for t, y in self.hist)
        return sxy / sxx if sxx > 0 else None

    def _rate_at(self, sp):
        """Planned setpoint rate when the setpoint is at sp."""
        remaining = abs(self.target_c - sp)
        d_taper = TAPER_DEAD_TIMES * abs(self.rate_c_per_s) * DEAD_TIME_S
        if d_taper <= 0.0 or remaining >= d_taper:
            return self.rate_c_per_s
        return self.rate_c_per_s * max(math.sqrt(remaining / d_taper), TAPER_FLOOR)

    def _tapered_rate(self):
        return self._rate_at(self.t_setpoint)

    def _future_setpoint(self, horizon_s):
        """Where the planned setpoint will be horizon_s from now, and its rate."""
        sp = self.t_setpoint
        if self.state != "RAMP" or sp is None:
            return sp, 0.0
        climbing = self.rate_c_per_s >= 0
        step = 3.0
        horizon_s -= self.start_delay_left      # setpoint is parked until then
        if horizon_s <= 0.0:
            return sp, 0.0
        for _ in range(int(horizon_s / step)):
            nxt = sp + self._rate_at(sp) * step
            if (climbing and nxt >= self.target_c) or (not climbing and nxt <= self.target_c):
                return self.target_c, 0.0
            sp = nxt
        return sp, self._rate_at(sp)

    @staticmethod
    def _honored_rate(distance_c, minutes, start_delay_s=0.0):
        """Rate that makes the ramp -- start delay and taper included -- take `minutes`."""
        t = minutes * 60.0
        moving = t - start_delay_s - TAPER_DEAD_TIMES * DEAD_TIME_S
        return distance_c / max(moving, 0.4 * t)

    def control_step(self, dt):
        error = self.t_setpoint - self.last_temp

        # Feedforward aims LOOKAHEAD_S ahead along the planned ramp, since that
        # is when power applied now actually lands on the sample.  Feedback
        # still compares against the setpoint now.
        if self.state == "RAMP":
            sp_ff, rate = self._future_setpoint(LOOKAHEAD_S)
        else:
            sp_ff, rate = self.t_setpoint, 0.0

        v_ff = feedforward_volts(
            sp_ff, rate, self.r_live, self.p["heat_capacity"]
        )
        kp_s, ki_s = scheduled_gains(self.t_setpoint, v_ff, self.p["heat_capacity"])
        kp = kp_s * self.p["kp"]
        ki = ki_s * self.p["ki"]
        self.kp_now, self.ki_now = kp, ki
        self.out_q.put(("gains", kp, ki))
        v_p = kp * error

        v_unclamped = v_ff + v_p + self.i_term
        v_target = max(0.0, min(V_MAX, v_unclamped))

        # slew limit, asymmetric: rising is limited, falling is nearly free.
        # A resistance that suddenly drops (a contact healing, 18 Sep) turns a
        # held voltage into a power spike; cutting back must not be slow.
        up = V_SLEW_UP_V_PER_MIN / 60.0 * dt
        down = V_SLEW_DOWN_V_PER_MIN / 60.0 * dt
        v_target = max(self.v_cmd - down, min(self.v_cmd + up, v_target))

        # Anti-windup.  Two ways the loop can lie to the integrator: the
        # command is clamped or slew-limited, or a leash is pinning the error.
        # While pinned, the error is an artifact of the limiter, so the
        # integrator may not ACCUMULATE from it -- but it must still be allowed
        # to UNWIND back toward zero.  Freezing it outright (the 16 Sep fix)
        # locked in -4 V of charge on 18 Sep and stalled the ramp for 17 min.
        saturated = abs(v_target - v_unclamped) > 1e-6
        if not saturated and self.state in ("RAMP", "SOAK"):
            step = ki * error * dt
            if self.leashed or self.rebased:
                if self.i_term * step < 0.0:            # moving toward zero
                    new = self.i_term + step
                    self.i_term = 0.0 if new * self.i_term < 0.0 else new
            else:
                self.i_term += step
            self.i_term = max(-I_TERM_CLAMP_V, min(I_TERM_CLAMP_V, self.i_term))

        self.v_cmd = v_target
        return self.psu.set_voltage(v_target)

    def check_psu(self, v_commanded):
        try:
            v_meas, i_meas, in_cc = self.psu.read()
            self.psu_fails = 0
        except Exception as exc:
            self.psu_fails += 1
            if self.psu_fails >= MAX_CONSECUTIVE_PSU_FAILS:
                self.fault(f"Supply unresponsive: {exc}")
            return None

        power = v_meas * i_meas

        # Track heater resistance, but only trust it when it is plausible.
        # On 18 Sep the circuit read 100-290 ohm for four minutes (a marginal
        # contact), feedforward trusted it and drove 19 V, and when the
        # contact healed that became a 16 W spike.  An implausible reading
        # means "something is wrong with the wiring", not "push harder".
        if self.psu.output_on and v_meas > 3.0 and self.last_temp is not None:
            r_model = resistance_ohm(self.last_temp)
            i_expected = v_meas / r_model
            if i_meas < i_expected / R_SANITY_FACTOR or i_meas > i_expected * R_SANITY_FACTOR:
                r_seen = v_meas / i_meas if i_meas > 0.005 else float("inf")
                self.r_live = None
                if not self.r_warned:
                    self.r_warned = True
                    self.say(
                        f"WARNING: heater reads {r_seen:.0f} ohm, expected ~{r_model:.0f}. "
                        f"Check heater connections. Using model resistance."
                    )
                    self.out_q.put(("r_alarm", r_seen, r_model))
            elif i_meas >= 0.10:
                self.r_live = v_meas / i_meas
                if self.r_warned:
                    self.r_warned = False
                    self.say(f"Heater resistance back to normal ({self.r_live:.1f} ohm).")

        if i_meas > I_TRIP_A:
            self.fault(f"Overcurrent: {i_meas:.2f} A readback")
        elif power > P_TRIP_W:
            self.fault(f"Overpower: {power:.1f} W readback")
        elif self.psu.output_on and abs(v_meas - v_commanded) > V_READBACK_TOLERANCE \
                and v_commanded > Bk1685b.V_FLOOR + 0.1:
            self.fault(
                f"Supply not following: commanded {v_commanded:.1f} V, "
                f"reading {v_meas:.2f} V"
            )
        elif in_cc and self.state in ("RAMP", "SOAK"):
            # current-limited means VOLT commands no longer do anything
            self.say("WARNING: supply is in constant-current mode.")

        return v_meas, i_meas, power, in_cc

    # ----------------------------------------------------------- the loop

    def run(self):
        self.last_temp = None
        rm = None
        try:
            rm = pyvisa.ResourceManager()
            self.dmm = rm.open_resource(DMM_ADDRESS)
            self.dmm.timeout = 10000
            idn = self.dmm.query("*IDN?").strip()
            self.dmm.write("*RST")
            self.dmm.write("*CLS")
            self.dmm.write("CONF:VOLT:DC 0.1,1E-6")
            self.dmm.write("VOLT:DC:NPLC 10")
            self.dmm.write("TRIG:SOUR IMM")
            self.say(f"DMM: {idn}")

            self.psu = Bk1685b(log=self.say)
            self.psu.configure()
            atexit.register(self.psu.shutdown)
            v0, i0, _ = self.psu.read()
            self.say(f"PSU connected. Output off, reading {v0:.2f} V / {i0:.2f} A.")

            self.state = "IDLE"
            self.out_q.put(("ready", None))
        except Exception as exc:
            self.fault(f"Connection failed: {exc}")
            return

        while not self.stop_flag.is_set():
            loop_start = time.time()
            self.drain_commands()

            # ---- read the thermocouple ----
            temp_c = float("nan")
            emf_mv = float("nan")
            volts_raw = float("nan")
            cj = self.p["cold_junction"]
            try:
                volts_raw = float(self.dmm.query("READ?"))
                emf_mv = volts_raw * 1000.0
                self.read_fails = 0
                try:
                    temp_c = measured_to_celsius(emf_mv, cj)
                    self.last_temp = temp_c
                except ValueError as exc:
                    if self.state in ("RAMP", "HOLD", "SOAK", "COOLDOWN"):
                        self.fault(f"Thermocouple out of range: {exc}")
            except Exception as exc:
                self.read_fails += 1
                self.say(f"Read error ({self.read_fails}): {exc}")
                if self.read_fails >= MAX_CONSECUTIVE_READ_FAILS:
                    if self.state in ("RAMP", "HOLD", "SOAK", "COOLDOWN"):
                        self.fault("Thermocouple lost for 3 consecutive reads.")

            # ---- control ----
            self.tick += 1
            psu_data = None
            v_commanded = self.v_cmd

            if self.state in ("RAMP", "HOLD", "SOAK", "COOLDOWN"):
                if time.time() > self.run_deadline:
                    self.fault(f"Run exceeded {MAX_RUN_HOURS} h limit.")
                elif not math.isnan(temp_c) and temp_c > self.target_c + T_OVERTEMP_MARGIN_C:
                    self.fault(
                        f"Overtemperature: {temp_c:.0f} C, "
                        f"{T_OVERTEMP_MARGIN_C:.0f} C above target."
                    )
                elif self.tick % CONTROL_EVERY_N == 0 and not math.isnan(temp_c):
                    now = time.time()
                    dt = now - self.last_control_t
                    self.last_control_t = now
                    self.advance_setpoint(dt)
                    v_commanded = self.control_step(dt)
                    psu_data = self.check_psu(v_commanded)

                    if self.state == "COOLDOWN":
                        # Once the command has sat at the supply floor for a
                        # while, passive cooling is faster than the requested
                        # descent and the heater is contributing nothing.  Turn
                        # it off and keep logging.
                        if v_commanded <= Bk1685b.V_FLOOR + 0.05:
                            self.floor_since = self.floor_since or now
                        else:
                            self.floor_since = None
                        at_floor = (self.floor_since is not None and
                                    now - self.floor_since >= COOLDOWN_OFF_AFTER_S)
                        done = (self.t_setpoint <= COOLDOWN_END_C and
                                temp_c <= COOLDOWN_END_C + 5)
                        if at_floor or done:
                            self.psu.shutdown()
                            self.state = "IDLE"
                            self.floor_since = None
                            self.out_q.put(("idle", None))
                            self.say(
                                "Cooldown complete. Output off." if done else
                                f"Heater no longer needed at {temp_c:.0f} C -- output off. "
                                f"Stage cooling on its own; still logging."
                            )

            self.out_q.put((
                "data", datetime.now(), volts_raw, emf_mv, cj, temp_c,
                self.t_setpoint, v_commanded,
                psu_data[0] if psu_data else float("nan"),
                psu_data[1] if psu_data else float("nan"),
                psu_data[2] if psu_data else float("nan"),
                self.state,
            ))

            remaining = READ_INTERVAL_S - (time.time() - loop_start)
            if remaining > 0:
                self.stop_flag.wait(remaining)

        if self.psu:
            self.psu.close()
        self.say("Worker stopped, output disabled.")


# ==================================================================== GUI


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("PBN Ramp Controller -- 34401A + BK1685B")
        self.geometry("1080x780")
        self.protocol("WM_DELETE_WINDOW", self.on_close)

        self.out_q = queue.Queue()
        self.cmd_q = queue.Queue()
        self.params = {
            "cold_junction": DEFAULT_COLD_JUNCTION_C,
            "kp": KP_DEFAULT,
            "ki": KI_DEFAULT,
            "heat_capacity": HEAT_CAPACITY_J_PER_K,
        }

        self.times, self.temps, self.setpoints = [], [], []
        self.v_cmds, self.powers = [], []
        self.csv_file = None
        self.csv_writer = None
        self.start_wall = None
        self.state = "CONNECTING"

        self._build_controls()
        self._build_plots()
        self.open_log()

        self.worker = RampWorker(self.out_q, self.cmd_q, self.params)
        self.worker.start()
        self.after(100, self.pump)

    # ------------------------------------------------------------ layout

    def _build_controls(self):
        bar = ttk.Frame(self, padding=8)
        bar.pack(side=tk.TOP, fill=tk.X)

        ttk.Label(bar, text="Target (C):").grid(row=0, column=0, sticky="e")
        self.e_target = ttk.Entry(bar, width=8)
        self.e_target.insert(0, "600")
        self.e_target.grid(row=0, column=1, padx=(2, 10))
        self.e_target.bind("<KeyRelease>", lambda _: self.update_preview())

        ttk.Label(bar, text="Ramp time (min):").grid(row=0, column=2, sticky="e")
        self.e_minutes = ttk.Entry(bar, width=8)
        self.e_minutes.insert(0, "45")
        self.e_minutes.grid(row=0, column=3, padx=(2, 10))
        self.e_minutes.bind("<KeyRelease>", lambda _: self.update_preview())

        ttk.Label(bar, text="Cold junction (C):").grid(row=0, column=4, sticky="e")
        self.e_cj = ttk.Entry(bar, width=8)
        self.e_cj.insert(0, str(DEFAULT_COLD_JUNCTION_C))
        self.e_cj.grid(row=0, column=5, padx=(2, 10))

        ttk.Label(bar, text="Kp \u00d7").grid(row=1, column=0, sticky="e", pady=(6, 0))
        self.e_kp = ttk.Entry(bar, width=8)
        self.e_kp.insert(0, str(KP_DEFAULT))
        self.e_kp.grid(row=1, column=1, padx=(2, 10), pady=(6, 0))

        ttk.Label(bar, text="Ki \u00d7").grid(row=1, column=2, sticky="e", pady=(6, 0))
        self.e_ki = ttk.Entry(bar, width=8)
        self.e_ki.insert(0, str(KI_DEFAULT))
        self.e_ki.grid(row=1, column=3, padx=(2, 10), pady=(6, 0))

        ttk.Label(bar, text="C (J/K):").grid(row=1, column=4, sticky="e", pady=(6, 0))
        self.e_cap = ttk.Entry(bar, width=8)
        self.e_cap.insert(0, str(HEAT_CAPACITY_J_PER_K))
        self.e_cap.grid(row=1, column=5, padx=(2, 10), pady=(6, 0))

        self.b_start = ttk.Button(bar, text="Start Ramp", command=self.on_start)
        self.b_start.grid(row=0, column=6, padx=4)
        self.b_apply = ttk.Button(bar, text="Apply Changes",
                                  command=self.on_apply, state="disabled")
        self.b_apply.grid(row=1, column=6, padx=4, pady=(6, 0))
        self.b_hold = ttk.Button(bar, text="HOLD", command=self.on_hold, state="disabled")
        self.b_hold.grid(row=0, column=7, padx=4)
        self.b_cool = ttk.Button(bar, text="Cooldown", command=self.on_cooldown,
                                 state="disabled")
        self.b_cool.grid(row=0, column=8, padx=4)
        self.b_stop = tk.Button(bar, text="EMERGENCY STOP", bg="#c0392b", fg="white",
                                font=("TkDefaultFont", 10, "bold"),
                                command=self.on_estop)
        self.b_stop.grid(row=0, column=9, rowspan=2, padx=(16, 4), ipadx=8, ipady=6)

        self.preview = ttk.Label(bar, text="", foreground="#555")
        self.preview.grid(row=2, column=0, columnspan=10, sticky="w", pady=(6, 0))

        self.readout = ttk.Label(self, text="Connecting...",
                                 font=("TkDefaultFont", 14, "bold"), padding=(8, 2))
        self.readout.pack(side=tk.TOP, anchor="w")
        self.status = ttk.Label(self, text="", foreground="#333", padding=(8, 0))
        self.status.pack(side=tk.TOP, anchor="w")
        self.gains_lbl = ttk.Label(self, text="", foreground="#777", padding=(8, 0))
        self.gains_lbl.pack(side=tk.TOP, anchor="w")

        self.update_preview()

    def _build_plots(self):
        fig = Figure(figsize=(10, 5.4), dpi=100)
        self.ax_t = fig.add_subplot(211)
        self.ax_v = fig.add_subplot(212, sharex=self.ax_t)
        self.ax_t.set_ylabel("Temperature (C)")
        self.ax_v.set_ylabel("Volts / Watts")
        self.ax_v.set_xlabel("Time")
        for ax in (self.ax_t, self.ax_v):
            ax.grid(alpha=0.3)
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
        fig.tight_layout()
        self.canvas = FigureCanvasTkAgg(fig, master=self)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # ----------------------------------------------------------- actions

    def update_preview(self):
        try:
            target = float(self.e_target.get())
            minutes = float(self.e_minutes.get())
        except ValueError:
            self.preview.config(text="")
            return
        base = self.temps[-1] if self.temps else 25.0
        rate = (target - base) / minutes if minutes > 0 else 0.0
        v_hold = feedforward_volts(target)
        try:
            cap = float(self.e_cap.get())
        except ValueError:
            cap = HEAT_CAPACITY_J_PER_K
        v_climb = feedforward_volts(base + 20.0, rate / 60.0, None, cap)
        head = V_MAX - v_hold
        limit = max_rate_c_per_min(target)
        note = (f"Hold {target:.0f} C needs about {v_hold:.1f} V "
                f"({head:.1f} V headroom); starting the climb needs about "
                f"{v_climb:.1f} V.  Requested {rate:.1f} C/min, "
                f"stage manages about {limit:.0f} C/min there.")
        self.preview.config(text=note)

    def on_start(self):
        try:
            target = float(self.e_target.get())
            minutes = float(self.e_minutes.get())
            self.params["cold_junction"] = float(self.e_cj.get())
            self.params["kp"] = float(self.e_kp.get())
            self.params["ki"] = float(self.e_ki.get())
            self.params["heat_capacity"] = float(self.e_cap.get())
        except ValueError:
            messagebox.showerror("Bad input", "Check the numeric fields.")
            return

        if not 0 < target <= T_TARGET_MAX_C:
            messagebox.showerror(
                "Target out of range",
                f"Target must be between 0 and {T_TARGET_MAX_C:.0f} C.\n\n"
                f"Above that the 40 V ceiling leaves the controller no headroom.",
            )
            return
        if minutes <= 0:
            messagebox.showerror("Bad input", "Ramp time must be positive.")
            return

        base = self.temps[-1] if self.temps else 25.0
        rate = (target - base) / minutes
        limit = max_rate_c_per_min(target)
        if rate > limit:
            ok = messagebox.askyesno(
                "Ramp probably too fast",
                f"You asked for {rate:.1f} C/min but the stage managed about "
                f"{limit:.0f} C/min near {target:.0f} C on 28 Aug.\n\n"
                f"The setpoint waits whenever it gets more than "
                f"{LEASH_BEHIND_C:.0f} C ahead of the stage, so the run will simply take longer than "
                f"{minutes:.0f} min.\n\nStart anyway?",
            )
            if not ok:
                return

        self.cmd_q.put(("start", target, minutes))
        self.b_start.config(state="disabled")
        self.b_apply.config(state="normal")
        self.b_hold.config(state="normal", text="HOLD")
        self.b_cool.config(state="normal")

    def on_apply(self):
        """Push edited fields into a run that is already going."""
        try:
            target = float(self.e_target.get())
            minutes = float(self.e_minutes.get())
            self.params["cold_junction"] = float(self.e_cj.get())
            self.params["kp"] = float(self.e_kp.get())
            self.params["ki"] = float(self.e_ki.get())
            self.params["heat_capacity"] = float(self.e_cap.get())
        except ValueError:
            messagebox.showerror("Bad input", "Check the numeric fields.")
            return
        if not 0 < target <= T_TARGET_MAX_C:
            messagebox.showerror(
                "Target out of range",
                f"Target must be between 0 and {T_TARGET_MAX_C:.0f} C.",
            )
            return
        if minutes <= 0:
            messagebox.showerror("Bad input", "Ramp time must be positive.")
            return
        self.cmd_q.put(("retarget", target, minutes))

    def on_hold(self):
        if self.state == "RAMP":
            self.cmd_q.put(("hold",))
            self.b_hold.config(text="RESUME")
        else:
            self.cmd_q.put(("resume",))
            self.b_hold.config(text="HOLD")

    def on_cooldown(self):
        self.cmd_q.put(("cooldown",))
        self.b_hold.config(state="disabled")
        self.b_apply.config(state="disabled")

    def on_estop(self):
        self.cmd_q.put(("estop",))

    # -------------------------------------------------------------- data

    def open_log(self):
        os.makedirs(SAVE_DIR, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        path = os.path.join(SAVE_DIR, f"ramp_log_{stamp}.csv")
        self.csv_file = open(path, "w", newline="")
        self.csv_writer = csv.writer(self.csv_file)
        # first six columns match the old tc_log format so existing viewers
        # keep working; epoch_s is there for aligning against RGA exports
        self.csv_writer.writerow([
            "iso_time", "elapsed_s", "volts", "emf_mv", "cold_junction_c", "temp_c",
            "epoch_s", "setpoint_c", "v_cmd", "v_meas", "i_meas",
            "power_w", "resistance_ohm", "state",
        ])
        self.status.config(text=f"Logging to {path}")

    def pump(self):
        try:
            while True:
                msg = self.out_q.get_nowait()
                kind = msg[0]

                if kind == "status":
                    warn = msg[1].startswith("WARNING")
                    self.status.config(text=msg[1],
                                       foreground="#c0392b" if warn else "#333")
                elif kind in ("ready", "idle"):
                    self.readout.config(text="Idle -- output off", foreground="black")
                    self.b_start.config(state="normal")
                    self.b_apply.config(state="disabled")
                    self.b_hold.config(state="disabled", text="HOLD")
                    self.b_cool.config(state="disabled")
                elif kind == "leash":
                    pass
                elif kind == "r_alarm":
                    self.bell()
                    self.status.config(foreground="#c0392b")
                elif kind == "gains":
                    self.gains_lbl.config(
                        text=f"Scheduled gains now:  Kp = {msg[1]:.3f} V/C,  "
                             f"Ki = {msg[2]:.2e} V/(C s)"
                    )
                elif kind == "fault":
                    self.state = "FAULT"
                    self.bell()
                    self.readout.config(text=f"FAULT -- {msg[1]}", foreground="#c0392b")
                    self.b_start.config(state="disabled")
                    self.b_apply.config(state="disabled")
                    self.b_hold.config(state="disabled")
                    self.b_cool.config(state="disabled")
                    messagebox.showerror(
                        "Fault -- heater power cut", f"{msg[1]}\n\nOutput disabled."
                    )
                elif kind == "data":
                    self.on_data(msg)
        except queue.Empty:
            pass

        self.redraw()
        self.after(250, self.pump)

    def on_data(self, msg):
        (_, when, volts_raw, emf_mv, cj, temp_c,
         setpoint, v_cmd, v_meas, i_meas, power, state) = msg

        self.state = state
        if self.start_wall is None:
            self.start_wall = when
        elapsed = (when - self.start_wall).total_seconds()

        self.times.append(when)
        self.temps.append(temp_c)
        self.setpoints.append(setpoint if setpoint is not None else float("nan"))
        self.v_cmds.append(v_cmd)
        self.powers.append(power)

        resistance = (v_meas / i_meas) if (i_meas and i_meas > 0.02) else float("nan")

        if self.csv_writer:
            self.csv_writer.writerow([
                when.isoformat(), f"{elapsed:.3f}", f"{volts_raw:.9f}",
                f"{emf_mv:.6f}", cj, f"{temp_c:.3f}", f"{when.timestamp():.3f}",
                "" if setpoint is None else f"{setpoint:.3f}",
                f"{v_cmd:.1f}", f"{v_meas:.2f}", f"{i_meas:.2f}",
                f"{power:.2f}", f"{resistance:.2f}", state,
            ])
            self.csv_file.flush()

        if state != "FAULT":
            sp = f" | SP {setpoint:.1f} C" if setpoint is not None else ""
            pw = f" | {power:.0f} W" if not math.isnan(power) else ""
            self.readout.config(
                text=f"{state} -- {temp_c:.1f} C{sp} | {v_cmd:.1f} V{pw}",
                foreground="black",
            )

    def redraw(self):
        if len(self.times) < 2:
            return
        step = max(1, len(self.times) // 2000)
        t = self.times[::step]

        self.ax_t.clear()
        self.ax_t.plot(t, self.temps[::step], color="#c0392b", lw=1.2, label="measured")
        self.ax_t.plot(t, self.setpoints[::step], color="#2980b9", lw=1.0,
                       ls="--", label="setpoint")
        self.ax_t.set_ylabel("Temperature (C)")
        self.ax_t.legend(loc="upper left", fontsize=8)
        self.ax_t.grid(alpha=0.3)

        self.ax_v.clear()
        self.ax_v.plot(t, self.v_cmds[::step], color="#27ae60", lw=1.2, label="V cmd")
        self.ax_v.plot(t, self.powers[::step], color="#8e44ad", lw=1.0, label="W")
        self.ax_v.axhline(V_MAX, color="#c0392b", lw=0.8, ls=":")
        self.ax_v.axhline(P_MAX_W, color="#8e44ad", lw=0.8, ls=":")
        self.ax_v.set_ylabel("Volts / Watts")
        self.ax_v.set_xlabel("Time")
        self.ax_v.legend(loc="upper left", fontsize=8)
        self.ax_v.grid(alpha=0.3)
        self.ax_v.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))

        self.canvas.draw_idle()

    # ------------------------------------------------------------- close

    def on_close(self):
        if self.state in ("RAMP", "HOLD", "SOAK", "COOLDOWN"):
            if not messagebox.askyesno(
                "Heater is live",
                "The heater is powered. Closing will cut power immediately.\n\nClose?",
            ):
                return
        self.cmd_q.put(("quit",))
        self.worker.stop_flag.set()
        self.worker.join(timeout=5.0)
        if self.csv_file:
            self.csv_file.close()
        self.destroy()


# ====================================================== notebook entry points

_APP = None


def emergency_off(address=PSU_ADDRESS):
    """
    Cut heater power no matter what else is running.  Safe to call at any
    time, including after a kernel restart has orphaned the supply.  Tries
    the live object first, then falls back to opening the port fresh.
    """
    global _LIVE_PSU
    if _LIVE_PSU is not None:
        try:
            _LIVE_PSU.shutdown()
            print("Output disabled via the running controller.")
            return
        except Exception as exc:
            print(f"Live handle failed ({exc}); opening the port directly...")

    rm = pyvisa.ResourceManager()
    dev = rm.open_resource(address)
    dev.baud_rate = 9600
    dev.write_termination = "\r"
    dev.read_termination = "\r"
    dev.timeout = 2000
    for text in ("VOLT010", "SOUT1"):
        dev.write(text)
        try:
            dev.read()
        except Exception:
            pass
    try:
        dev.write("GETD")
        print("GETD after shutdown:", dev.read())
    except Exception:
        pass
    dev.close()
    print("Output disabled.")


def launch():
    """
    Open the controller window.  Blocks this cell until the window closes --
    that is normal, the cell will show as running the whole time.

    Close the window to stop.  If you interrupt the cell instead, the finally
    block still cuts power.  If you RESTART THE KERNEL while the heater is
    live, nothing here runs -- use emergency_off() in a fresh cell.
    """
    global _APP
    if _APP is not None:
        try:
            _APP.on_close()
        except Exception:
            pass
        _APP = None

    _APP = App()
    try:
        _APP.mainloop()
    finally:
        try:
            _APP.worker.stop_flag.set()
            if _APP.worker.psu:
                _APP.worker.psu.shutdown()
        except Exception:
            pass
        emergency_off()


if __name__ == "__main__":
    launch()

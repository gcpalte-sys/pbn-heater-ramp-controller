# PBN Heater Ramp Controller

Closed-loop temperature ramp controller for a pyrolytic boron nitride (PBN) heating stage in an ultra-high vacuum system. It reads a type-K thermocouple through an HP/Agilent 34401A multimeter over GPIB and drives a BK Precision 1685B DC supply over serial. It follows a linear ramp to a target temperature, soaks there, then runs a controlled cooldown, with a live Tkinter GUI, real-time plots, and CSV logging.

<!-- Add a screenshot of the GUI during a run here: -->
<!-- ![Controller GUI during a ramp](screenshot.png) -->

## Features

- **Feedforward + gain-scheduled PI control**, built from a physical model of the stage fit to measured heating runs
- **Look-ahead feedforward** that compensates for the 40–80 s thermal dead time, cutting overshoot near the target
- **Honored ramp time**: the time you enter is the time the ramp takes, including the start delay and final taper
- **Adaptive start**: the setpoint stays parked until the stage is measurably responding, then rebases onto the real temperature
- **Live retargeting**: change the target or remaining time mid-run without restarting or bumping the output
- **Layered hardware safety**: voltage, current, and power limits enforced against the supply's own readback, plus fault trips and an emergency stop
- **ITS-90 type-K conversion** with cold-junction compensation
- **CSV logs** with epoch timestamps for aligning against RGA and other instrument data

## Control approach

Steady-state heat loss is radiation-dominated, and the pyrolytic graphite element's resistance falls by nearly half as it heats up. The controller therefore computes feedforward in **power** and converts to voltage at the end:

$$P_{ff} = P_{loss}(T) + C\,\frac{dT}{dt}, \qquad V_{ff} = \sqrt{P_{ff}\,R(T)}$$

where $P_{loss}(T)$ is interpolated from measured steady-state points, $C$ is the stage heat capacity, and $R(T)$ is a fitted resistance curve (or the live value from the supply readback, when it passes a sanity check).

The PI loop only trims the residual. Because the stage's time constant spans roughly 46–900 s across the working range, the gains are scheduled with SIMC tuning, linearized about the current operating point:

$$K_p = \frac{C\,R}{2V(\tau_c + L)}, \qquad K_i = \frac{K_p}{\min\!\big(\tau,\ 4(\tau_c + L)\big)}$$

where $L$ is the dead time, $\tau_c$ is the target closed-loop time constant, and $\tau = C / (dP_{loss}/dT)$. There is no derivative term, since with this much dead time it would mostly amplify noise.

The integrator has anti-windup that blocks accumulation while the command is saturated or the setpoint is leashed, but still lets the term unwind back toward zero.

## Safety

| Protection | Behavior |
|---|---|
| Voltage / current / power limits | 40 V, 3.00 A, 120 W (element nameplate), power checked against `GETD` readback every loop |
| Overcurrent / overpower trip | Output cut above 3.10 A or 126 W measured |
| Overtemperature trip | Output cut 30 °C above target |
| Supply tracking | Fault if measured voltage disagrees with the command by more than 1 V |
| Heater wiring check | Implausible resistance readings are rejected and flagged instead of trusted |
| Asymmetric slew limit | Voltage rises slowly but can drop fast |
| Sensor loss | Fault after 3 consecutive failed thermocouple reads |
| Run timeout | Fault after 8 hours |
| Shutdown paths | E-stop button, window close, `atexit` hook, and `emergency_off()` for orphaned sessions |

## Hardware

- HP/Agilent 34401A DMM (GPIB), reading a type-K thermocouple
- BK Precision 1685B DC power supply (USB/serial)
- PBN heating stage with a pyrolytic graphite element

## Installation

```bash
pip install pyvisa matplotlib
```

You also need a VISA backend (NI-VISA or Keysight IO Libraries) with a working GPIB interface.

## Usage

Set the instrument addresses at the top of `ramp_controller.py`:

```python
DMM_ADDRESS = "GPIB0::6::INSTR"
PSU_ADDRESS = "ASRL7::INSTR"
```

Then run:

```bash
python ramp_controller.py
```

Enter a target temperature and ramp time, and press **Start Ramp**. The preview line estimates the voltage needed and warns if the requested rate is faster than the stage can manage.

From a Jupyter notebook:

```python
from ramp_controller import launch, emergency_off
launch()          # opens the GUI; the cell blocks until the window closes
```

If the kernel is restarted while the heater is live, run `emergency_off()` in a fresh cell to cut power.

## Log format

Logs are saved to `~/tc_logs/ramp_log_YYYYMMDD_HHMMSS.csv` with the columns:

`iso_time, elapsed_s, volts, emf_mv, cold_junction_c, temp_c, epoch_s, setpoint_c, v_cmd, v_meas, i_meas, power_w, resistance_ohm, state`

## Note

The model constants (loss power, resistance curve, heat capacity, dead time) are fit to one specific stage. On different hardware they need to be re-measured; the PI will cover small errors but not a different system.

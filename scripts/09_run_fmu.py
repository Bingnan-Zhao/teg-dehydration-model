"""
09_run_fmu.py  -  run the exported FMU (modelica/TEGUnit.fmu) with the same step test
and compare it with the Python twin (results/08_twin_steptest.csv).

Needs fmpy (installed by setup_env.bat). A Dymola FMU exported without the
"Binary Model Export" licence option needs a Dymola licence on the machine
that runs it.

Outputs: results/09_fmu_steptest.csv, results/09_fmu_vs_twin.png, results/09_summary.txt
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from fmpy import read_model_description, simulate_fmu

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FMU = os.path.join(PROJECT, "modelica", "TEGUnit.fmu")
OUTPUTS = ["dryGasWater", "T_reboiler", "leanTEG", "surgeLevel", "richTEG", "boilUp"]

md = read_model_description(FMU)
print("FMU:", md.modelName, "| FMI", md.fmiVersion,
      "| inputs:", [v.name for v in md.modelVariables if v.causality == "input"])

# piecewise-constant input table (same steps as StepTest)
t_end = 26 * 3600
times = [0, 2 * 3600, 8 * 3600, 14 * 3600, 20 * 3600]
rows = []
for k, t in enumerate(times):
    q = 1.0e5 - (1.5e4 if t >= 2 * 3600 else 0)
    m = 0.473 - (0.15 if t >= 8 * 3600 else 0)
    f = 1.42 + (0.3 if t >= 14 * 3600 else 0)
    tg = 40 + (5 if t >= 20 * 3600 else 0)
    if k > 0:                                   # duplicate time stamp -> step in the input table
        prev = rows[-1]
        rows.append((t,) + prev[1:])
    rows.append((t, q, m, f, tg))
rows.append((t_end,) + rows[-1][1:])
signals = np.array(rows, dtype=[("time", np.float64), ("Q_reb", np.float64), ("m_circ", np.float64),
                                ("F_gas", np.float64), ("T_gas", np.float64)])

fmi_type = "CoSimulation" if md.coSimulation is not None else "ModelExchange"
res = simulate_fmu(FMU, fmi_type=fmi_type, start_time=0, stop_time=t_end, output_interval=60,
                   input=signals, output=OUTPUTS)
fmu = pd.DataFrame({"time_h": res["time"] / 3600, **{k: res[k] for k in OUTPUTS}})
fmu.to_csv(os.path.join(PROJECT, "results", "09_fmu_steptest.csv"), index=False)

twin_path = os.path.join(PROJECT, "results", "08_twin_steptest.csv")
lines = [f"FMU run ({fmi_type}): {len(fmu)} points over {t_end / 3600:.0f} h"]
fig, ax = plt.subplots(4, 1, figsize=(9, 9), sharex=True)
cols = ["dryGasWater", "T_reboiler", "leanTEG", "surgeLevel"]
twin = pd.read_csv(twin_path) if os.path.exists(twin_path) else None
for a, col in zip(ax, cols):
    a.plot(fmu.time_h, fmu[col], label="FMU (Dymola)")
    if twin is not None:
        a.plot(twin.time_h, twin[col], "--", label="Python twin")
        diff = np.interp(fmu.time_h, twin.time_h, twin[col]) - fmu[col]
        lines.append(f"{col:12s} max |FMU - twin| = {np.abs(diff).max():.3g}")
    a.set_ylabel(col)
    a.grid(alpha=0.3)
ax[0].legend(fontsize=8)
ax[-1].set_xlabel("Time (h)")
fig.tight_layout()
fig.savefig(os.path.join(PROJECT, "results", "09_fmu_vs_twin.png"), dpi=150)
text = "\n".join(lines)
print(text)
open(os.path.join(PROJECT, "results", "09_summary.txt"), "w").write(text + "\n")

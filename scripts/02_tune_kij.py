"""
02_tune_kij.py  -  Step 2 of the TEG dehydration case.

DWSIM's default Peng-Robinson methane-water kij (0.5) under-predicts the water
content of the raw gas by about 14 %. This script fits a single methane-water
kij so that the saturated water content matches the Bukacek correlation over
20-50 C and 3-10 MPa, then saves it to results/02_kij.json. All later scripts
read that file through teg_common.py.

Outputs: 02_kij.json, 02_kij_fit.csv, 02_kij_fit.png, 02_summary.txt
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

import teg_common as tc

GRID = [("T", t, tc.P_BASE_KPA) for t in (20, 25, 30, 35, 40, 45, 50)] + \
       [("P", tc.T_BASE_C, p) for p in (3000, 4000, 5000, 6000, 7000, 8000, 10000)]
REF = np.array([tc.bukacek_mg_sm3(t, p) for _, t, p in GRID])


def predict(kij):
    return np.array([tc.saturated_water(t, p, kij) for _, t, p in GRID])


def rms_rel(kij):
    return float(np.sqrt(np.mean((predict(kij) / REF - 1) ** 2)))


print("Fitting methane-water kij (about 1-2 minutes)...")
fit = minimize_scalar(rms_rel, bounds=(0.0, 0.6), method="bounded", options={"xatol": 1e-3})
kij = round(float(fit.x), 4)

default, tuned = predict(0.5), predict(kij)
df = pd.DataFrame(GRID, columns=["sweep", "T_C", "P_kPa"])
df["bukacek"], df["pr_default"], df["pr_tuned"] = REF, default, tuned
df["dev_default_pct"] = (df.pr_default / df.bukacek - 1) * 100
df["dev_tuned_pct"] = (df.pr_tuned / df.bukacek - 1) * 100

out = os.path.join(tc.PROJECT, "results")
df.round(2).to_csv(os.path.join(out, "02_kij_fit.csv"), index=False)
json.dump({"kij_methane_water": kij, "package": "Peng-Robinson",
           "fitted_to": "Bukacek (1955), 20-50 C at 7000 kPa and 3-10 MPa at 40 C",
           "rms_rel_dev_pct": round(fit.fun * 100, 2)},
          open(os.path.join(out, "02_kij.json"), "w"), indent=2)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
for ax, sweep, x, xlabel in ((axes[0], "T", "T_C", "Temperature (°C), 7000 kPa abs"),
                             (axes[1], "P", "P_kPa", "Pressure (kPa abs), 40 °C")):
    d = df[df.sweep == sweep]
    ax.plot(d[x], d.bukacek, "k--", lw=1.5, label="Bukacek correlation")
    ax.plot(d[x], d.pr_default, "o-", color="tab:gray", ms=4, label="PR, DWSIM default kij = 0.5")
    ax.plot(d[x], d.pr_tuned, "o-", color="tab:blue", ms=4, label=f"PR, tuned kij = {kij}")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Saturated water content (mg/Sm³)")
    ax.grid(alpha=0.3)
axes[0].legend(fontsize=8)
fig.suptitle(f"Tuning one kij cuts the error from {df.dev_default_pct.abs().mean():.0f} % "
             f"to {df.dev_tuned_pct.abs().mean():.1f} %")
fig.tight_layout()
fig.savefig(os.path.join(out, "02_kij_fit.png"), dpi=150)

text = "\n".join([
    f"Tuned Peng-Robinson methane-water kij = {kij} (DWSIM default 0.5)",
    f"Mean |deviation| from Bukacek: default {df.dev_default_pct.abs().mean():.1f} %, "
    f"tuned {df.dev_tuned_pct.abs().mean():.1f} % (max {df.dev_tuned_pct.abs().max():.1f} %)",
    f"Base case 40 C, 7000 kPa: Bukacek {tc.bukacek_mg_sm3(tc.T_BASE_C, tc.P_BASE_KPA):.0f}, "
    f"tuned PR {tc.saturated_water(tc.T_BASE_C, tc.P_BASE_KPA, kij):.0f} mg/Sm3",
])
print(text)
open(os.path.join(out, "02_summary.txt"), "w").write(text + "\n")

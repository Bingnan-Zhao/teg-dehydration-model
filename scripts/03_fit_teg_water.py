"""
03_fit_teg_water.py  -  Step 3 of the TEG dehydration case.

With DWSIM's default Peng-Robinson parameters, TEG barely holds water: gas in
equilibrium with 99 wt% TEG at 45 C keeps ~900 mg/Sm3 of water. This script
fits the TEG-water kij to open experimental VLE data (85 kPa, NIST ThermoML,
data/teg_water_vle_85kPa.csv) and checks what the tuned model predicts at
contactor conditions.

The fit weights the TEG-rich points (x_water <= 0.41), the region that matters
for both the contactor (~0.08) and the regenerator reboiler (~0.05-0.2).

Outputs: 03_kij_teg.json, 03_teg_water_fit.csv, 03_teg_water_fit.png,
03_summary.txt
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
from teg_common import Q

data = pd.read_csv(os.path.join(tc.PROJECT, "data", "teg_water_vle_85kPa.csv"), comment="#")
FIT = data[data.x_water <= 0.41]


def bubble_T(x_w, kij, p_kpa=85.0):
    fs = tc.new_flowsheet("bp", ["Water", tc.TEG], {(tc.TEG, "Water"): kij})
    s = (fs.AddMaterialStream("s").WithPressure(Q.Pascal(p_kpa * 1000)).WithVaporFraction(0.0)
         .WithMolarFlow(Q.MolPerSecond(1.0)).WithComposition(tc.composition({"Water": x_w, tc.TEG: 1 - x_w})))
    fs.Solve()
    return s.TemperatureK


def rms_T(kij):
    return float(np.sqrt(np.mean([(bubble_T(x, kij) - t) ** 2 for x, t in zip(FIT.x_water, FIT.T_K)])))


def contactor_equilibrium(kij_teg, t_c=45.0, p_kpa=7000.0, teg_wt=99.0):
    """Water content of gas in equilibrium with lean TEG (mg/Sm3)."""
    kijs = {(tc.TEG, "Water"): kij_teg}
    if tc.KIJ_CH4_H2O is not None:
        kijs[("Methane", "Water")] = tc.KIJ_CH4_H2O
    fs = tc.new_flowsheet("eq", list(tc.GAS) + ["Water", tc.TEG], kijs)
    mT, mW = teg_wt / 150.17, (100 - teg_wt) / 18.015
    xw = mW / (mT + mW)
    fr = {k: v * 0.5 for k, v in tc.GAS.items()}           # half gas, half lean TEG
    fr[tc.TEG], fr["Water"] = 0.5 * (1 - xw), 0.5 * xw
    f = (fs.AddMaterialStream("f").At(Q.Kelvin(t_c + 273.15), Q.Pascal(p_kpa * 1000))
         .WithMolarFlow(Q.MolPerSecond(100.0)).WithComposition(tc.composition(fr)))
    v, l = fs.AddMaterialStream("v"), fs.AddMaterialStream("l")
    fs.AddSeparator("S").ConnectFeed(f, 0).ConnectProduct(v, 0).ConnectProduct(l, 1)
    fs.Solve()
    return tc.y_to_mg_sm3(v.OverallMoleFraction("Water"))


print("Fitting TEG-water kij to ThermoML data (about 1 minute)...")
fit = minimize_scalar(rms_T, bounds=(-0.3, 0.1), method="bounded", options={"xatol": 1e-3})
kij = round(float(fit.x), 4)

data["T_default"] = [bubble_T(x, 0.0) for x in data.x_water]
data["T_tuned"] = [bubble_T(x, kij) for x in data.x_water]
out = os.path.join(tc.PROJECT, "results")
data.round(3).to_csv(os.path.join(out, "03_teg_water_fit.csv"), index=False)
json.dump({"kij_teg_water": kij, "package": "Peng-Robinson",
           "fitted_to": "Water + TEG bubble temperatures at 85 kPa, x_water <= 0.41, "
                        "J. Chem. Eng. Data 2009, 54, 876 (doi:10.1021/je800675u), via NIST ThermoML",
           "rms_T_K": round(fit.fun, 2)},
          open(os.path.join(out, "03_kij_teg.json"), "w"), indent=2)

eq_default = contactor_equilibrium(0.0)
eq_tuned = contactor_equilibrium(kij)

fig, ax = plt.subplots(figsize=(6.5, 4.5))
ax.plot(data.x_water, data.T_K - 273.15, "ks", ms=5, label="Experiment (ThermoML, 85 kPa)")
xs = np.linspace(0.05, 0.95, 19)
ax.plot(xs, [bubble_T(x, 0.0) - 273.15 for x in xs], "-", color="tab:gray", label="PR, default kij = 0")
ax.plot(xs, [bubble_T(x, kij) - 273.15 for x in xs], "-", color="tab:blue", label=f"PR, tuned kij = {kij}")
ax.axvspan(0.0, 0.41, color="tab:blue", alpha=0.06, label="Fitted region (TEG-rich)")
ax.set_xlabel("Water mole fraction in liquid")
ax.set_ylabel("Bubble temperature at 85 kPa (°C)")
ax.set_title("Default PR makes water far too volatile over TEG")
ax.grid(alpha=0.3)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(out, "03_teg_water_fit.png"), dpi=150)

fit_pts = data.x_water <= 0.41
text = "\n".join([
    f"Tuned Peng-Robinson TEG-water kij = {kij} (DWSIM default 0)",
    f"Bubble T error, TEG-rich points: default {np.sqrt(np.mean((data.T_default - data.T_K)[fit_pts] ** 2)):.1f} K RMS, "
    f"tuned {fit.fun:.1f} K RMS",
    f"Bubble T error, all points: tuned {np.sqrt(np.mean((data.T_tuned - data.T_K) ** 2)):.1f} K RMS",
    f"Gas in equilibrium with 99 wt% TEG at 45 C, 7000 kPa: default {eq_default:.0f} mg/Sm3, "
    f"tuned {eq_tuned:.0f} mg/Sm3 (spec {tc.SPEC_MG_SM3:.0f})",
    "Note: fitted at 98-148 C; contactor (40-45 C) is an extrapolation, to be checked against "
    "equilibrium dew point data (Bahadori & Vuthaluru 2009).",
])
print(text)
open(os.path.join(out, "03_summary.txt"), "w").write(text + "\n")

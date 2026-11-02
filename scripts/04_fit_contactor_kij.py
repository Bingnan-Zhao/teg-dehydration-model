"""
04_fit_contactor_kij.py  -  Step 4: calibrate TEG-water for contactor conditions.

Step 3 fitted the TEG-water kij to regenerator-range data (98-148 C). Here the
same kij is checked and refitted against equilibrium water dew points of gas
over TEG solutions at 10-60 C (data/teg_dewpoint_points.csv, from GPSA /
Parrish 1986 / Herskowitz & Gottlieb 1984 as reported by Bahadori & Vuthaluru 2009).

Model dew point: gas at the contactor temperature in equilibrium with TEG of the
given strength -> water content y -> temperature at which water-saturated gas
holds the same water (metastable liquid water, tuned PR).

Outputs: 04_kij_contactor.json, 04_dewpoint_fit.csv, 04_dewpoint_fit.png, 04_summary.txt
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

P_KPA = 7000.0
pts = pd.read_csv(os.path.join(tc.PROJECT, "data", "teg_dewpoint_points.csv"), comment="#")

# saturation curve of the gas (tuned PR), used to turn water content into a dew point
TD_GRID = np.arange(200.0, 336.0, 5.0)
SAT = np.log([tc.saturated_water(t - 273.15, P_KPA, tc.KIJ_CH4_H2O) for t in TD_GRID])


def dew_point_K(mg_sm3):
    return float(np.interp(np.log(mg_sm3), SAT, TD_GRID))


def equilibrium_water(kij_teg, t_k, teg_wt):
    kijs = {(tc.TEG, "Water"): kij_teg, ("Methane", "Water"): tc.KIJ_CH4_H2O}
    fs = tc.new_flowsheet("eq", list(tc.GAS) + ["Water", tc.TEG], kijs)
    mT, mW = teg_wt / 150.17, (100 - teg_wt) / tc.MW_W
    xw = mW / (mT + mW)
    fr = {k: v * 0.5 for k, v in tc.GAS.items()}
    fr[tc.TEG], fr["Water"] = 0.5 * (1 - xw), 0.5 * xw
    f = (fs.AddMaterialStream("f").At(Q.Kelvin(t_k), Q.Pascal(P_KPA * 1000))
         .WithMolarFlow(Q.MolPerSecond(100.0)).WithComposition(tc.composition(fr)))
    v, l, l2 = fs.AddMaterialStream("v"), fs.AddMaterialStream("l"), fs.AddMaterialStream("l2")
    fs.AddSeparator("S").ConnectFeed(f, 0).ConnectProduct(v, 0).ConnectProduct(l, 1).ConnectProduct(l2, 2)
    fs.Solve()
    return tc.y_to_mg_sm3(v.OverallMoleFraction("Water"))


def model_dew(kij):
    return np.array([dew_point_K(equilibrium_water(kij, t, w))
                     for w, t in zip(pts.teg_wt_pct, pts.T_contactor_K)])


def rms(kij):
    return float(np.sqrt(np.mean((model_dew(kij) - pts.T_dew_K) ** 2)))


print("Fitting TEG-water kij to contactor dew point data (about 1 minute)...")
fit = minimize_scalar(rms, bounds=(-0.5, 0.05), method="bounded", options={"xatol": 1e-3})
kij_c = round(float(fit.x), 4)
kij_regen = tc.KIJ_TEG_H2O

pts["Td_default"] = model_dew(0.0)
pts["Td_regen_kij"] = model_dew(kij_regen)
pts["Td_contactor_kij"] = model_dew(kij_c)
out = os.path.join(tc.PROJECT, "results")
pts.round(2).to_csv(os.path.join(out, "04_dewpoint_fit.csv"), index=False)
json.dump({"kij_teg_water_contactor": kij_c, "package": "Peng-Robinson",
           "fitted_to": "Equilibrium water dew points over TEG, 10-60 C, 90-99.9 wt% "
                        "(Bahadori & Vuthaluru 2009, Table 2 and example)",
           "rms_dew_K": round(fit.fun, 2)},
          open(os.path.join(out, "04_kij_contactor.json"), "w"), indent=2)

fig, ax = plt.subplots(figsize=(6, 5))
lim = [225, 300]
ax.plot([l - 273.15 for l in lim], [l - 273.15 for l in lim], "k:", lw=1, label="Perfect agreement")
ax.plot(pts.T_dew_K - 273.15, pts.Td_default - 273.15, "s", color="tab:gray", label="PR default kij = 0")
ax.plot(pts.T_dew_K - 273.15, pts.Td_regen_kij - 273.15, "^", color="tab:orange",
        label=f"PR kij = {kij_regen} (regenerator fit)")
ax.plot(pts.T_dew_K - 273.15, pts.Td_contactor_kij - 273.15, "o", color="tab:blue",
        label=f"PR kij = {kij_c} (contactor fit)")
ax.set_xlabel("Reported equilibrium dew point (°C)")
ax.set_ylabel("Model equilibrium dew point (°C)")
ax.set_title("Contactor-range calibration of TEG-water")
ax.grid(alpha=0.3)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(out, "04_dewpoint_fit.png"), dpi=150)


def r(col):
    return np.sqrt(np.mean((pts[col] - pts.T_dew_K) ** 2))


text = "\n".join([
    f"Contactor-range TEG-water kij = {kij_c} (regenerator fit {kij_regen}, DWSIM default 0)",
    f"RMS dew point error: default {r('Td_default'):.1f} K, regenerator kij {r('Td_regen_kij'):.1f} K, "
    f"contactor kij {r('Td_contactor_kij'):.1f} K",
    "Point by point (reported vs contactor-kij model, C):",
] + [f"  {w:6.2f} wt% at {t - 273.15:5.1f} C: {a - 273.15:6.1f} vs {b - 273.15:6.1f}"
     for w, t, a, b in zip(pts.teg_wt_pct, pts.T_contactor_K, pts.T_dew_K, pts.Td_contactor_kij)])
print(text)
open(os.path.join(out, "04_summary.txt"), "w").write(text + "\n")

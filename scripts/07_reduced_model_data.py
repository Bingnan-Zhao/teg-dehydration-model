"""
07_reduced_model_data.py  -  Step 7: data for the reduced-order dynamic model (Modelica/FMU).

The Modelica model cannot call DWSIM, so two relations are taken from the
calibrated DWSIM model and fitted to simple formulas:
  1. Bubble temperature of TEG-water at the still pressure vs water mass fraction
     (regenerator kij)  ->  T_bub = c0 + c1*ln(w) + c2*ln(w)^2
  2. Contactor performance map: dry gas water vs lean TEG purity, TEG/water
     ratio and gas temperature (contactor kij, 4 stages)
     ->  ln(dry) = quadratic polynomial in (lean wt%, ln ratio, T_gas)

Outputs: data/reduced_model_coeffs.json, results/07_bubble_fit.png,
results/07_contactor_map.csv, results/07_summary.txt   (about 3-5 minutes)
"""
import itertools
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import teg_common as tc
from teg_common import Q

P_STILL = 101.3
P_CONTACTOR = 7000.0
GAS_MOLS = 693.0
MW_TEG = 150.17
COMPS = list(tc.GAS) + ["Water", tc.TEG]
out = os.path.join(tc.PROJECT, "results")


def kpa(p):
    return Q.Pascal(p * 1000.0)


def mole_from_wt(wt_teg):
    mT, mW = wt_teg / MW_TEG, (1 - wt_teg) / tc.MW_W
    return {tc.TEG: mT / (mT + mW), "Water": mW / (mT + mW)}


# ---- 1. bubble curve ----------------------------------------------------------
def bubble_T(w_water):
    fs = tc.new_flowsheet("bp", ["Water", tc.TEG], tc.tuned_kijs("regenerator"))
    s = (fs.AddMaterialStream("s").WithPressure(kpa(P_STILL)).WithVaporFraction(0.0)
         .WithMolarFlow(Q.MolPerSecond(1.0)).WithComposition(tc.composition(mole_from_wt(1 - w_water))))
    fs.Solve()
    return s.TemperatureK - 273.15


W = np.array([0.003, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.04, 0.06, 0.08, 0.10])
TB = np.array([bubble_T(w) for w in W])
cb = np.polyfit(np.log(W), TB, 2)[::-1]           # c0, c1, c2
TB_fit = cb[0] + cb[1] * np.log(W) + cb[2] * np.log(W) ** 2

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(W * 100, TB, "ks", label="DWSIM (tuned PR, 101.3 kPa)")
wf = np.linspace(W.min(), W.max(), 100)
ax.plot(wf * 100, cb[0] + cb[1] * np.log(wf) + cb[2] * np.log(wf) ** 2, "-", label="fit used in Modelica")
ax.set_xlabel("Water in TEG (wt%)")
ax.set_ylabel("Bubble temperature (°C)")
ax.grid(alpha=0.3)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(out, "07_bubble_fit.png"), dpi=150)


# ---- 2. contactor map -----------------------------------------------------------
def contactor(lean_wt, ratio, t_gas):
    """Dry gas water (mg/Sm3) for lean TEG of lean_wt (fraction), TEG mass = ratio x water in wet gas."""
    fs = tc.new_flowsheet("contactor", COMPS, tc.tuned_kijs("contactor"))
    raw = {k: v * 0.98 for k, v in tc.GAS.items()}
    raw["Water"] = 0.02
    rg = (fs.AddMaterialStream("raw").At(Q.Kelvin(t_gas + 273.15), kpa(P_CONTACTOR))
          .WithMolarFlow(Q.MolPerSecond(GAS_MOLS / 0.98)).WithComposition(tc.composition(raw)))
    wet, fw = fs.AddMaterialStream("wet"), fs.AddMaterialStream("fw")
    fs.AddSeparator("V-100").ConnectFeed(rg, 0).ConnectProduct(wet, 0).ConnectProduct(fw, 1)
    fs.Solve()
    water_in = wet.MolarFlowMolPerSecond * wet.OverallMoleFraction("Water") * tc.MW_W / 1000  # kg/s
    lean = (fs.AddMaterialStream("lean").At(Q.Kelvin(t_gas + 5 + 273.15), kpa(P_CONTACTOR))
            .WithMassFlow(Q.KgPerSecond(ratio * water_in / lean_wt))
            .WithComposition(tc.composition(mole_from_wt(lean_wt))))
    dry, rich = fs.AddMaterialStream("dry"), fs.AddMaterialStream("rich")
    (fs.AddAbsorptionColumn("T-100").WithNumberOfStages(4).WithFeed(lean, 0).WithFeed(wet, 3)
       .WithTopProduct(dry).WithBottoms(rich).WithTopPressure(kpa(P_CONTACTOR - 50))
       .WithColumnPressureDrop(kpa(50)))
    fs.Solve()
    return (tc.y_to_mg_sm3(dry.OverallMoleFraction("Water")),
            tc.y_to_mg_sm3(wet.OverallMoleFraction("Water")))


rows = []
print("Contactor map (about 100 runs)...")
for lw, r, tg in itertools.product((0.985, 0.988, 0.991, 0.994, 0.997), (8, 12, 18, 26, 36), (30, 35, 40, 45)):
    try:
        dry, wet = contactor(lw, r, tg)
        rows.append(dict(lean_wt_pct=lw * 100, ratio=r, T_gas_C=tg, dry_mg_Sm3=dry, wet_mg_Sm3=wet))
    except Exception as e:
        print("failed", lw, r, tg, str(e).splitlines()[0][:60])
mp = pd.DataFrame(rows)
mp.round(3).to_csv(os.path.join(out, "07_contactor_map.csv"), index=False)


def features(lw, r, tg):
    a, b, c = (np.asarray(lw) - 99.0), np.log(np.asarray(r) / 18.0), (np.asarray(tg) - 40.0) / 5.0
    one = np.ones_like(a, dtype=float)
    return np.column_stack([one, a, b, c, a * a, b * b, c * c, a * b, a * c, b * c])


X = features(mp.lean_wt_pct, mp.ratio, mp.T_gas_C)
yv = np.log(mp.dry_mg_Sm3.values)
coef, *_ = np.linalg.lstsq(X, yv, rcond=None)
pred = np.exp(X @ coef)
err = (pred / mp.dry_mg_Sm3 - 1) * 100

# wet gas water vs T_gas (fit ln w = d0 + d1*T + d2*T^2) so the FMU can take T_gas as a disturbance
wg = mp.drop_duplicates("T_gas_C")[["T_gas_C", "wet_mg_Sm3"]]
cw = np.polyfit(wg.T_gas_C, np.log(wg.wet_mg_Sm3), 2)[::-1]

coeffs = {
    "bubble_T_C": {"form": "c0 + c1*ln(w) + c2*ln(w)^2, w = water mass fraction, 101.3 kPa",
                   "c": [float(v) for v in cb]},
    "contactor_ln_dry_mg_Sm3": {
        "form": "sum(k_i * f_i); a = leanWt% - 99, b = ln(ratio/18), c = (T_gas - 40)/5; "
                "f = [1, a, b, c, a^2, b^2, c^2, a*b, a*c, b*c]; ratio = kg TEG / kg water in wet gas",
        "k": [float(v) for v in coef],
        "valid": "lean 98.5-99.7 wt%, ratio 8-36, T_gas 30-45 C, 7000 kPa, lean TEG 5 C above gas"},
    "wet_gas_ln_mg_Sm3": {"form": "d0 + d1*T_gas + d2*T_gas^2 (7000 kPa)", "d": [float(v) for v in cw]},
}
json.dump(coeffs, open(os.path.join(tc.PROJECT, "data", "reduced_model_coeffs.json"), "w"), indent=2)

text = "\n".join([
    f"Bubble curve fit: max error {np.abs(TB_fit - TB).max():.2f} K over 0.3-10 wt% water; "
    f"T_bub at 0.85 wt% water = {cb[0] + cb[1] * math.log(0.0085) + cb[2] * math.log(0.0085) ** 2:.1f} C",
    f"Contactor map: {len(mp)} DWSIM runs, fit error mean {err.abs().mean():.1f} %, max {err.abs().max():.1f} %",
    f"Wet gas water: {', '.join(f'{t:.0f} C {w:.0f}' for t, w in zip(wg.T_gas_C, wg.wet_mg_Sm3))} mg/Sm3",
    "Coefficients: data/reduced_model_coeffs.json",
])
print(text)
open(os.path.join(out, "07_summary.txt"), "w").write(text + "\n")

"""
01_water_content.py  -  Step 1 of the TEG dehydration case.

Question: how much water does the raw gas carry into the TEG contactor,
and which DWSIM property package predicts it best?

Method
- Saturate a representative Queensland coal-seam gas with water in a flash
  drum (excess water) at contactor conditions.
- Compare the predicted water content with the Bukacek (1955) correlation,
  the standard hand method for lean sweet gas (about +/-5 %).
- Sweep temperature at fixed pressure and pressure at fixed temperature,
  for three cubic equations of state.

Outputs (results folder): 01_water_content.csv, 01_water_content.png,
01_summary.txt
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from dwsim_setup import load_dwsim, PROJECT

Flowsheet, PP, Q = load_dwsim()
from System import Action
from DWSIM.Automation.FluentAPI.Builders import CompositionBuilder

# ---- design basis (representative CSG, dry basis, mole fraction) ----------
GAS = {"Methane": 0.970, "Nitrogen": 0.015, "Carbon dioxide": 0.010, "Ethane": 0.005}
T_BASE_C, P_BASE_KPA = 40.0, 7000.0          # contactor inlet, kPa abs
SPEC_MG_SM3 = 112.0                          # AS 4564 / pipeline limit
V_STD = 8.314462 * 288.15 / 101325.0         # m3/mol at 15 C, 101.325 kPa
MW_W = 18.015                                # g/mol

PACKAGES = {"Peng-Robinson": PP.PengRobinson,
            "SRK": PP.SoaveRedlichKwong,
            "PRSV2-VL": PP.PRSV2VL}


def bukacek_mg_sm3(t_c, p_kpa):
    """Bukacek (1955) water content of sweet natural gas, mg/Sm3."""
    t_f = t_c * 9 / 5 + 32
    p_psia = p_kpa / 6.894757
    pv_kpa = 0.61121 * math.exp((18.678 - t_c / 234.5) * (t_c / (257.14 + t_c)))  # Buck
    pv_psia = pv_kpa / 6.894757
    b = 10 ** (-3083.87 / (459.6 + t_f) + 6.69449)
    w_lb_mmscf = 47484 * pv_psia / p_psia + b
    return w_lb_mmscf * 16.018


def dwsim_mg_sm3(pkg, t_c, p_kpa):
    """Saturate the gas with water in a flash drum; return water in vapour, mg/Sm3."""
    fs = Flowsheet.Create("sat")
    for c in list(GAS) + ["Water"]:
        fs = fs.WithCompound(c)
    fs = fs.WithPropertyPackage(pkg)

    def comp(c):
        for k, v in GAS.items():
            c.Mole(k, v * 0.95)
        c.Mole("Water", 0.05)          # excess water so a free-water phase exists

    feed = (fs.AddMaterialStream("feed").At(Q.Kelvin(t_c + 273.15), Q.Pascal(p_kpa * 1000))
            .WithMolarFlow(Q.MolPerSecond(100.0))
            .WithComposition(Action[CompositionBuilder](comp)))
    vap, liq = fs.AddMaterialStream("gas"), fs.AddMaterialStream("water")
    fs.AddSeparator("V-1").ConnectFeed(feed, 0).ConnectProduct(vap, 0).ConnectProduct(liq, 1)
    fs.Solve()
    y = vap.OverallMoleFraction("Water")
    return y / V_STD * MW_W * 1000.0


rows = []
for name, pkg in PACKAGES.items():
    for t in (20, 25, 30, 35, 40, 45, 50):
        rows.append(dict(sweep="T", package=name, T_C=t, P_kPa=P_BASE_KPA,
                         dwsim=dwsim_mg_sm3(pkg, t, P_BASE_KPA)))
    for p in (3000, 4000, 5000, 6000, 7000, 8000, 10000):
        rows.append(dict(sweep="P", package=name, T_C=T_BASE_C, P_kPa=p,
                         dwsim=dwsim_mg_sm3(pkg, T_BASE_C, p)))
df = pd.DataFrame(rows)
df["bukacek"] = [bukacek_mg_sm3(t, p) for t, p in zip(df.T_C, df.P_kPa)]
df["dev_pct"] = (df.dwsim / df.bukacek - 1) * 100

out = os.path.join(PROJECT, "results")
os.makedirs(out, exist_ok=True)
df.round(2).to_csv(os.path.join(out, "01_water_content.csv"), index=False)

# ---- plot ------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
for ax, sweep, x, xlabel in ((axes[0], "T", "T_C", "Temperature (°C), 7000 kPa"),
                             (axes[1], "P", "P_kPa", "Pressure (kPa abs), 40 °C")):
    d = df[df.sweep == sweep]
    first = d[d.package == "Peng-Robinson"]
    ax.plot(first[x], first.bukacek, "k--", lw=1.5, label="Bukacek correlation")
    for name in PACKAGES:
        if name == "PRSV2-VL":         # off the scale (missing parameters); reported in the summary
            continue
        dd = d[d.package == name]
        ax.plot(dd[x], dd.dwsim, marker="o", ms=4, label=f"DWSIM {name}")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Saturated water content (mg/Sm³)")
    ax.grid(alpha=0.3)
axes[0].legend(fontsize=8)
fig.suptitle("Raw gas carries far more water than the 112 mg/Sm³ pipeline limit")
fig.tight_layout()
fig.savefig(os.path.join(out, "01_water_content.png"), dpi=150)

# ---- summary ---------------------------------------------------------------
base = df[(df.sweep == "T") & (df.T_C == T_BASE_C)]
lines = [f"Contactor inlet {T_BASE_C:.0f} C, {P_BASE_KPA:.0f} kPa abs, gas {GAS}",
         f"Bukacek: {base.bukacek.iloc[0]:.0f} mg/Sm3; spec {SPEC_MG_SM3:.0f} mg/Sm3",
         "Mean |deviation| from Bukacek over all points:"]
for name in PACKAGES:
    d = df[df.package == name]
    b = base[base.package == name].iloc[0]
    lines.append(f"  {name:14s} base {b.dwsim:7.0f} mg/Sm3 ({b.dev_pct:+.1f} %), "
                 f"mean |dev| {d.dev_pct.abs().mean():.1f} %, max |dev| {d.dev_pct.abs().max():.1f} %")
best = min(PACKAGES, key=lambda n: df[df.package == n].dev_pct.abs().mean())
need = base.bukacek.iloc[0] - SPEC_MG_SM3
lines.append(f"Best match: {best}")
lines.append(f"Water to remove at base case: about {need:.0f} mg/Sm3 "
             f"({need / base.bukacek.iloc[0] * 100:.0f} % of inlet water)")
text = "\n".join(lines)
print(text)
open(os.path.join(out, "01_summary.txt"), "w").write(text + "\n")

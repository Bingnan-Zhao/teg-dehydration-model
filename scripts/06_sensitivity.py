"""
06_sensitivity.py  -  Step 6: which operating point meets spec at the lowest energy?

Uses the closed-loop unit model from 05_teg_unit.py and varies, one at a time
around the base case (25 L/kg, 4 stages, 200 C):
  A. TEG circulation, L TEG per kg water absorbed
  B. reboiler temperature
  C. contactor equilibrium stages
then maps circulation x reboiler temperature to find the cheapest point that
meets the 112 mg/Sm3 spec with a design margin (target <= 90 mg/Sm3).

Outputs: 06_sensitivity.csv, 06_sensitivity.png, 06_operating_window.png,
06_summary.txt  (about 5-10 minutes)
"""
import importlib.util
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("teg_unit", os.path.join(HERE, "05_teg_unit.py"))
unit = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(unit)
tc = unit.tc

BASE = dict(l_per_kg=25.0, stages=4, t_reboiler=200.0)
DESIGN_TARGET = 90.0          # mg/Sm3, ~20 % margin below the 112 mg/Sm3 spec

runs = []


def run(study, **kw):
    p = {**BASE, **kw}
    t0 = time.time()
    try:
        *_, m = unit.simulate(verbose=False, **p)
        m["status"] = "ok"
    except Exception as e:
        m = {"l_per_kg": p["l_per_kg"], "stages": p["stages"], "t_reboiler_C": p["t_reboiler"],
             "status": "failed: " + str(e).splitlines()[0][:80]}
    m["study"] = study
    runs.append(m)
    print(f"{study:8s} L={p['l_per_kg']:5.1f} N={p['stages']} Treb={p['t_reboiler']:5.1f} -> "
          f"{m.get('dry_gas_mg_Sm3', float('nan')):6.1f} mg/Sm3, "
          f"reboiler {m.get('reboiler_kW', float('nan')):5.1f} kW  ({time.time() - t0:.0f} s)")


for L in (10, 15, 20, 25, 30, 40):
    run("A_circ", l_per_kg=L)
for T in (180, 185, 190, 195, 200, 204):
    run("B_treb", t_reboiler=T)
for N in (4, 5, 6, 8):
    run("C_stages", stages=N)
for L in (10, 15, 20, 30):
    for T in (185, 195, 204):
        run("D_map", l_per_kg=L, t_reboiler=T)

df = pd.DataFrame(runs)
out = os.path.join(tc.PROJECT, "results")
df.round(3).to_csv(os.path.join(out, "06_sensitivity.csv"), index=False)
ok = df[df.status == "ok"]

# ---- one-at-a-time sensitivities ---------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(13, 4))
panels = (("A_circ", "l_per_kg", "TEG circulation (L TEG / kg water)"),
          ("B_treb", "t_reboiler_C", "Reboiler temperature (°C)"),
          ("C_stages", "stages", "Contactor equilibrium stages"))
for ax, (study, x, label) in zip(axes, panels):
    d = ok[(ok.study == study) | ((ok.study != "D_map") & (ok.l_per_kg == BASE["l_per_kg"])
                                  & (ok.stages == BASE["stages"]) & (ok.t_reboiler_C == BASE["t_reboiler"]))]
    d = d.drop_duplicates(subset=[x]).sort_values(x)
    ax.plot(d[x], d.dry_gas_mg_Sm3, "o-", color="tab:blue")
    ax.axhline(tc.SPEC_MG_SM3, color="tab:red", ls="--", lw=1)
    ax.text(d[x].min(), tc.SPEC_MG_SM3 * 1.04, "spec 112", color="tab:red", fontsize=8)
    ax.set_xlabel(label)
    ax.set_ylabel("Dry gas water (mg/Sm³)", color="tab:blue")
    ax.grid(alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(d[x], d.reboiler_kW, "s--", color="tab:orange", ms=4)
    ax2.set_ylabel("Reboiler duty (kW)", color="tab:orange")
fig.suptitle("Lean TEG purity (reboiler temperature) drives dry gas water; circulation mostly drives energy")
fig.tight_layout()
fig.savefig(os.path.join(out, "06_sensitivity.png"), dpi=150)

# ---- operating window -------------------------------------------------------
d = ok[ok.study.isin(["D_map", "A_circ", "B_treb"]) & (ok.stages == BASE["stages"])] \
    .drop_duplicates(subset=["l_per_kg", "t_reboiler_C"])
fig, ax = plt.subplots(figsize=(7, 5))
for T, g in d.groupby("t_reboiler_C"):
    if len(g) < 2:
        continue
    g = g.sort_values("l_per_kg")
    ax.plot(g.reboiler_kW, g.dry_gas_mg_Sm3, "o-", label=f"reboiler {T:.0f} °C")
    for _, r in g.iterrows():
        ax.annotate(f"{r.l_per_kg:.0f}", (r.reboiler_kW, r.dry_gas_mg_Sm3), fontsize=7,
                    xytext=(3, 3), textcoords="offset points")
ax.axhline(tc.SPEC_MG_SM3, color="tab:red", ls="--", lw=1, label="spec 112 mg/Sm³")
ax.axhline(DESIGN_TARGET, color="tab:red", ls=":", lw=1, label=f"design target {DESIGN_TARGET:.0f}")
ax.set_xlabel("Reboiler duty (kW)")
ax.set_ylabel("Dry gas water (mg/Sm³)")
ax.set_title("Operating window (labels: L TEG / kg water)")
ax.grid(alpha=0.3)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(out, "06_operating_window.png"), dpi=150)

# ---- summary ------------------------------------------------------------------
base = ok[(ok.study == "A_circ") & (ok.l_per_kg == BASE["l_per_kg"])].iloc[0]
feasible = d[d.dry_gas_mg_Sm3 <= DESIGN_TARGET].sort_values("reboiler_kW")
lines = [f"Base case (25 L/kg, 4 stages, 200 C): {base.dry_gas_mg_Sm3:.0f} mg/Sm3, "
         f"reboiler {base.reboiler_kW:.0f} kW, circulation {base.teg_circulation_kg_h:.0f} kg/h"]
if len(feasible):
    best = feasible.iloc[0]
    lines.append(f"Cheapest point meeting the {DESIGN_TARGET:.0f} mg/Sm3 design target: "
                 f"{best.l_per_kg:.0f} L/kg at {best.t_reboiler_C:.0f} C -> {best.dry_gas_mg_Sm3:.0f} mg/Sm3, "
                 f"reboiler {best.reboiler_kW:.0f} kW ({(1 - best.reboiler_kW / base.reboiler_kW) * 100:.0f} % "
                 f"less than base), circulation {best.teg_circulation_kg_h:.0f} kg/h")
else:
    lines.append("No mapped point meets the design target.")
for study, x in (("A_circ", "l_per_kg"), ("B_treb", "t_reboiler_C"), ("C_stages", "stages")):
    g = ok[ok.study == study].sort_values(x)
    lines.append(f"{study}: " + ", ".join(f"{v:g}->{w:.0f} mg ({k:.0f} kW)"
                                         for v, w, k in zip(g[x], g.dry_gas_mg_Sm3, g.reboiler_kW)))
failed = df[df.status != "ok"]
if len(failed):
    lines.append(f"{len(failed)} run(s) failed to converge; see 06_sensitivity.csv")
text = "\n".join(lines)
print(text)
open(os.path.join(out, "06_summary.txt"), "w").write(text + "\n")

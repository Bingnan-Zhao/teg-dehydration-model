"""
05_teg_unit.py  -  Step 5: the complete TEG dehydration unit (base case).

Flowsheet
  raw gas -> V-100 inlet scrubber -> T-100 contactor (4 equilibrium stages)
  rich TEG -> VLV-100 -> E-100 preheat -> V-101 flash drum (500 kPa)
           -> VLV-101 -> E-101 lean/rich exchanger (cold side) -> E-102 reboiler
           -> V-102 still (atmospheric, reboiler temperature) -> lean TEG
  lean TEG -> E-103 lean/rich exchanger (hot side) -> P-100 pump -> E-104 cooler
           -> back to the contactor top

The lean TEG loop is closed in Python (successive substitution on the lean
TEG composition), so every DWSIM solve is a simple once-through flowsheet.
Thermodynamics: Peng-Robinson with tuned kij values. TEG-water kij is strongly
temperature dependent, so the unit is solved as two DWSIM flowsheets:
  contactor section  (10-60 C)   TEG-water kij from step 4
  regeneration section (65-200 C) TEG-water kij from step 3

Outputs: flowsheets/05_contactor.dwxmz and flowsheets/05_regeneration.dwxmz
(open in DWSIM), results/05_contactor.png, results/05_regeneration.png,
results/05_streams.csv, results/05_summary.txt
"""
import os

import pandas as pd

import teg_common as tc
from teg_common import Q

# ---- design basis ----------------------------------------------------------
GAS_MOLS = 693.0             # mol/s dry gas, about 50 MMscfd (1.42 million Sm3/d)
P_CONTACTOR = 7000.0         # kPa abs
T_GAS, T_LEAN = 40.0, 45.0   # C; lean TEG kept ~5 C above the gas
STAGES = 4                   # equilibrium stages (DWSIM minimum is 4)
L_PER_KG_WATER = 25.0        # TEG circulation, L TEG per kg water absorbed (~3 US gal/lb)
P_FLASH, T_FLASH = 500.0, 65.0
T_RICH_TO_STILL = 150.0      # rich TEG leaving the lean/rich exchanger
T_REBOILER = 200.0           # kept below ~204 C to limit TEG degradation
P_STILL = 101.3
T_LEAN_AFTER_EXCH = 90.0
T_REFLUX = 102.0            # reflux condenser outlet, C (water stays vapour at 101 kPa)
RHO_TEG = 1.12               # kg/L
MW_TEG = 150.17

COMPS = list(tc.GAS) + ["Water", tc.TEG]
KIJ_C, KIJ_R = tc.tuned_kijs("contactor"), tc.tuned_kijs("regenerator")


def kpa(p):
    return Q.Pascal(p * 1000.0)


def degc(t):
    return Q.Kelvin(t + 273.15)


def teg_frac_to_mole(wt_teg):
    mT, mW = wt_teg / MW_TEG, (1 - wt_teg) / tc.MW_W
    return {tc.TEG: mT / (mT + mW), "Water": mW / (mT + mW)}


def build(lean_wt_teg, teg_kg_s, reflux):
    """reflux: (kg/s, TEG mass fraction) of the condensate returned to the still."""
    fs = fs1 = tc.new_flowsheet("TEG contactor section", COMPS, KIJ_C)
    S, U = {}, {}

    def stream(tag):
        S[tag] = fs.AddMaterialStream(tag)
        return S[tag]

    raw = {k: v * 0.98 for k, v in tc.GAS.items()}
    raw["Water"] = 0.02                                    # excess -> free water knocked out
    stream("raw gas").At(degc(T_GAS), kpa(P_CONTACTOR)).WithMolarFlow(
        Q.MolPerSecond(GAS_MOLS / 0.98)).WithComposition(tc.composition(raw))
    stream("wet gas"), stream("free water")
    fs.AddSeparator("V-100").ConnectFeed(S["raw gas"], 0).ConnectProduct(S["wet gas"], 0) \
        .ConnectProduct(S["free water"], 1)

    lean_mass = teg_kg_s / lean_wt_teg
    stream("lean TEG").At(degc(T_LEAN), kpa(P_CONTACTOR)).WithMassFlow(
        Q.KgPerSecond(lean_mass)).WithComposition(tc.composition(teg_frac_to_mole(lean_wt_teg)))
    stream("dry gas"), stream("rich TEG")
    (fs.AddAbsorptionColumn("T-100").WithNumberOfStages(STAGES)
       .WithFeed(S["lean TEG"], 0).WithFeed(S["wet gas"], STAGES - 1)
       .WithTopProduct(S["dry gas"]).WithBottoms(S["rich TEG"])
       .WithTopPressure(kpa(P_CONTACTOR - 50)).WithColumnPressureDrop(kpa(50)))

    fs1.AutoLayout()
    fs1.Solve()

    # regeneration section: same rich TEG, regenerator-range thermodynamics
    rich = S["rich TEG"]
    fs = fs2 = tc.new_flowsheet("TEG regeneration section", COMPS, KIJ_R)
    stream("rich TEG to regen").At(Q.Kelvin(rich.TemperatureK), Q.Pascal(rich.PressurePa)) \
        .WithMassFlow(Q.KgPerSecond(rich.MassFlowKgPerSecond)) \
        .WithComposition(tc.composition({c: rich.OverallMoleFraction(c) for c in COMPS}))
    stream("rich 1"), stream("rich 2"), stream("flash gas"), stream("rich 3")
    stream("rich 4"), stream("rich 5"), stream("rich 6")
    fs.AddValve("VLV-100").WithOutletPressure(kpa(P_FLASH + 20)).ConnectFeed(S["rich TEG to regen"]).ConnectProduct(S["rich 1"])
    U["E-100"] = fs.AddHeater("E-100").WithOutletTemperature(degc(T_FLASH)).WithPressureDrop(kpa(20)) \
        .ConnectFeed(S["rich 1"]).ConnectProduct(S["rich 2"])
    fs.AddSeparator("V-101").ConnectFeed(S["rich 2"], 0).ConnectProduct(S["flash gas"], 0) \
        .ConnectProduct(S["rich 3"], 1)
    fs.AddValve("VLV-101").WithOutletPressure(kpa(P_STILL + 30)).ConnectFeed(S["rich 3"]).ConnectProduct(S["rich 4"])
    U["E-101"] = fs.AddHeater("E-101").WithOutletTemperature(degc(T_RICH_TO_STILL)).WithPressureDrop(kpa(10)) \
        .ConnectFeed(S["rich 4"]).ConnectProduct(S["rich 5"])
    stream("reflux in").At(degc(T_REFLUX), kpa(P_STILL + 20)).WithMassFlow(Q.KgPerSecond(reflux[0])) \
        .WithComposition(tc.composition(teg_frac_to_mole(reflux[1])))
    stream("rich 5b")
    fs.AddMixer("MIX-100").ConnectFeed(S["rich 5"], 0).ConnectFeed(S["reflux in"], 1).ConnectProduct(S["rich 5b"], 0)
    U["E-102"] = fs.AddHeater("E-102").WithOutletTemperature(degc(T_REBOILER)).WithPressureDrop(kpa(20)) \
        .ConnectFeed(S["rich 5b"]).ConnectProduct(S["rich 6"])
    stream("still overhead"), stream("lean hot")
    fs.AddSeparator("V-102").ConnectFeed(S["rich 6"], 0).ConnectProduct(S["still overhead"], 0) \
        .ConnectProduct(S["lean hot"], 1)
    # reflux condenser: cool the overhead just above water's boiling point so TEG condenses
    stream("overhead 2"), stream("still vent"), stream("reflux out")
    U["E-105"] = fs.AddCooler("E-105").WithOutletTemperature(degc(T_REFLUX)).ConnectFeed(S["still overhead"]) \
        .ConnectProduct(S["overhead 2"])
    fs.AddSeparator("V-103").ConnectFeed(S["overhead 2"], 0).ConnectProduct(S["still vent"], 0) \
        .ConnectProduct(S["reflux out"], 1)

    stream("lean 1"), stream("lean 2"), stream("lean return")
    U["E-103"] = fs.AddCooler("E-103").WithOutletTemperature(degc(T_LEAN_AFTER_EXCH)).ConnectFeed(S["lean hot"]) \
        .ConnectProduct(S["lean 1"])
    U["P-100"] = fs.AddPump("P-100").WithOutletPressure(kpa(P_CONTACTOR + 50)).WithEfficiencyPercent(70) \
        .ConnectFeed(S["lean 1"]).ConnectProduct(S["lean 2"])
    U["E-104"] = fs.AddCooler("E-104").WithOutletTemperature(degc(T_LEAN)).WithPressureDrop(kpa(50)) \
        .ConnectFeed(S["lean 2"]).ConnectProduct(S["lean return"])
    fs2.AutoLayout()
    fs2.Solve()
    return (fs1, fs2), S, U


def wt_teg(s):
    return s.OverallMassFraction(tc.TEG)


# ---- close the lean TEG loop -------------------------------------------------
def duty(b):
    for name in ("HeatDutyKW", "HeatRemovedKW", "DutyKW"):
        if hasattr(b, name):
            return getattr(b, name)
    return b.Object.DeltaQ          # DWSIM object, kW


def simulate(l_per_kg=None, stages=None, t_reboiler=None, verbose=True):
    """Solve the unit with the lean TEG and reflux loops closed. Returns (fs, S, U, metrics)."""
    global L_PER_KG_WATER, STAGES, T_REBOILER
    if l_per_kg is not None:
        L_PER_KG_WATER = l_per_kg
    if stages is not None:
        STAGES = stages
    if t_reboiler is not None:
        T_REBOILER = t_reboiler
    lean_wt, teg_kg_s, reflux = 0.99, 0.5, (0.01, 0.9)
    for it in range(15):
        fs, S, U = build(lean_wt, teg_kg_s, reflux)
        water_removed = (S["wet gas"].MolarFlowMolPerSecond * S["wet gas"].OverallMoleFraction("Water")
                         - S["dry gas"].MolarFlowMolPerSecond * S["dry gas"].OverallMoleFraction("Water")) \
            * tc.MW_W / 1000.0                                   # kg/s
        teg_new = L_PER_KG_WATER * water_removed * RHO_TEG       # kg/s pure TEG
        lean_new = wt_teg(S["lean return"])
        r = S["reflux out"]
        reflux_new = (max(r.MassFlowKgPerSecond, 1e-5), min(max(wt_teg(r), 0.0), 0.999))
        if verbose:
            print(f"iter {it + 1}: lean TEG {lean_wt * 100:.3f} -> {lean_new * 100:.3f} wt%, "
                  f"TEG {teg_kg_s:.4f} -> {teg_new:.4f} kg/s")
        done = (abs(lean_new - lean_wt) < 2e-5 and abs(teg_new - teg_kg_s) / teg_kg_s < 1e-3
                and abs(reflux_new[0] - reflux[0]) < 1e-5)
        lean_wt, teg_kg_s, reflux = lean_new, teg_new, reflux_new
        if done:
            break
    fs, S, U = build(lean_wt, teg_kg_s, reflux)
    m = {
        "l_per_kg": L_PER_KG_WATER, "stages": STAGES, "t_reboiler_C": T_REBOILER,
        "dry_gas_mg_Sm3": tc.y_to_mg_sm3(S["dry gas"].OverallMoleFraction("Water")),
        "lean_teg_wt_pct": wt_teg(S["lean TEG"]) * 100,
        "rich_teg_wt_pct": wt_teg(S["rich TEG"]) * 100,
        "teg_circulation_kg_h": S["lean TEG"].MassFlowKgPerSecond * 3600,
        "reboiler_kW": U["E-102"].HeatDutyKW,
        "rich_heating_kW": U["E-100"].HeatDutyKW + U["E-101"].HeatDutyKW,
        "pump_kW": U["P-100"].PowerKW,
        "teg_loss_kg_d": sum(S[t].MassFlowKgPerSecond * S[t].OverallMassFraction(tc.TEG)
                             for t in ("dry gas", "flash gas", "still vent")) * 86400,
    }
    return fs, S, U, m


if __name__ == "__main__":
    fs, S, U, _ = simulate()
    # ---- results -------------------------------------------------------------------
    proj = tc.PROJECT
    os.makedirs(os.path.join(proj, "flowsheets"), exist_ok=True)
    for f, name in zip(fs, ("05_contactor", "05_regeneration")):
        f.Save(os.path.join(proj, "flowsheets", name + ".dwxmz"), True)
        try:
            f.SaveScreenshot(os.path.join(proj, "results", name + ".png"))
        except Exception as e:                                      # image is a nice-to-have
            print("flowsheet image not saved:", str(e).splitlines()[0])

    rows = []
    for tag, s in S.items():
        gas_like = tag in ("raw gas", "wet gas", "dry gas", "flash gas", "still vent", "still overhead", "overhead 2")
        rows.append({
            "stream": tag,
            "T_C": round(s.TemperatureK - 273.15, 1),
            "P_kPa": round(s.PressurePa / 1000, 1),
            "mass_kg_h": round(s.MassFlowKgPerSecond * 3600, 2),
            "water_mg_Sm3": round(tc.y_to_mg_sm3(s.OverallMoleFraction("Water")), 1) if gas_like else None,
            "TEG_wt_pct": None if gas_like else round(wt_teg(s) * 100, 3),
        })
    streams = pd.DataFrame(rows)
    streams.to_csv(os.path.join(proj, "results", "05_streams.csv"), index=False)

    dry = tc.y_to_mg_sm3(S["dry gas"].OverallMoleFraction("Water"))
    wet = tc.y_to_mg_sm3(S["wet gas"].OverallMoleFraction("Water"))
    teg_loss_kg_d = sum(S[t].MassFlowKgPerSecond * S[t].OverallMassFraction(tc.TEG)
                        for t in ("dry gas", "flash gas", "still vent")) * 86400
    text = "\n".join([
        f"Gas {GAS_MOLS * tc.V_STD * 86400 / 1e6:.2f} million Sm3/d at {P_CONTACTOR:.0f} kPa, {T_GAS:.0f} C; "
        f"{STAGES} stages; circulation {L_PER_KG_WATER:.0f} L TEG/kg water",
        f"Wet gas {wet:.0f} mg/Sm3 -> dry gas {dry:.0f} mg/Sm3 (spec {tc.SPEC_MG_SM3:.0f}): "
        f"{'MEETS SPEC' if dry <= tc.SPEC_MG_SM3 else 'OFF SPEC'}",
        f"Lean TEG {wt_teg(S['lean TEG']) * 100:.2f} wt% (reboiler {T_REBOILER:.0f} C, {P_STILL:.0f} kPa), "
        f"rich TEG {wt_teg(S['rich TEG']) * 100:.2f} wt%",
        f"TEG circulation {S['lean TEG'].MassFlowKgPerSecond * 3600:.0f} kg/h; "
        f"TEG losses {teg_loss_kg_d:.1f} kg/d",
        f"Duties kW: rich preheat E-100 {U['E-100'].HeatDutyKW:.0f}, lean/rich cold side E-101 {U['E-101'].HeatDutyKW:.0f}, "
        f"reboiler E-102 {U['E-102'].HeatDutyKW:.0f}, lean/rich hot side E-103 {duty(U['E-103']):.0f}, "
        f"lean cooler E-104 {duty(U['E-104']):.0f}, reflux condenser E-105 {duty(U['E-105']):.0f}; pump {U['P-100'].PowerKW:.1f} kW",
        f"Flowsheets saved: flowsheets/05_contactor.dwxmz, flowsheets/05_regeneration.dwxmz",
    ])
    print(text)
    open(os.path.join(proj, "results", "05_summary.txt"), "w").write(text + "\n")

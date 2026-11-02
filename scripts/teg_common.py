"""
teg_common.py  -  shared design basis and helpers for the TEG dehydration case.
"""
import json
import math
import os

from dwsim_setup import load_dwsim, PROJECT

Flowsheet, PP, Q = load_dwsim()
from System import Action
from DWSIM.Automation.FluentAPI.Builders import CompositionBuilder
from DWSIM.Automation.FluentAPI import PropertyPackageBuilder, PRConfig, PropertyPackageConfigExtensions

# ---- design basis ----------------------------------------------------------
# Representative Queensland coal-seam gas, dry basis, mole fraction
GAS = {"Methane": 0.970, "Nitrogen": 0.015, "Carbon dioxide": 0.010, "Ethane": 0.005}
TEG = "Triethylene glycol"
T_BASE_C, P_BASE_KPA = 40.0, 7000.0          # contactor inlet, kPa abs
SPEC_MG_SM3 = 112.0                          # pipeline limit (QGP spec, AS 4564)
V_STD = 8.314462 * 288.15 / 101325.0         # m3/mol at 15 C, 101.325 kPa
MW_W = 18.015                                # g/mol

# Methane-water kij for Peng-Robinson, tuned by 02_tune_kij.py (DWSIM default is 0.5)
_KIJ_FILE = os.path.join(PROJECT, "results", "02_kij.json")
KIJ_CH4_H2O = json.load(open(_KIJ_FILE))["kij_methane_water"] if os.path.exists(_KIJ_FILE) else None


def bukacek_mg_sm3(t_c, p_kpa):
    """Bukacek (1955) water content of sweet natural gas, mg/Sm3."""
    t_f = t_c * 9 / 5 + 32
    p_psia = p_kpa / 6.894757
    pv_kpa = 0.61121 * math.exp((18.678 - t_c / 234.5) * (t_c / (257.14 + t_c)))  # Buck
    b = 10 ** (-3083.87 / (459.6 + t_f) + 6.69449)
    return (47484 * (pv_kpa / 6.894757) / p_psia + b) * 16.018


def y_to_mg_sm3(y):
    return y / V_STD * MW_W * 1000.0


_TEG_FILE = os.path.join(PROJECT, "results", "03_kij_teg.json")
KIJ_TEG_H2O = json.load(open(_TEG_FILE))["kij_teg_water"] if os.path.exists(_TEG_FILE) else None


_TEGC_FILE = os.path.join(PROJECT, "results", "04_kij_contactor.json")
KIJ_TEG_H2O_CONTACTOR = (json.load(open(_TEGC_FILE))["kij_teg_water_contactor"]
                         if os.path.exists(_TEGC_FILE) else None)


def tuned_kijs(zone="regenerator"):
    """kij overrides fitted so far, as {(c1, c2): kij}.
    zone: "contactor" (10-60 C, step 4 fit) or "regenerator" (100-200 C, step 3 fit)."""
    k = {}
    if KIJ_CH4_H2O is not None:
        k[("Methane", "Water")] = KIJ_CH4_H2O
    teg = KIJ_TEG_H2O_CONTACTOR if zone == "contactor" and KIJ_TEG_H2O_CONTACTOR is not None else KIJ_TEG_H2O
    if teg is not None:
        k[(TEG, "Water")] = teg
    return k


def new_flowsheet(name, compounds, kij=None):
    """Flowsheet with Peng-Robinson.
    kij: None (DWSIM defaults), a number (methane-water kij) or {(c1, c2): kij}."""
    fs = Flowsheet.Create(name)
    for c in compounds:
        fs = fs.WithCompound(c)
    if kij is None:
        return fs.WithPropertyPackage(PP.PengRobinson)
    pairs = kij if isinstance(kij, dict) else {("Methane", "Water"): kij}
    pairs = {p: v for p, v in pairs.items() if p[0] in compounds and p[1] in compounds}

    def pr(c):
        for (a, b), v in pairs.items():
            c.WithKij(a, b, float(v))

    def cfg(pp):
        pp.ConfigurePR(Action[PRConfig](pr))
    return PropertyPackageConfigExtensions.WithPropertyPackage(
        fs, PP.PengRobinson, Action[PropertyPackageBuilder](cfg))


def composition(fracs):
    def comp(c):
        for k, v in fracs.items():
            c.Mole(k, v)
    return Action[CompositionBuilder](comp)


def saturated_water(t_c, p_kpa, kij=KIJ_CH4_H2O):
    """Water content (mg/Sm3) of the gas saturated with free water at t_c, p_kpa."""
    fs = new_flowsheet("sat", list(GAS) + ["Water"], kij)
    fr = {k: v * 0.95 for k, v in GAS.items()}
    fr["Water"] = 0.05
    feed = (fs.AddMaterialStream("feed").At(Q.Kelvin(t_c + 273.15), Q.Pascal(p_kpa * 1000))
            .WithMolarFlow(Q.MolPerSecond(100.0)).WithComposition(composition(fr)))
    vap, liq = fs.AddMaterialStream("gas"), fs.AddMaterialStream("water")
    fs.AddSeparator("V-1").ConnectFeed(feed, 0).ConnectProduct(vap, 0).ConnectProduct(liq, 1)
    fs.Solve()
    return y_to_mg_sm3(vap.OverallMoleFraction("Water"))

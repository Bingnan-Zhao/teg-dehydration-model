"""
08_python_twin.py  -  Python twin of the Modelica model TEGDehydration.TEGUnit.

Same equations and parameters as modelica/TEGDehydration/package.mo, solved with
scipy, running the same step test as TEGDehydration.Examples.StepTest. Use it
to check the Dymola/FMU results (they should match to plotting accuracy) and as
a reference that runs without Modelica.

Outputs: results/08_twin_steptest.csv, results/08_twin_steptest.png, results/08_summary.txt
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
C = json.load(open(os.path.join(PROJECT, "data", "reduced_model_coeffs.json")))
cBub = C["bubble_T_C"]["c"]
kMap = C["contactor_ln_dry_mg_Sm3"]["k"]
dWet = C["wet_gas_ln_mg_Sm3"]["d"]

P = dict(cp=2600.0, dHvap=2.26e6, T_richIn=150.0, M_weir=550.0, k_weir=0.01, k_evap=0.05,
         M_surgeMax=3000.0, tau_contactor=180.0, tau_rich=600.0, eps=0.05)


def softplus(x, eps):
    return 0.5 * (x + np.sqrt(x * x + eps * eps))


def algebra(x, u):
    M_reb, w_reb, T_reb, M_surge, w_surge, w_lean, m_rich, w_rich = x
    Q, m_circ, F_gas, T_gas = u
    T_bub = cBub[0] + cBub[1] * np.log(w_reb) + cBub[2] * np.log(w_reb) ** 2
    V = P["k_evap"] * softplus(T_reb - T_bub, P["eps"])
    m_over = P["k_weir"] * softplus(M_reb - P["M_weir"], P["eps"])
    wet = np.exp(dWet[0] + dWet[1] * T_gas + dWet[2] * T_gas ** 2)
    water_in = F_gas * 1e6 / 86400 * wet * 1e-6
    ratio = m_circ * (1 - w_lean) / max(water_in, 1e-6)
    a = min(max(100 * (1 - w_lean) - 99, -1.0), 0.9)
    b = np.log(min(max(ratio, 6), 45) / 18)
    c = (T_gas - 40) / 5
    f = [1, a, b, c, a * a, b * b, c * c, a * b, a * c, b * c]
    dry = min(np.exp(np.dot(kMap, f)), wet)
    water_abs = F_gas * 1e6 / 86400 * (wet - dry) * 1e-6
    m_now = m_circ + water_abs
    w_now = (m_circ * w_lean + water_abs) / m_now
    return dict(T_bub=T_bub, V=V, m_over=m_over, wet=wet, dry=dry, m_now=m_now, w_now=w_now)


def rhs(t, x, ufun):
    M_reb, w_reb, T_reb, M_surge, w_surge, w_lean, m_rich, w_rich = x
    Q, m_circ, F_gas, T_gas = u = ufun(t)
    a = algebra(x, u)
    return [
        m_rich - a["m_over"] - a["V"],
        (m_rich * (w_rich - w_reb) - a["V"] * (1 - w_reb)) / M_reb,
        (m_rich * P["cp"] * (P["T_richIn"] - T_reb) + Q - a["V"] * P["dHvap"]) / (M_reb * P["cp"]),
        a["m_over"] - m_circ,
        a["m_over"] * (w_reb - w_surge) / M_surge,
        (w_surge - w_lean) / P["tau_contactor"],
        (a["m_now"] - m_rich) / P["tau_rich"],
        (a["w_now"] - w_rich) / P["tau_rich"],
    ]


def steady_state(u, level=0.5):
    M_s = level * P["M_surgeMax"]

    def res(z):
        x = [z[0], z[1], z[2], M_s, z[3], z[4], z[5], z[6]]
        d = rhs(0, x, lambda t: u)
        return [d[0], d[1] * 1e3, d[2], d[4] * 1e3, d[5] * 1e3, d[6], d[7] * 1e3]

    z0 = [560, 0.0085, 200, 0.0085, 0.0085, 0.49, 0.04]
    z = fsolve(res, z0, xtol=1e-12)
    return np.array([z[0], z[1], z[2], M_s, z[3], z[4], z[5], z[6]])


def u_step(t):
    Q = 1.0e5 - (1.5e4 if t >= 2 * 3600 else 0)
    m = 0.473 - (0.15 if t >= 8 * 3600 else 0)
    F = 1.42 + (0.3 if t >= 14 * 3600 else 0)
    T = 40 + (5 if t >= 20 * 3600 else 0)
    return Q, m, F, T


if __name__ == "__main__":
    x0 = steady_state(u_step(0))
    t_end = 26 * 3600
    t_eval = np.arange(0, t_end + 1, 60.0)
    sol = solve_ivp(rhs, (0, t_end), x0, args=(u_step,), method="BDF", t_eval=t_eval,
                    max_step=120, rtol=1e-7, atol=1e-9)
    rows = []
    for t, x in zip(sol.t, sol.y.T):
        u = u_step(t)
        a = algebra(x, u)
        rows.append(dict(time_h=t / 3600, Q_kW=u[0] / 1e3, m_circ=u[1], F_gas=u[2], T_gas=u[3],
                         dryGasWater=a["dry"], T_reboiler=x[2], leanTEG=100 * (1 - x[5]),
                         surgeLevel=x[3] / P["M_surgeMax"], richTEG=100 * (1 - x[7]), boilUp=a["V"]))
    df = pd.DataFrame(rows)
    df.round(5).to_csv(os.path.join(PROJECT, "results", "08_twin_steptest.csv"), index=False)

    fig, ax = plt.subplots(5, 1, figsize=(9, 11), sharex=True)
    ax[0].plot(df.time_h, df.dryGasWater)
    ax[0].axhline(112, color="tab:red", ls="--", lw=1)
    ax[0].set_ylabel("Dry gas water\n(mg/Sm³)")
    ax[1].plot(df.time_h, df.T_reboiler)
    ax[1].set_ylabel("Reboiler T (°C)")
    ax[2].plot(df.time_h, df.leanTEG, label="lean")
    ax[2].plot(df.time_h, df.richTEG, label="rich")
    ax[2].set_ylabel("TEG (wt%)")
    ax[2].legend(fontsize=8)
    ax[3].plot(df.time_h, df.surgeLevel * 100)
    ax[3].set_ylabel("Surge level (%)")
    ax[4].step(df.time_h, df.Q_kW, where="post", label="reboiler duty (kW)")
    ax[4].step(df.time_h, df.m_circ * 200, where="post", label="circulation (kg/s x200)")
    ax[4].step(df.time_h, df.F_gas * 50, where="post", label="gas flow (MSm³/d x50)")
    ax[4].step(df.time_h, df.T_gas * 2, where="post", label="gas T (°C x2)")
    ax[4].set_ylabel("Inputs (scaled)")
    ax[4].legend(fontsize=7, ncol=2)
    ax[4].set_xlabel("Time (h)")
    for a_ in ax:
        a_.grid(alpha=0.3)
        for s in (2, 8, 14, 20):
            a_.axvline(s, color="gray", lw=0.6, ls=":")
    fig.suptitle("Step test: -15 kW duty at 2 h, -0.15 kg/s circulation at 8 h,\n"
                 "+0.3 MSm³/d gas at 14 h, +5 K gas temperature at 20 h")
    fig.tight_layout()
    fig.savefig(os.path.join(PROJECT, "results", "08_twin_steptest.png"), dpi=150)

    def at(h):
        return df.iloc[(df.time_h - h).abs().idxmin()]
    lines = [f"Initial steady state: T_reb {x0[2]:.1f} C, lean {100 * (1 - x0[5]):.2f} wt%, "
             f"dry gas {at(0).dryGasWater:.0f} mg/Sm3, boil-up {at(0).boilUp * 3600:.1f} kg/h"]
    for h, label in ((7.9, "after -15 kW duty"), (13.9, "after -0.15 kg/s circulation"),
                     (19.9, "after +0.3 MSm3/d gas"), (26, "after +5 K gas temperature")):
        r = at(h)
        lines.append(f"{label:30s}: dry {r.dryGasWater:5.0f} mg/Sm3, T_reb {r.T_reboiler:5.1f} C, "
                     f"lean {r.leanTEG:5.2f} wt%, surge {r.surgeLevel * 100:4.1f} %")
    text = "\n".join(lines)
    print(text)
    open(os.path.join(PROJECT, "results", "08_summary.txt"), "w").write(text + "\n")

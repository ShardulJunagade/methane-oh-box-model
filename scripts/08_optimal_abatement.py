"""Step 8 - Optimal methane abatement as a dynamic optimisation problem (L19, GEKKO; L18 economics).

State:    dC/dt = [E_b(t)(1 - u) + E_nat]/beta - k P(t) C/(a C/C_ref + 1 - a) - k_o C   (best-fit model)
Control:  u(t) in [0, 0.7], the fraction of baseline (SSP2-4.5) anthropogenic emissions abated.
Cost:     linear marginal abatement cost MAC(u) = m u ($/t), so the yearly cost is 0.5 m u^2 E_b.
          m = 4000 $/t: 30 % abatement costs 1200 $/t at the margin (illustrative; UNEP/CCAC 2021 finds
          ~30 % achievable this decade, much of it cheap).

Problem A (cost-benefit, like the L19 fishery):  min  sum_t e^{-rt} [ 0.5 m u^2 E_b + D C ]
          with D ($ per ppb per year) set so that one tonne emitted today costs the social cost of
          methane SC: D = SC * beta * 1e6 * (r + 1/tau_pert).
Problem B (target):  min  sum_t e^{-rt} 0.5 m u^2 E_b   subject to  C(t) <= C_GMP(2050) for t >= 2050.
Both are solved with GEKKO (IMODE 6, IPOPT), with r = 2 %/yr.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings

from gekko import GEKKO

warnings.filterwarnings("ignore", category=DeprecationWarning)

from analysis import common as cm
from analysis import models as M
from analysis import scenarios as SC

print("08 optimal abatement (GEKKO)")
bf = cm.load_json("04_best_fit")
p = M.Params(k_oh=bf["k_oh"], e_nat=bf["e_nat"], a=bf["a"], gN=bf["gN"], gC=bf["gC"], q10=bf["q10"],
             arrhenius=bf["arrhenius"])
tau_p = cm.load_json("03_oh_feedback")["lifetimes"]["2025"]["tau_perturbation"]
R, MS, UMAX = 0.02, 4000.0, 0.7
years = np.arange(2025, 2101)
t = years - 2025.0

e_b, drv = SC.scenario("ssp245")
Eb, P = M.yearly_coeffs(p, years, e_b, drv)
Eb = Eb - p.e_nat   # yearly_coeffs adds the natural source; keep anthropogenic baseline only (q10 = 1)
_, jan1 = M.simulate(p, 1984, 2024, M.noaa_jan1(1984), *SC.scenario("ssp245"))
C2025 = float(jan1[-1])
disc = np.exp(-R * t)


def solve(problem, SCm=1600.0, target=None, u_fixed=None):
    """problem 'A' or 'B'; or 'sim' to simulate a given abatement path with the same discretisation."""
    m = GEKKO(remote=False)
    m.time = t
    E = m.Param(value=Eb)
    Pp = m.Param(value=P)
    dsc = m.Param(value=disc)
    if problem == "sim":
        u = m.Param(value=u_fixed)
    else:
        u = m.MV(value=0.0, lb=0, ub=UMAX)   # first point = 2025, no abatement yet
        u.STATUS = 1
    C = m.Var(value=C2025)
    m.Equation(C.dt() == (E * (1 - u) + p.e_nat) / cm.BETA
               - p.k_oh * Pp * C / (p.a * C / M.C_REF + 1 - p.a) - cm.K_OTHER * C)
    cost = m.Intermediate(0.5 * MS * u ** 2 * E * 1e-3)                 # $ billion per year
    if problem == "sim":
        m.options.IMODE = 4
        m.options.NODES = 3
        m.solve(disp=False)
        return np.array(u.value), np.array(C.value), np.array(cost.value)
    if problem == "A":
        D = SCm * cm.BETA * 1e6 * (R + 1 / tau_p) * 1e-9                # $ billion per ppb per year
        m.Minimize(dsc * (cost + D * C))
    else:
        mask = m.Param(value=(years >= 2050).astype(float))
        m.Equation(mask * (C - target) <= 0)
        m.Minimize(dsc * cost)
    m.options.IMODE = 6
    m.options.SOLVER = 3
    m.options.NODES = 3
    m.solve(disp=False)
    return np.array(u.value), np.array(C.value), np.array(cost.value)


def path_cost(u):
    return float((disc * 0.5 * MS * np.asarray(u) ** 2 * Eb * 1e-3).sum())


out = {"settings": {"discount_rate": R, "MAC_slope_usd_per_t": MS, "u_max": UMAX, "tau_pert_yr": tau_p,
                    "baseline": "SSP2-4.5"}}

# GMP path expressed as abatement fraction of the baseline
gmp = SC.gmp_path()
u_gmp = np.r_[0.0, np.clip(1 - gmp.values / Eb[1:], 0, UMAX)]
_, C_gmp, _ = solve("sim", u_fixed=u_gmp)   # same discretisation as the optimiser
s_gmp = pd.Series(C_gmp, index=years)
out["gmp"] = {"u_2030": float(u_gmp[5]), "u_2050": float(u_gmp[25]), "C_2050": float(s_gmp.loc[2050]),
              "pv_cost_usd_bn": path_cost(u_gmp)}

resA = {}
for SCm in (800.0, 1600.0, 3200.0):
    u, C, _ = solve("A", SCm)
    resA[SCm] = (u, C)
    out.setdefault("problem_A", {})[f"SC_{SCm:.0f}"] = {
        "u_2030": float(u[5]), "u_2050": float(u[25]), "u_2090": float(u[65]), "u_2100": float(u[-1]),
        "analytic_u_star": min(SCm / MS, UMAX), "C_2050": float(C[25]), "pv_cost_usd_bn": path_cost(u)}
    print(f"  A  SC-CH4 ${SCm:.0f}/t: u(2030)={u[5]:.2f}, u(2050)={u[25]:.2f}, u(2100)={u[-1]:.2f} "
          f"(analytic {min(SCm / MS, UMAX):.2f}), C2050={C[25]:.0f}")

target = out["gmp"]["C_2050"]
uB, CB, _ = solve("B", target=target)
out["problem_B"] = {"target_C": target, "u_2030": float(uB[5]), "u_2040": float(uB[15]), "u_2050": float(uB[25]),
                    "u_2070": float(uB[45]), "u_peak": float(uB.max()), "u_peak_year": int(years[uB.argmax()]),
                    "C_2050": float(CB[25]), "pv_cost_usd_bn": path_cost(uB),
                    "cost_ratio_gmp_over_optimal": out["gmp"]["pv_cost_usd_bn"] / path_cost(uB)}
print(f"  B  reach C(2050) = {target:.0f}: u(2030)={uB[5]:.2f}, u(2050)={uB[25]:.2f}; "
      f"GMP costs {out['problem_B']['cost_ratio_gmp_over_optimal']:.2f}x the optimum")
cm.save_json(out, "08_optimal_abatement")

fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.4))
for i, SCm in enumerate(resA):
    ax[0].plot(years, resA[SCm][0], color=cm.SERIES[i], label=r"Cost-benefit, social cost \$" + f"{SCm:.0f}/t")
ax[0].plot(years, uB, color=cm.SERIES[3], label="Cheapest way to meet GMP's 2050 level")
ax[0].plot(years, u_gmp, color=cm.INK, ls="--", lw=1.2, label="Global Methane Pledge")
ax[0].set_ylabel("Abated fraction u(t)")
ax[0].set_title("Optimal abatement paths (GEKKO)")
ax[0].legend(fontsize=6.8)
_, base, _ = solve("sim", u_fixed=np.zeros(len(years)))
ax[1].plot(years, base, color=cm.MUTED, label="No abatement (SSP2-4.5)")
for i, SCm in enumerate(resA):
    ax[1].plot(years, resA[SCm][1], color=cm.SERIES[i], lw=1.3)
ax[1].plot(years, CB, color=cm.SERIES[3], lw=1.3)
ax[1].plot(s_gmp.index, s_gmp.values, color=cm.INK, ls="--", lw=1.2, label="Global Methane Pledge")
ax[1].set_ylabel("CH$_4$ (ppb)")
ax[1].set_title("Resulting methane")
ax[1].legend(fontsize=7)
cm.savefig(fig, "08_optimal_abatement")

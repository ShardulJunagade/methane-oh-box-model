"""Step 3 - Theory of the CH4-OH system (L04 stability/feedback/tipping, L15 mass-action kinetics).

Two-variable mass-action model in normalised units (c = C/C_ref, X = OH/OH_ref):

    dC/dt = E/beta - k X C - k_o C                      (CH4 + OH -> products, plus other sinks)
    dX/dt = [P - X (a C/C_ref + 1 - a)] / tau_X          (OH made at rate P, destroyed by CH4 and by CO/VOC)

tau_X ~ 1 s, so the system is extremely stiff: X relaxes 3e8 times faster than C. Setting dX/dt = 0
(quasi-steady state, as in enzyme kinetics) gives the one-equation model with a saturating sink.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit

from analysis import common as cm
from analysis import models as M

print("03 CH4-OH feedback theory")
k, a, ko, Cr, beta = M.K_OH_PRIOR[0], cm.A_FEEDBACK, cm.K_OTHER, M.C_REF, cm.BETA
TAU_X = 1.0 / 3.156e7  # yr  (1 second)
out = {"k_oh": k, "a": a, "k_other": ko, "C_ref": Cr}


def sink(C, P=1.0, a=a, ko=ko):           # ppb/yr, quasi-steady OH
    return k * P * C / (a * C / Cr + 1 - a) + ko * C


def dsink(C, P=1.0, a=a, ko=ko):
    return k * P * (1 - a) / (a * C / Cr + 1 - a) ** 2 + ko


def full_rhs(t, y, E, P=1.0, tau_x=TAU_X):  # noqa: D401
    C, X = y
    return [E / beta - k * X * C - ko * C, (P - X * (a * C / Cr + 1 - a)) / tau_x]


# --- (a) lifetimes and the feedback factor ------------------------------------------------
C_now = float(cm.noaa_annual().loc[2025])
for lab, C in [("2010 (reference)", Cr), ("2025", C_now)]:
    tau_g, tau_p = C / sink(C), 1 / dsink(C)
    out.setdefault("lifetimes", {})[lab] = {"C": C, "tau_global": tau_g, "tau_perturbation": tau_p,
                                            "feedback_factor": tau_p / tau_g,
                                            "OH_relative": 1 / (a * C / Cr + 1 - a)}
    print(f"  {lab}: tau={tau_g:.2f} yr, tau_pert={tau_p:.2f} yr, f={tau_p / tau_g:.3f}")

# --- (b) Jacobian of the full system at equilibrium (L04) -----------------------------------
E_now = beta * sink(C_now)        # emission that would hold C at today's level
X_eq = 1 / (a * C_now / Cr + 1 - a)
J = np.array([[-k * X_eq - ko, -k * C_now],
              [-a / Cr * X_eq / TAU_X, -(a * C_now / Cr + 1 - a) / TAU_X]])
ev = np.sort(np.linalg.eigvals(J).real)
out["jacobian"] = {"E_equilibrium_Tg": E_now, "eigenvalues_per_yr": ev.tolist(),
                   "fast_timescale_s": float(-1 / ev[0] * 3.156e7), "slow_timescale_yr": float(-1 / ev[1]),
                   "stiffness_ratio": float(ev[0] / ev[1])}
print(f"  eigenvalues {ev}: fast {-1 / ev[0] * 3.156e7:.2f} s, slow {-1 / ev[1]:.2f} yr")

# QSSA vs full stiff system after a 20 % emission cut
sol = solve_ivp(full_rhs, (0, 40), [C_now, X_eq], args=(0.8 * E_now,), method="Radau",
                rtol=1e-8, atol=[1e-6, 1e-10], dense_output=True)
tt = np.linspace(0, 40, 401)
q = solve_ivp(lambda t, y: [0.8 * E_now / beta - sink(y[0])], (0, 40), [C_now], t_eval=tt, rtol=1e-9)
out["qssa_vs_full_max_diff_ppb"] = float(np.abs(sol.sol(tt)[0] - q.y[0]).max())

# --- (c) pulse experiment --------------------------------------------------------------------
pulse = 100 / beta   # 100 Tg
res = {}
for lab, aa in [("with OH feedback", a), ("constant OH", 0.0)]:
    kk = k if aa else k   # same reference lifetime at C_ref
    E_eq = beta * sink(C_now, a=aa)
    s = solve_ivp(lambda t, y: [E_eq / beta - sink(y[0], a=aa)], (0, 60), [C_now + pulse], t_eval=tt * 1.5, rtol=1e-9)
    ex = s.y[0] - C_now
    (efold,), _ = curve_fit(lambda t, T: pulse * np.exp(-t / T), s.t, ex, p0=[10])
    res[lab] = {"t": s.t, "excess": ex, "efold": efold}
    out.setdefault("pulse_efold_yr", {})[lab] = float(efold)
print("  pulse e-folding:", {k_: round(v["efold"], 2) for k_, v in res.items()})

# --- (d) equilibrium vs emissions: the OH-runaway threshold ------------------------------
E_grid = np.linspace(100, 2500, 600)
curves = {}
for lab, aa, ko_, P in [("Realistic (OH feedback + other sinks)", a, ko, 1.0),
                        ("OH sink only (no strat./soil/Cl sinks)", a, 0.0, 1.0),
                        ("OH sink only, 20 % less NOx-driven OH source", a, 0.0, 0.8),
                        ("Constant OH (no feedback)", 0.0, ko, 1.0)]:
    Cs = []
    for E in E_grid:
        f = lambda C: E / beta - sink(C, P=P, a=aa, ko=ko_)
        try:
            Cs.append(brentq(f, 1, 1e6))
        except ValueError:
            Cs.append(np.nan)
    curves[lab] = np.array(Cs)
E_crit = beta * k * Cr / a           # OH-only system: no equilibrium above this
E_crit_lowP = beta * k * 0.8 * Cr / a
out["threshold"] = {"E_crit_OH_only_Tg": E_crit, "E_crit_OH_only_lowNOx_Tg": E_crit_lowP,
                    "today_total_emission_Tg": 575.0, "ratio": E_crit / 575.0,
                    "Ceq_at_today_emission_realistic": float(brentq(lambda C: 575 / beta - sink(C), 1, 1e6)),
                    "Ceq_at_today_emission_noFeedback": float(brentq(lambda C: 575 / beta - sink(C, a=0), 1, 1e6)),
                    "dCeq_dE_ppb_per_Tg_today": float(1 / (beta * dsink(C_now)))}
print(f"  OH-only runaway threshold {E_crit:.0f} Tg/yr = {E_crit / 575:.1f} x today's emissions")

# --- figures ----------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.3))
for i, (lab, r) in enumerate(res.items()):
    ax[0].plot(r["t"], r["excess"], color=cm.SERIES[i], label=f"{lab} (e-fold {r['efold']:.1f} yr)")
ax[0].set_xlabel("Years after pulse")
ax[0].set_ylabel("Excess CH$_4$ (ppb)")
ax[0].set_title("A 100 Tg pulse decays slower with feedback")
ax[0].legend()
Cg = np.linspace(500, 6000, 300)
ax[1].plot(Cg, Cg / sink(Cg), color=cm.SERIES[0], label="Global lifetime C/S(C)")
ax[1].plot(Cg, 1 / dsink(Cg), color=cm.SERIES[1], label="Perturbation lifetime 1/S'(C)")
ax[1].axvline(C_now, color=cm.MUTED, lw=0.8, ls=":")
ax[1].text(C_now + 60, 9.6, "2025", fontsize=7.5, color=cm.INK2)
ax[1].set_xlabel("CH$_4$ (ppb)")
ax[1].set_ylabel("Years")
ax[1].set_title("Methane lengthens its own lifetime")
ax[1].legend()
cm.savefig(fig, "03_pulse_lifetime")

fig, ax = plt.subplots(figsize=(7.0, 4.0))
for i, (lab, Cs) in enumerate(curves.items()):
    ax.plot(E_grid, Cs, color=cm.SERIES[i], label=lab)
ax.axvline(E_crit, color=cm.SERIES[1], lw=0.8, ls="--")
ax.axvline(E_crit_lowP, color=cm.SERIES[2], lw=0.8, ls="--")
ax.text(E_crit + 25, 11300, f"OH runaway\nthreshold\n{E_crit:.0f} Tg/yr", fontsize=7, color=cm.INK2, va="top")
ax.text(E_crit_lowP - 25, 11300, f"with 20 %\nless OH source\n{E_crit_lowP:.0f} Tg/yr", fontsize=7,
        color=cm.INK2, va="top", ha="right")
ax.plot([575], [out["threshold"]["Ceq_at_today_emission_realistic"]], "o", color=cm.INK, ms=6, mec="white",
        label="Today's emissions (575 Tg/yr)")
ax.set_ylim(0, 12000)
ax.set_xlabel("Total emissions E (Tg/yr)")
ax.set_ylabel("Equilibrium CH$_4$ (ppb)")
ax.set_title("Equilibrium methane vs emissions: the OH-runaway asymptote")
ax.legend(loc="upper left", fontsize=7.5)
cm.savefig(fig, "03_equilibrium_threshold")

# phase portrait with OH artificially slowed (tau_X = 0.5 yr) so the fast-slow structure is visible:
# trajectories first jump vertically onto the OH nullcline, then slide slowly along it to the node.
tx = 0.5
fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.4), sharey=True)
for j, (E, ttl) in enumerate([(E_now, "Emissions held at today's level"), (2.0 * E_now, "Emissions doubled")]):
    Ce = brentq(lambda C: E / beta - sink(C), 1, 1e6)
    cmax = max(4000, 1.25 * Ce)
    for c0 in np.linspace(900, cmax * 0.95, 5):
        for x0 in (0.45, 1.55):
            tr = solve_ivp(full_rhs, (0, 200), [c0, x0], args=(E, 1.0, tx), method="Radau", dense_output=True)
            ts = np.r_[np.linspace(0, 3, 60), np.linspace(3, 200, 300)]
            yy = tr.sol(ts)
            ax[j].plot(yy[0], yy[1], color=cm.MUTED, lw=0.7)
            k_ = 20
            ax[j].annotate("", xy=(yy[0][k_ + 1], yy[1][k_ + 1]), xytext=(yy[0][k_], yy[1][k_]),
                           arrowprops=dict(arrowstyle="->", color=cm.MUTED, lw=0.7))
    cc = np.linspace(800, cmax, 300)
    ax[j].plot(cc, (E / beta - ko * cc) / (k * cc), color=cm.SERIES[0], label="dC/dt = 0")
    ax[j].plot(cc, 1 / (a * cc / Cr + 1 - a), color=cm.SERIES[1], label="dX/dt = 0 (quasi-steady OH)")
    ax[j].plot(Ce, 1 / (a * Ce / Cr + 1 - a), "o", color=cm.INK, ms=6, mec="white", zorder=5,
               label=f"Stable node ({Ce:.0f} ppb)")
    ax[j].set_ylim(0.4, 1.6)
    ax[j].set_xlim(800, cmax)
    ax[j].set_xlabel("CH$_4$ (ppb)")
    ax[j].set_title(ttl)
    ax[j].legend(loc="upper right", fontsize=6.8)
ax[0].set_ylabel("OH (relative to 2010)")
cm.savefig(fig, "03_phase_portrait")
cm.save_json(out, "03_oh_feedback")

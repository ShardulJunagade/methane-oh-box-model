"""Step 2 - The one-box model dC/dt = E/beta - C/tau (L11 residence time, L02 confronting data).

* analytic solution for exponential emissions (as in the L11 two-box derivation)
* fit to NOAA 1984-2025 with the lifetime fixed (independent knowledge) and free
* cost surface: why lifetime and natural emissions cannot both be fitted (identifiability)
* top-down mass balance: total emissions implied by the observed record
* spin-up from 1750: does inventory + constant natural source reproduce the ice-core history?
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analysis import common as cm
from analysis import models as M

print("02 one-box model")
obs = cm.noaa_annual()
e_anth = cm.anthro_emissions("edgar")
drv = cm.drivers()
out = {}

# --- (a) analytic solution check (L11) ---------------------------------------------------
tau = 9.1
yrs = np.arange(1984, 2026)
q, alpha = np.polyfit(yrs - 1984, np.log(e_anth.loc[1984:2025].values), 1)[::-1]
q = np.exp(q)
E_nat_a = 200.0
C0 = M.noaa_jan1(1984)
t = np.linspace(0, 42, 421)
analytic = (C0 * np.exp(-t / tau) + (E_nat_a / cm.BETA) * tau * (1 - np.exp(-t / tau))
            + (q / cm.BETA) / (alpha + 1 / tau) * (np.exp(alpha * t) - np.exp(-t / tau)))
p_a = M.Params(k_oh=1 / tau - cm.K_OTHER, e_nat=E_nat_a)
expo = pd.Series(q * np.exp(alpha * (np.arange(1750, 2101) - 1984 + 0.5)), index=np.arange(1750, 2101))
num, jan1 = M.simulate(p_a, 1984, 2025, C0, expo, drv, steps_per_year=48)
err = np.abs(jan1 - np.interp(np.arange(43), t, analytic)).max()
out["analytic_check"] = {"q_Tg": q, "alpha_per_yr": alpha, "max_abs_diff_ppb": float(err)}
print(f"  analytic vs RK4 max difference {err:.3f} ppb (emission growth {100 * alpha:.2f} %/yr)")

# --- (b) fits ---------------------------------------------------------------------------
k_fix = M.K_OH_PRIOR[0]
fixed = M.fit_variant("one-box, lifetime fixed (9.1 yr)", ["e_nat"], {"a": 0.0, "k_oh": k_fix})
free = M.fit_variant("one-box, lifetime free", ["k_oh", "e_nat"], {"a": 0.0}, prior=False)
for f in (fixed, free):
    out[f.name] = f.summary()
    print(f"  {f.name}: tau={f.params.lifetime():.2f} yr, E_nat={f.params.e_nat:.0f} Tg, RMSE={f.rmse:.1f} ppb")

# out-of-sample: fit 1984-2005, predict 2006-2025
early = M.fit_variant("fixed, fit 1984-2005", ["e_nat"], {"a": 0.0, "k_oh": k_fix}, y1=2005)
pred, _ = M.simulate(early.params, 1984, 2025, C0, e_anth, drv)
out["holdout_fit1984_2005"] = {"rmse_2006_2025": float(np.sqrt(((pred.loc[2006:] - obs.loc[2006:2025]) ** 2).mean())),
                               "error_2025_ppb": float(pred.loc[2025] - obs.loc[2025])}

fig, ax = plt.subplots(2, 1, figsize=(7.2, 5.4), sharex=True)
ax[0].plot(obs.index, obs.values, "o", ms=3.5, color=cm.INK, label="NOAA annual mean", zorder=3)
ax[0].plot(fixed.model.index, fixed.model.values, color=cm.SERIES[0], label=f"Lifetime fixed 9.1 yr (RMSE {fixed.rmse:.0f} ppb)")
ax[0].plot(free.model.index, free.model.values, color=cm.SERIES[1],
           label=f"Lifetime free -> {free.params.lifetime():.1f} yr (RMSE {free.rmse:.0f} ppb)")
ax[0].plot(pred.index, pred.values, color=cm.SERIES[2], ls="--", label="Fit to 1984-2005 only, then predicted")
ax[0].set_ylabel("CH$_4$ (ppb)")
ax[0].set_title("One-box model vs observations (EDGAR emissions + constant natural source)")
ax[0].legend(loc="upper left")
g = obs.diff().dropna()
ax[1].bar(g.index, g.values, color=cm.GRID, edgecolor=cm.MUTED, lw=0.4, label="Observed growth")
ax[1].plot(fixed.model.diff().index, fixed.model.diff().values, color=cm.SERIES[0], label="Lifetime fixed")
ax[1].plot(free.model.diff().index, free.model.diff().values, color=cm.SERIES[1], label="Lifetime free")
ax[1].set_ylabel("Growth (ppb/yr)")
ax[1].legend(loc="upper left", ncols=3)
cm.savefig(fig, "02_one_box_fit")

# --- (c) cost surface -------------------------------------------------------------------
taus = np.linspace(5, 12, 36)
enats = np.linspace(50, 550, 41)
yrs = obs.loc[1984:2025]
J = np.empty((len(enats), len(taus)))
for i, en in enumerate(enats):
    for j, tt in enumerate(taus):
        p = M.Params(k_oh=1 / tt - cm.K_OTHER, e_nat=en)
        s, _ = M.simulate(p, 1984, 2025, C0, e_anth, drv, steps_per_year=4)
        J[i, j] = np.sqrt(((s.values - yrs.values) ** 2).mean())
fig, ax = plt.subplots(figsize=(5.6, 4.0))
cs = ax.contourf(taus, enats, np.log10(J), levels=20, cmap="Blues_r")
ax.contour(taus, enats, J, levels=[5, 10, 20, 40], colors=cm.INK2, linewidths=0.6)
ax.plot(free.params.lifetime(), free.params.e_nat, "o", color=cm.SERIES[1], ms=7, mec="white", label="Best fit (free)")
ax.axvspan(9.1 - 0.9, 9.1 + 0.9, color=cm.SERIES[3], alpha=0.25, lw=0, label="Literature lifetime (9.1 $\\pm$ 0.9 yr)")
ax.set_xlabel("Methane lifetime $\\tau$ (yr)")
ax.set_ylabel("Natural emissions (Tg/yr)")
ax.set_title("RMSE (ppb): a valley, not a point")
fig.colorbar(cs, ax=ax, label="log$_{10}$ RMSE (ppb)")
ax.legend(loc="lower left")
cm.savefig(fig, "02_cost_surface")

# --- (d) top-down mass balance ----------------------------------------------------------
C = obs.loc[1984:2025]
dCdt = (C.shift(-1) - C.shift(1)) / 2
E_td = (cm.BETA * (dCdt + C / 9.1)).dropna()
E_nat_implied = E_td - e_anth.reindex(E_td.index)
gcp_total = 575.0
out["topdown"] = {"tau_yr": 9.1,
                  "E_total_mean_2010_2019_Tg": float(E_td.loc[2010:2019].mean()),
                  "gcp_topdown_total_2010_2019_Tg": gcp_total,
                  "implied_nonEDGAR_mean_Tg": {"1985-1999": float(E_nat_implied.loc[1985:1999].mean()),
                                               "2000-2006": float(E_nat_implied.loc[2000:2006].mean()),
                                               "2007-2019": float(E_nat_implied.loc[2007:2019].mean()),
                                               "2020-2022": float(E_nat_implied.loc[2020:2022].mean()),
                                               "2023-2024": float(E_nat_implied.loc[2023:2024].mean())}}
fig, ax = plt.subplots(figsize=(7.2, 3.4))
ax.plot(E_td.index, E_td.values, color=cm.INK, label="Total emissions implied by NOAA record ($\\tau$ = 9.1 yr)")
ax.plot(e_anth.loc[1984:2025].index, e_anth.loc[1984:2025].values, color=cm.SERIES[0], label="EDGAR anthropogenic")
ax.plot(E_nat_implied.index, E_nat_implied.values, color=cm.SERIES[1], label="Difference = natural + missing sources")
ax.errorbar([2014.5], [gcp_total], yerr=[[gcp_total - 553], [586 - gcp_total]], fmt="s", color=cm.SERIES[2],
            ms=5, label="GCP top-down total 2010-2019")
ax.set_ylabel("Tg CH$_4$/yr")
ax.set_title("Top-down mass balance: E = $\\beta$(dC/dt + C/$\\tau$)")
ax.legend(loc="center left", fontsize=7.5)
cm.savefig(fig, "02_topdown_emissions")

# --- (e) spin-up from 1750 ---------------------------------------------------------------
hist = cm.ci_concentrations().CH4
tau_f = fixed.params.lifetime()
runs = {}
for label, en in [("Natural source = pre-industrial balance", cm.BETA * cm.C_PI / tau_f),
                  ("Natural source = 1984-2025 fit", fixed.params.e_nat)]:
    p = M.Params(k_oh=k_fix, e_nat=en)
    c0 = cm.C_PI if "pre" in label else (en + e_anth.loc[1750]) / cm.BETA * tau_f
    s, _ = M.simulate(p, 1750, 2025, c0, e_anth, drv, steps_per_year=4)
    runs[label] = s
    out.setdefault("spinup", {})[label] = {"E_nat": float(en), "C_1850": float(s.loc[1850]),
                                           "C_1950": float(s.loc[1950]), "C_2025": float(s.loc[2025])}
fig, ax = plt.subplots(figsize=(7.2, 3.4))
ax.plot(hist.loc[1850:].index, hist.loc[1850:].values, color=cm.INK, lw=2.2, label="Observed (ice cores + NOAA/AGAGE)")
for i, (k, s) in enumerate(runs.items()):
    ax.plot(s.index, s.values, color=cm.SERIES[i], label=k)
ax.set_xlim(1750, 2026)
ax.set_ylabel("CH$_4$ (ppb)")
ax.set_title("Spin-up from 1750 with a 9.1-yr lifetime")
ax.legend(loc="upper left")
cm.savefig(fig, "02_spinup_1750")
out["spinup_obs"] = {"C_1850": float(hist.loc[1850]), "C_1950": float(hist.loc[1950]), "C_2025": float(hist.loc[2025])}
cm.save_json(out, "02_one_box")

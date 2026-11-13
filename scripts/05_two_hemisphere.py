"""Step 5 - Interhemispheric two-box model (the L11 two-box model with exchange, applied to the atmosphere).

    dN/dt = 2 f E/beta - k N - (N - S)/tau_ex
    dS/dt = 2 (1-f) E/beta - k S + (N - S)/tau_ex          (each hemisphere = half the atmosphere)
    => gradient G = N - S obeys dG/dt = 2 (2f-1) E/beta - (k + 2/tau_ex) G

1. Methyl chloroform (emitted almost only in the NH) gives the exchange time tau_ex.
2. With tau_ex known, the CH4 gradient gives the NH share f of methane emissions, year by year.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from analysis import common as cm

print("05 two-hemisphere model")
out = {}


def two_box(N0, S0, E, f, k, tau_ex, beta, steps=12):
    """RK4, yearly-constant emissions; returns annual-mean N and S arrays."""
    h = 1 / steps
    N, S = N0, S0
    an, as_ = [], []

    def d(N, S, e):
        G = N - S
        return 2 * f * e / beta - k * N - G / tau_ex, 2 * (1 - f) * e / beta - k * S + G / tau_ex

    for e in E:
        vn, vs = [N], [S]
        for _ in range(steps):
            a1 = d(N, S, e)
            a2 = d(N + h / 2 * a1[0], S + h / 2 * a1[1], e)
            a3 = d(N + h / 2 * a2[0], S + h / 2 * a2[1], e)
            a4 = d(N + h * a3[0], S + h * a3[1], e)
            N += h / 6 * (a1[0] + 2 * a2[0] + 2 * a3[0] + a4[0])
            S += h / 6 * (a1[1] + 2 * a2[1] + 2 * a3[1] + a4[1])
            vn.append(N)
            vs.append(S)
        an.append(np.mean(vn))
        as_.append(np.mean(vs))
    return np.array(an), np.array(as_)


# --- 1. MCF: exchange time --------------------------------------------------------------------
mcf = cm.mcf_annual().loc[1993:2015]
E_m = cm.mcf_emissions().reindex(mcf.index).values


F_MCF = 0.98   # MCF was an industrial solvent used almost entirely in the NH; f is not identifiable
               # from 1993+ data (fit with f free hits the bound f = 1), so it is fixed. Sensitivity below.


def mcf_fit(f):
    def res(v):
        tau_ex, k = v
        n, s = two_box(mcf.nh.iloc[0], mcf.sh.iloc[0], E_m, f, k, tau_ex, cm.BETA_MCF)
        return np.r_[np.log(n) - np.log(mcf.nh.values), np.log(s) - np.log(mcf.sh.values)]
    r = least_squares(res, [1.0, 0.2], bounds=([0.2, 0.05], [5, 0.5]))
    cov = np.linalg.inv(r.jac.T @ r.jac) * (r.fun ** 2).sum() / (len(r.fun) - 2)
    return r, np.sqrt(np.diag(cov))


r, se = mcf_fit(F_MCF)
tau_ex, k_mcf = r.x
f_mcf = F_MCF
sens = {str(f): float(mcf_fit(f)[0].x[0]) for f in (0.95, 1.0)}
nM, sM = two_box(mcf.nh.iloc[0], mcf.sh.iloc[0], E_m, f_mcf, k_mcf, tau_ex, cm.BETA_MCF)
out["mcf_fit"] = {"tau_ex_yr": [tau_ex, se[0]], "NH_emission_share_fixed": f_mcf,
                  "tau_ex_if_NH_share": sens, "loss_rate_per_yr": [k_mcf, se[1]], "lifetime_yr": 1 / k_mcf,
                  "rms_log_error": float(np.sqrt((r.fun ** 2).mean()))}
print(f"  MCF: tau_ex = {tau_ex:.2f} +- {se[0]:.2f} yr (NH share fixed {f_mcf}), lifetime {1 / k_mcf:.2f} yr")

# --- 2. CH4 hemispheric means from the 23 sites ------------------------------------------------
z = cm.zonal_means(cm.site_annual())
z = z.loc[1985:2025]
noaa = cm.noaa_annual()
diff = (z.glob - noaa.reindex(z.index)).dropna()
out["site_global_minus_noaa_ppb"] = {"mean": float(diff.mean()), "max_abs": float(diff.abs().max())}
G = z.nh - z.sh

# top-down total emissions with tau = 9.15 yr, as in step 2
tau = 1 / (1 / 11.2 + cm.K_OTHER)
C = noaa.loc[1984:2025]
E_tot = (cm.BETA * ((C.shift(-1) - C.shift(1)) / 2 + C / tau)).dropna()
k = 1 / tau
dG = (G.shift(-1) - G.shift(1)) / 2
fN = (0.5 + cm.BETA * (dG + G * (k + 2 / tau_ex)) / (4 * E_tot)).dropna()
tau_lo, tau_hi = min(sens.values()) - 2 * se[0], max(sens.values()) + 2 * se[0]   # f range and +-2 s.e.
fN_low = (0.5 + cm.BETA * (dG + G * (k + 2 / tau_hi)) / (4 * E_tot)).dropna()
fN_high = (0.5 + cm.BETA * (dG + G * (k + 2 / tau_lo)) / (4 * E_tot)).dropna()

# EDGAR anthropogenic NH share (countries split by the fraction of their emissions south of the equator;
# approximate: whole-country weights, see RESULTS.md)
SH = {"AUS": 1, "NZL": 1, "ARG": 1, "CHL": 1, "URY": 1, "PRY": 1, "BOL": 1, "PER": 0.9, "BRA": 0.9, "ZAF": 1,
      "NAM": 1, "BWA": 1, "ZWE": 1, "ZMB": 1, "MOZ": 1, "MDG": 1, "MWI": 1, "AGO": 1, "LSO": 1, "SWZ": 1,
      "MUS": 1, "TZA": 1, "COD": 0.6, "COG": 0.6, "GAB": 0.7, "RWA": 1, "BDI": 1, "KEN": 0.5, "IDN": 0.7,
      "PNG": 1, "TLS": 1, "ECU": 0.5, "FJI": 1, "NCL": 1, "SLB": 1, "VUT": 1, "REU": 1, "COM": 1, "SYC": 1,
      "UGA": 0.2, "SOM": 0.3, "STP": 0.5, "WSM": 1, "TON": 1}
ct = cm.edgar_country_totals()
w_sh = pd.Series({c: SH.get(c, 0.0) for c in ct.index})
anth_sh = (ct.mul(w_sh, axis=0)).sum() / ct.sum()
f_anth = 1 - anth_sh
e_anth = cm.anthro_emissions("edgar")
yrs = fN.index.intersection(f_anth.index)
f_nat = (fN.loc[yrs] * E_tot.loc[yrs] - f_anth.loc[yrs] * e_anth.loc[yrs]) / (E_tot.loc[yrs] - e_anth.loc[yrs])

out["ch4"] = {"gradient_ppb": {str(y): float(G.loc[y]) for y in [1990, 2000, 2006, 2010, 2015, 2020, 2024] if y in G},
              "NH_share_total": {"1986-1999": float(fN.loc[1986:1999].mean()), "2000-2006": float(fN.loc[2000:2006].mean()),
                                 "2007-2019": float(fN.loc[2007:2019].mean()), "2020-2024": float(fN.loc[2020:].mean())},
              "NH_share_total_range_mean": [float(fN_low.mean()), float(fN_high.mean())], "tau_ex_range": [tau_lo, tau_hi],
              "NH_share_EDGAR_anthropogenic": {"1990": float(f_anth.loc[1990]), "2024": float(f_anth.loc[2024])},
              "NH_share_implied_natural": {"1990-1999": float(f_nat.loc[1990:1999].mean()),
                                           "2007-2019": float(f_nat.loc[2007:2019].mean()),
                                           "2020-2024": float(f_nat.loc[2020:].mean())},
              "gradient_trend_ppb_per_decade": float(10 * np.polyfit(G.index, G.values, 1)[0])}
print("  CH4 NH share of total emissions:", {k_: round(v, 3) for k_, v in out["ch4"]["NH_share_total"].items()})
cm.save_json(out, "05_two_hemisphere")

# --- figures ---------------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.3))
ax[0].semilogy(mcf.index, mcf.nh, "o", ms=3.5, color=cm.SERIES[0], label="NH observed")
ax[0].semilogy(mcf.index, mcf.sh, "o", ms=3.5, color=cm.SERIES[1], label="SH observed")
ax[0].semilogy(mcf.index, nM, color=cm.SERIES[0], lw=1.3, label="NH model")
ax[0].semilogy(mcf.index, sM, color=cm.SERIES[1], lw=1.3, label="SH model")
ax[0].set_ylabel("CH$_3$CCl$_3$ (ppt)")
ax[0].set_title(f"MCF fit: exchange time {tau_ex:.2f} yr")
ax[0].legend()
ax[1].plot(G.index, G.values, color=cm.INK, label="NH - SH CH$_4$ (sites)")
ax[1].set_ylabel("Gradient (ppb)")
ax[1].set_title("Methane north-south gradient")
ax[1].legend()
cm.savefig(fig, "05_mcf_gradient")

fig, ax = plt.subplots(figsize=(7.2, 3.3))
ax.fill_between(fN_low.index, fN_low.values, fN_high.values, color=cm.SERIES[0], alpha=0.18, lw=0,
                label="Range of exchange time")
ax.plot(fN.index, fN.values, color=cm.SERIES[0], label="All sources (from gradient)")
ax.plot(f_anth.index, f_anth.values, color=cm.SERIES[1], label="EDGAR anthropogenic")
ax.plot(f_nat.index, f_nat.values, color=cm.SERIES[2], label="Implied natural + missing")
ax.set_ylabel("Share emitted in the NH")
ax.set_ylim(0.3, 1.0)
ax.set_title("Where is methane emitted? NH share from the two-box model")
ax.legend(loc="lower left", ncols=2)
cm.savefig(fig, "05_nh_share")

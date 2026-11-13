"""Step 7 - Futures: SSP scenarios, the Global Methane Pledge, and the climate-methane feedback loop.

Model: best fit from step 4 (CH4-OH feedback + NOx/CO drivers + Arrhenius OH rate).
Future drivers: anthropogenic CH4, NOx and CO follow each RCMIP/SSP scenario as a ratio to its 2025
value, applied to our 2025 values; global temperature follows the IPCC AR6 SPM assessed warming of
each SSP (2021-40, 2041-60, 2081-2100 means).  Warming from methane is computed with the two-layer
EBM of step 6 (linear, so scenario differences can be superposed).

Global Methane Pledge (GMP): anthropogenic CH4 down 30 % from 2020 by 2030 (linear from 2025), then held.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analysis import common as cm
from analysis import models as M
from analysis import scenarios as SC

print("07 scenarios and the pledge")
bf = cm.load_json("04_best_fit")
ebm = cm.load_json("06_forcing_warming")["ebm_two_layer"]
P0 = M.Params(k_oh=bf["k_oh"], e_nat=bf["e_nat"], a=bf["a"], gN=bf["gN"], gC=bf["gC"], q10=bf["q10"],
              arrhenius=bf["arrhenius"])
Y0, Y1 = 1984, 2100
years = np.arange(Y0, Y1 + 1)
e_hist = cm.anthro_emissions("edgar")
drv_hist = cm.drivers()
obs = cm.noaa_annual()
N2O_2025 = float(cm.ci_concentrations().N2O.loc[2025])

E2020 = float(e_hist.loc[2020])
fut_years = SC.FUT
gmp = SC.gmp_path()
cases = {
    "SSP1-1.9": ("ssp119", None, None), "SSP1-2.6": ("ssp126", None, None), "SSP2-4.5 (baseline)": ("ssp245", None, None),
    "SSP3-7.0": ("ssp370", None, None), "SSP3-7.0-lowNTCF": ("ssp370", None, "ssp370-lowNTCF-gidden"),
    "SSP5-8.5": ("ssp585", None, None),
    "Global Methane Pledge (-30 % by 2030, then held)": ("ssp245", gmp, None),
    "Emissions frozen at 2025": ("ssp245", pd.Series(e_hist.loc[2025], index=fut_years), None),
    "All anthropogenic emissions stop in 2026": ("ssp245", pd.Series(0.0, index=fut_years), None),
}


two_layer = SC.two_layer


def run(p, ssp, override, ntcf):
    e, d = SC.scenario(ssp, override, ntcf)
    s, _ = M.simulate(p, Y0, Y1, M.noaa_jan1(Y0), e, d)
    return s


runs = {k: run(P0, *v) for k, v in cases.items()}
F = {k: M.erf_ch4_etminan(s.values, N2O_2025) for k, s in runs.items()}
base = "SSP2-4.5 (baseline)"
i25 = list(years).index(2025)
dT = {k: two_layer(F[k][i25:] - F[base][i25:], ebm["C_Wyr_m2K"], ebm["lambda_Wm2K"], ebm["gamma_Wm2K"])
      for k in runs}
out = {"best_fit": bf, "E2020_Tg": E2020, "scenarios": {}}
for k, s in runs.items():
    fut = s.loc[2025:]
    out["scenarios"][k] = {"C_2030": float(s.loc[2030]), "C_2050": float(s.loc[2050]), "C_2100": float(s.loc[2100]),
                           "C_peak": float(fut.max()), "peak_year": int(fut.idxmax()),
                           "dF_vs_baseline_2050": float(F[k][list(years).index(2050)] - F[base][list(years).index(2050)]),
                           "dT_vs_baseline_2050_K": float(dT[k][2050 - 2025]),
                           "dT_vs_baseline_2100_K": float(dT[k][2100 - 2025])}
    print(f"  {k:52s} C2050={s.loc[2050]:6.0f}  C2100={s.loc[2100]:6.0f}  dT2050={dT[k][25]:+.3f} K")

# e-folding after a stop
stop = runs["All anthropogenic emissions stop in 2026"].loc[2025:]
eq = stop.iloc[-1]
efold = float(np.interp(np.exp(-1), ((stop - eq) / (stop.iloc[0] - eq)).values[::-1], stop.index.values[::-1]) - 2025)
out["stop_case"] = {"C_2100": float(eq), "e_folding_yr": efold}

# GMP baseline comparison against MAGICC (RCMIP) concentrations
magicc = cm.rcmip_series("Atmospheric Concentrations|CH4", "ssp245", kind="concentrations")
out["magicc_ssp245"] = {"C_2050": float(magicc.loc[2050]), "C_2100": float(magicc.loc[2100])}

# --- Monte Carlo on the pledge benefit ----------------------------------------------------------
rng = np.random.default_rng(1)
mc = []
for _ in range(200):
    p = M.Params(**{k: getattr(P0, k) for k in ["k_oh", "e_nat", "a", "gN", "gC", "q10", "arrhenius"]})
    draw = rng.multivariate_normal([getattr(P0, k_) for k_ in bf["free"]], np.array(bf["cov"]))
    for k_, v_ in zip(bf["free"], draw):
        setattr(p, k_, max(v_, 1e-3) if k_ in ("k_oh", "e_nat") else v_)
    p.a = rng.uniform(0.25, 0.37)
    lam = rng.normal(ebm["lambda_Wm2K"], 0.25)
    b = run(p, "ssp245", None, None)
    gm = run(p, "ssp245", gmp, None)
    dF = M.erf_ch4_etminan(gm.values, N2O_2025) - M.erf_ch4_etminan(b.values, N2O_2025)
    t = two_layer(dF[i25:], ebm["C_Wyr_m2K"], max(lam, 0.6), ebm["gamma_Wm2K"])
    mc.append([gm.loc[2050] - b.loc[2050], t[25], t[75]])
mc = np.array(mc)
out["gmp_monte_carlo"] = {"n": len(mc), "dC_2050_p5_p50_p95": np.percentile(mc[:, 0], [5, 50, 95]).tolist(),
                          "dT_2050_p5_p50_p95": np.percentile(mc[:, 1], [5, 50, 95]).tolist(),
                          "dT_2100_p5_p50_p95": np.percentile(mc[:, 2], [5, 50, 95]).tolist()}
print("  GMP avoided warming 2050 (5/50/95 %):", np.round(out["gmp_monte_carlo"]["dT_2050_p5_p50_p95"], 3))

# --- climate-methane feedback loop gain (L04) ------------------------------------------------------
C_now = float(obs.loc[2025])
lam = ebm["lambda_Wm2K"]
dFdC = float((M.erf_ch4_etminan(C_now + 1, N2O_2025) - M.erf_ch4_etminan(C_now - 1, N2O_2025)) / 2)
tau_p = cm.load_json("03_oh_feedback")["lifetimes"]["2025"]["tau_perturbation"]
oh_sink_Tg = cm.BETA * C_now / cm.load_json("03_oh_feedback")["lifetimes"]["2025"]["tau_global"] * 0.82
dE_arrh = -oh_sink_Tg * cm.E_ARRH / cm.T_TROP ** 2
gain = {}
for q10 in (1.0, 2.0, 3.0, 4.0):
    dE = 165.0 * np.log(q10) / 10 + dE_arrh            # Tg/yr per K (GCP wetlands 165 Tg/yr)
    g = (1 / lam) * dFdC * (tau_p / cm.BETA) * dE
    gain[str(q10)] = {"dE_dT_Tg_per_K": dE, "loop_gain": g, "amplification": 1 / (1 - g)}
dE_needed = 1 / ((1 / lam) * dFdC * (tau_p / cm.BETA))
out["loop_gain"] = {"terms": {"dT/dF_K_per_Wm2": 1 / lam, "dF/dC_W_per_ppb": dFdC, "dCeq/dE_ppb_per_Tg": tau_p / cm.BETA,
                              "Arrhenius_dE_equiv_Tg_per_K": dE_arrh},
                    "by_q10": gain, "dE_dT_needed_for_runaway_Tg_per_K": dE_needed,
                    "wetland_q10_needed_for_runaway": float(np.exp(10 * (dE_needed - dE_arrh) / 165.0))}
print("  loop gain by Q10:", {k: round(v["loop_gain"], 4) for k, v in gain.items()},
      f"runaway needs dE/dT = {dE_needed:.0f} Tg/yr/K")

# warming-driven wetland emissions in SSP5-8.5 (forward run with Q10 = 2 and 3, no refit)
wet = {}
cm.WETLAND_SHARE = 165.0 / P0.e_nat   # Q10 acts on GCP wetlands (165 Tg/yr), not on the whole fitted natural term
for q10 in (1.0, 2.0, 3.0):
    p = M.Params(**{k: getattr(P0, k) for k in ["k_oh", "e_nat", "a", "gN", "gC", "arrhenius"]}, q10=q10)
    p.e_nat = P0.e_nat   # natural source at the reference temperature unchanged
    wet[str(q10)] = float(run(p, "ssp585", None, None).loc[2100])
cm.WETLAND_SHARE = 0.80
out["ssp585_wetland_feedback_C2100"] = wet
cm.save_json(out, "07_scenarios")

# --- figures ---------------------------------------------------------------------------------
show = ["SSP1-1.9", "SSP2-4.5 (baseline)", "SSP3-7.0", "SSP5-8.5", "Global Methane Pledge (-30 % by 2030, then held)",
        "All anthropogenic emissions stop in 2026"]
fig, ax = plt.subplots(figsize=(7.2, 4.0))
ax.plot(obs.index, obs.values, color=cm.INK, lw=2.2, label="Observed (NOAA)")
for i, k in enumerate(show):
    s = runs[k].loc[2025:]
    ax.plot(s.index, s.values, color=cm.SERIES[i], label=k)
ax.plot(magicc.loc[2025:2100].index, magicc.loc[2025:2100].values, color=cm.SERIES[1], ls="--", lw=1.0,
        label="SSP2-4.5, MAGICC7 (RCMIP)")
ax.set_xlim(1984, 2100)
ax.set_ylabel("CH$_4$ (ppb)")
ax.set_title("Methane to 2100: scenarios and the Global Methane Pledge")
ax.legend(loc="upper left", fontsize=7)
cm.savefig(fig, "07_scenarios")

fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.3))
for i, k in enumerate(["Global Methane Pledge (-30 % by 2030, then held)", "SSP1-1.9", "SSP3-7.0-lowNTCF",
                       "All anthropogenic emissions stop in 2026"]):
    ax[0].plot(np.arange(2025, 2101), dT[k] if k != "SSP3-7.0-lowNTCF" else
               two_layer(F[k][i25:] - F["SSP3-7.0"][i25:], ebm["C_Wyr_m2K"], ebm["lambda_Wm2K"], ebm["gamma_Wm2K"]),
               color=cm.SERIES[i], label=k.replace(" (-30 % by 2030, then held)", "") +
               (" vs SSP3-7.0" if "NTCF" in k else ""))
ax[0].axhline(0, color=cm.MUTED, lw=0.8)
ax[0].set_ylabel("Warming vs SSP2-4.5 (K)")
ax[0].set_title("Warming avoided by methane cuts")
ax[0].legend(fontsize=6.8)
q = list(gain)
ax[1].bar(range(len(q)), [gain[x]["loop_gain"] for x in q], color=cm.SERIES[0], width=0.6)
ax[1].set_xticks(range(len(q)), [f"Q10 = {float(x):g}" for x in q])
ax[1].set_ylabel("Loop gain g (runaway at g = 1)")
ax[1].set_title("Climate-wetland-methane loop")
cm.savefig(fig, "07_pledge_feedback")

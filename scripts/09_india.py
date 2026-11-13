"""Step 9 - India: where its methane comes from, intensity, and what a 30 % cut would do (L18).

Kaya-style identity (L18 production functions / intensity):  E = Pop x (GDP/Pop) x (E/GDP)
=> ln(E1/E0) = ln(Pop1/Pop0) + ln(y1/y0) + ln(I1/I0)  splits emission growth into population,
affluence and methane intensity of the economy.
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

print("09 India")
ct = cm.edgar_country_totals()
owid = pd.read_csv(cm.RAW / "owid/owid-co2-data.csv", low_memory=False)
countries = {"IND": "India", "CHN": "China", "USA": "United States", "BRA": "Brazil"}
world_e = cm.edgar_global()
out = {}

# --- sectors -----------------------------------------------------------------------------
sec = cm.edgar_sectors("India")
order = ["Livestock", "Rice", "Waste", "Coal", "Oil & gas", "Other (fuel use, burning, industry)"]
out["india_2024_Tg"] = sec.loc[2024].round(2).to_dict()
out["india_total_Tg"] = {"1970": float(sec.loc[1970].sum()), "2000": float(sec.loc[2000].sum()),
                         "2020": float(sec.loc[2020].sum()), "2024": float(sec.loc[2024].sum())}
out["india_share_of_world_pct"] = float(100 * sec.loc[2024].sum() / world_e.loc[2024])
out["india_agri_share_2024_pct"] = float(100 * sec.loc[2024, ["Livestock", "Rice"]].sum() / sec.loc[2024].sum())

# --- intensities --------------------------------------------------------------------------
rows = {}
for iso, name in countries.items():
    o = owid[owid.iso_code == iso].set_index("year")
    e = ct.loc[iso]
    for y in (1990, 2022):
        rows.setdefault(name, {})[y] = {"Tg": float(e.loc[y]), "kg_per_person": float(e.loc[y] * 1e9 / o.population.loc[y]),
                                        "g_per_usd": float(e.loc[y] * 1e12 / o.gdp.loc[y])}
ow = owid[owid.country == "World"].set_index("year")
for y in (1990, 2022):
    rows.setdefault("World", {})[y] = {"Tg": float(world_e.loc[y]), "kg_per_person": float(world_e.loc[y] * 1e9 / ow.population.loc[y]),
                                       "g_per_usd": float(world_e.loc[y] * 1e12 / ow.gdp.loc[y])}
out["intensity"] = rows

# --- Kaya decomposition for India 1990 -> 2022 ---------------------------------------------
o = owid[owid.iso_code == "IND"].set_index("year")
E0, E1 = ct.loc["IND", 1990], ct.loc["IND", 2022]
P0, P1 = o.population.loc[1990], o.population.loc[2022]
G0, G1 = o.gdp.loc[1990], o.gdp.loc[2022]
kaya = {"emissions": np.log(E1 / E0), "population": np.log(P1 / P0), "GDP per person": np.log(G1 / P1 / (G0 / P0)),
        "methane per GDP": np.log(E1 / G1 / (E0 / G0))}
out["kaya_log_changes_1990_2022"] = {k: float(v) for k, v in kaya.items()}
out["kaya_pct_changes_1990_2022"] = {k: float(100 * (np.exp(v) - 1)) for k, v in kaya.items()}
print("  India 1990->2022 (% change):", {k: round(v) for k, v in out["kaya_pct_changes_1990_2022"].items()})

# --- what if India cuts 30 % of its 2020 emissions by 2030? ----------------------------------
bf = cm.load_json("04_best_fit")
p = M.Params(k_oh=bf["k_oh"], e_nat=bf["e_nat"], a=bf["a"], gN=bf["gN"], gC=bf["gC"], q10=bf["q10"],
             arrhenius=bf["arrhenius"])
e_b, drv = SC.scenario("ssp245")
cut = 0.3 * ct.loc["IND", 2020]
ramp = np.clip((SC.FUT - 2025) / 5, 0, 1)
e_ind = e_b.copy()
e_ind.loc[2026:2100] = e_b.loc[2026:2100] - cut * ramp
base, _ = M.simulate(p, 1984, 2100, M.noaa_jan1(1984), e_b, drv)
ind, _ = M.simulate(p, 1984, 2100, M.noaa_jan1(1984), e_ind, drv)
ebm = cm.load_json("06_forcing_warming")["ebm_two_layer"]
N2O = float(cm.ci_concentrations().N2O.loc[2025])
dF = M.erf_ch4_etminan(ind.loc[2025:].values, N2O) - M.erf_ch4_etminan(base.loc[2025:].values, N2O)
dT = SC.two_layer(dF, ebm["C_Wyr_m2K"], ebm["lambda_Wm2K"], ebm["gamma_Wm2K"])
gmp_global_cut = 0.3 * world_e.loc[2020]
out["india_30pct_cut"] = {"cut_Tg": float(cut), "share_of_global_GMP_cut_pct": float(100 * cut / gmp_global_cut),
                          "dC_2050_ppb": float(ind.loc[2050] - base.loc[2050]), "dT_2050_K": float(dT[25]),
                          "dT_2100_K": float(dT[-1])}
print(f"  India -30 %: {cut:.1f} Tg/yr -> {ind.loc[2050] - base.loc[2050]:.0f} ppb, {dT[25]:+.3f} K in 2050")
cm.save_json(out, "09_india")

# --- figures --------------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.4))
ax[0].stackplot(sec.index, *[sec[c] for c in order], labels=order, colors=cm.SERIES[:6], edgecolor="white", linewidth=0.5)
ax[0].set_ylabel("Tg CH$_4$/yr")
ax[0].set_title("India's methane by source (EDGAR)")
ax[0].legend(loc="upper left", fontsize=6.5)
names = list(rows)
x = np.arange(len(names))
ax[1].bar(x - 0.18, [rows[n][1990]["kg_per_person"] for n in names], width=0.34, color=cm.SERIES[0], label="1990")
ax[1].bar(x + 0.18, [rows[n][2022]["kg_per_person"] for n in names], width=0.34, color=cm.SERIES[1], label="2022")
ax[1].set_xticks(x, names, fontsize=7.5)
ax[1].set_ylabel("kg CH$_4$ per person per year")
ax[1].set_title("Methane per person")
ax[1].legend()
cm.savefig(fig, "09_india")

"""Step 1 - What is in the data? Record, growth rate, inventories, budget, OH proxy, gradients.

Course link: L02 (quantify the object before modelling it), L11 (sources, sinks, residence time).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analysis import common as cm

print("01 explore data")
m = cm.noaa_monthly()
ann = cm.noaa_annual()
gr = cm.noaa_growth()
grj = cm.noaa_growth_jan_dec()

# --- Figure 1: the record and its growth rate ---------------------------------------------
fig, ax = plt.subplots(2, 1, figsize=(7.2, 5.6), sharex=True, gridspec_kw={"height_ratios": [1.2, 1]})
ax[0].plot(m.t, m.ch4, color=cm.MUTED, lw=0.8, label="Monthly mean")
ax[0].plot(m.t, m.trend, color=cm.INK, lw=1.8, label="Deseasonalised trend")
ax[0].set_ylabel("CH$_4$ (ppb)")
ax[0].set_title("Global methane, NOAA marine surface network")
ax[0].legend(loc="upper left")
# growth from the deseasonalised trend (12-month difference) and annual
tr = m.set_index("t").trend.dropna()
g12 = (tr - tr.shift(12)).dropna()
ax[1].axvspan(2000, 2007, color=cm.GRID, alpha=0.8, lw=0, zorder=0)
ax[1].axvspan(2020, 2022, color=cm.GRID, alpha=0.8, lw=0, zorder=0)
ax[1].bar(gr.index + 0.5, gr.values, width=0.8, color=cm.SERIES[0], label="Difference of annual means", zorder=2)
ax[1].plot(g12.index, g12.values, color=cm.INK, lw=1.2, label="12-month change of trend", zorder=3)
for x, txt in [(2003.5, "stall"), (2007.6, "restart"), (2021, "surge")]:
    ax[1].text(x, 18.6, txt, ha="center", fontsize=8, color=cm.INK2)
ax[1].axvline(2007, color=cm.MUTED, lw=0.8, ls=":", zorder=1)
ax[1].axhline(0, color=cm.MUTED, lw=0.8)
ax[1].set_ylabel("Growth (ppb/yr)")
ax[1].set_ylim(-5, 20)
ax[1].legend(loc="upper left", ncols=1, fontsize=7.5)
cm.savefig(fig, "01_record_growth")

# --- Figure 2: emission inventories ---------------------------------------------------------
edgar = cm.edgar_global()
ci = cm.ci_slcf().CH4.dropna()
hist = cm.rcmip_series("Emissions|CH4", "historical")
fig, ax = plt.subplots(figsize=(7.2, 3.4))
ax.plot(edgar.index, edgar.values, color=cm.SERIES[0], label="EDGAR 2025 (no wildfires)")
ax.plot(ci.index, ci.values, color=cm.SERIES[1], label="CEDS + fires (ClimateIndicator)")
ax.plot(hist.loc[1950:].index, hist.loc[1950:].values, color=cm.SERIES[2], label="RCMIP/CEDS historical (CMIP6)")
ax.set_ylabel("Tg CH$_4$/yr")
ax.set_title("Anthropogenic methane emissions: three inventories")
ax.legend(loc="upper left")
cm.savefig(fig, "01_inventories")

# --- Figure 3: sectors ----------------------------------------------------------------------
sec = cm.edgar_sectors()
order = ["Livestock", "Oil & gas", "Coal", "Waste", "Rice", "Other (fuel use, burning, industry)"]
fig, ax = plt.subplots(figsize=(7.2, 3.4))
ax.stackplot(sec.index, *[sec[c] for c in order], labels=order, colors=cm.SERIES[:6],
             edgecolor="white", linewidth=0.6)
ax.set_ylabel("Tg CH$_4$/yr")
ax.set_title("EDGAR anthropogenic methane by source")
ax.legend(loc="upper left", ncols=2)
cm.savefig(fig, "01_sectors")

# --- Figure 4: MCF and the latitudinal gradient -------------------------------------------
mcf = cm.mcf_monthly()
sites = cm.site_annual()
lat = cm.site_table().latitude
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.2))
ax[0].semilogy(mcf.year + (mcf.month - 0.5) / 12, mcf.glob, color=cm.INK, label="Global")
ax[0].semilogy(mcf.year + (mcf.month - 0.5) / 12, mcf.nh, color=cm.SERIES[0], lw=1.0, label="NH")
ax[0].semilogy(mcf.year + (mcf.month - 0.5) / 12, mcf.sh, color=cm.SERIES[1], lw=1.0, label="SH")
ax[0].set_ylabel("CH$_3$CCl$_3$ (ppt, log)")
ax[0].set_title("Methyl chloroform tracer (step 5)")
ax[0].legend()
for i, yr in enumerate([1990, 2005, 2020]):
    r = sites.loc[yr].dropna()
    ax[1].plot(lat.loc[r.index], r.values, "o", ms=4, color=cm.SERIES[i], label=str(yr),
               markeredgecolor="white", markeredgewidth=0.6)
ax[1].set_xlabel("Latitude")
ax[1].set_ylabel("CH$_4$ (ppb)")
ax[1].set_title("North-south gradient")
ax[1].legend()
cm.savefig(fig, "01_mcf_gradient")

# --- GCP budget summary -----------------------------------------------------------------
td = cm.gcp_budget("2010-2019 TD Budget").groupby("category").Global.mean()
bu = cm.gcp_budget("2010-2019 BU Budget").groupby("category").Global.mean()

summary = {
    "noaa_monthly_range": [f"{int(m.year.iloc[0])}-{int(m.month.iloc[0]):02d}", f"{int(m.year.iloc[-1])}-{int(m.month.iloc[-1]):02d}"],
    "ch4_1984_ppb": float(ann.loc[1984]), "ch4_2025_ppb": float(ann.loc[2025]),
    "ch4_1750_ppb": cm.C_PI, "increase_since_1750_pct": float(100 * (ann.loc[2025] / cm.C_PI - 1)),
    "growth_mean_ppb_yr": {"1984-1999": float(gr.loc[1985:1999].mean()), "2000-2006": float(gr.loc[2000:2006].mean()),
                           "2007-2013": float(gr.loc[2007:2013].mean()), "2014-2019": float(gr.loc[2014:2019].mean()),
                           "2020-2022": float(gr.loc[2020:2022].mean()), "2023-2025": float(gr.loc[2023:2025].mean())},
    "growth_top3_jan_dec": {int(k): float(v) for k, v in grj.nlargest(3).items()},
    "inventories_tg": {str(y): {"edgar": float(edgar.get(y, np.nan)), "ceds_fires": float(ci.get(y, np.nan)),
                                "rcmip": float(hist.get(y, np.nan))} for y in [1990, 2000, 2010, 2014, 2020, 2024]},
    "edgar_sectors_2024_tg": sec.loc[2024].round(1).to_dict(),
    "edgar_growth_2000_2006_tg_per_yr": float(np.polyfit(range(2000, 2007), edgar.loc[2000:2006], 1)[0]),
    "edgar_growth_2007_2019_tg_per_yr": float(np.polyfit(range(2007, 2020), edgar.loc[2007:2019], 1)[0]),
    "gcp_topdown_2010_2019": {k: float(td[k]) for k in ["Total sources", "Total anthropogenic sources",
                                                        "Total natural sources", "Wetlands", "Total chemical loss",
                                                        "Soil uptake", "Total sinks"]},
    "gcp_bottomup_2010_2019": {k: float(bu[k]) for k in ["Total anthropogenic sources", "Wetlands", "Freshwaters",
                                                         "Geological", "Biomass burning"]},
    "mcf_ppt": {"first": float(mcf.glob.dropna().iloc[0]), "last": float(mcf.glob.dropna().iloc[-1])},
    "ns_gradient_ppb": {str(y): float(sites.loc[y, ["alt", "brw"]].mean() - sites.loc[y, ["spo", "cgo"]].mean())
                        for y in [1990, 2000, 2010, 2020, 2024]},
}
cm.save_json(summary, "01_data_summary")
print("  growth means:", {k: round(v, 1) for k, v in summary["growth_mean_ppb_yr"].items()})

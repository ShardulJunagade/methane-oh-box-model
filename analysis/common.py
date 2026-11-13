"""Shared paths, constants, data loaders and plotting style for the methane project.

Every loader returns tidy pandas objects indexed by (integer) year unless stated.
All data come from data/raw/ (see data/README.md for sources and quirks).
"""
from __future__ import annotations

import json
import warnings
from functools import lru_cache
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=UserWarning, module="openpyxl")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
RES = ROOT / "results"
FIG = ROOT / "figures"
RES.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------------------------
# Physical constants and literature values (sources in RESULTS.md)
# ---------------------------------------------------------------------------------------------
BETA = 2.75            # Tg CH4 per ppb (global burden per ppb of surface mixing ratio; GCP convention)
BETA_MCF = 23.6        # Gg CH3CCl3 per ppt: 1.773e20 mol air x 1e-12 x 133.4 g/mol
C_PI = 729.2           # ppb CH4 in 1750 (AR6 / ClimateIndicator)
N2O_PI = 270.1         # ppb N2O in 1750
TAU_STRAT, TAU_SOIL, TAU_CL = 120.0, 150.0, 200.0   # yr: non-OH CH4 sinks (Prather et al. 2012)
K_OTHER = 1 / TAU_STRAT + 1 / TAU_SOIL + 1 / TAU_CL  # 0.0200 / yr
A_FEEDBACK = 0.31      # d ln OH / d ln CH4 = -0.31 (Holmes et al. 2013; IPCC AR6 Ch. 6)
E_ARRH = 1775.0        # K: activation temperature of CH4 + OH, k = 2.45e-12 exp(-1775/T) (JPL)
T_TROP = 272.0         # K: CH4-weighted tropospheric temperature
WETLAND_SHARE = 0.80   # wetlands / total natural sources, GCP 2010-2019 top-down (165/206)
REF_YEAR = 2010        # reference state for the OH model (X = 1 here)

# ---------------------------------------------------------------------------------------------
# Plot style (palette: validated categorical order, light mode; observations in ink)
# ---------------------------------------------------------------------------------------------
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 160, "font.size": 9, "axes.titlesize": 10,
    "axes.labelsize": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": INK2, "ytick.color": INK2, "lines.linewidth": 1.6, "legend.frameon": False,
    "legend.fontsize": 8, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def savefig(fig, name: str) -> None:
    fig.tight_layout()
    fig.savefig(FIG / f"{name}.png", bbox_inches="tight")
    plt.close(fig)
    print(f"  figure  figures/{name}.png")


def save_json(obj, name: str) -> None:
    def conv(o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))
    (RES / f"{name}.json").write_text(json.dumps(obj, indent=2, default=conv))
    print(f"  results results/{name}.json")


def load_json(name: str) -> dict:
    return json.loads((RES / f"{name}.json").read_text())


# ---------------------------------------------------------------------------------------------
# NOAA global methane
# ---------------------------------------------------------------------------------------------
def _noaa(path, names):
    return pd.read_csv(RAW / path, comment="#", sep=r"\s+", header=None, names=names)


@lru_cache
def noaa_monthly() -> pd.DataFrame:
    d = _noaa("noaa_gml/ch4_mm_gl.txt", ["year", "month", "t", "ch4", "ch4_unc", "trend", "trend_unc"])
    return d.replace(-9.99, np.nan)


@lru_cache
def noaa_annual() -> pd.Series:
    return _noaa("noaa_gml/ch4_annmean_gl.txt", ["year", "ch4", "unc"]).set_index("year")["ch4"]


@lru_cache
def noaa_growth() -> pd.Series:
    """Growth as difference of annual means (ppb/yr), labelled by the later year."""
    return noaa_annual().diff().dropna()


@lru_cache
def noaa_growth_jan_dec() -> pd.Series:
    return _noaa("noaa_gml/ch4_gr_gl.txt", ["year", "inc", "unc"]).set_index("year")["inc"]


# ---------------------------------------------------------------------------------------------
# Emissions
# ---------------------------------------------------------------------------------------------
@lru_cache
def edgar_table(sheet: str = "TOTALS BY COUNTRY") -> pd.DataFrame:
    """EDGAR sheet as a DataFrame. Reading the xlsx takes ~40 s, so a CSV copy is cached
    in data/processed/ (rebuilt automatically if the xlsx is newer)."""
    src = RAW / "edgar/EDGAR_CH4_1970_2024.xlsx"
    cache = ROOT / "data" / "processed" / f"edgar_{sheet.replace(' ', '_').lower()}.csv"
    if cache.exists() and cache.stat().st_mtime > src.stat().st_mtime:
        return pd.read_csv(cache, keep_default_na=False, na_values=[""])
    d = pd.read_excel(src, sheet_name=sheet, header=9)
    cache.parent.mkdir(exist_ok=True)
    d.to_csv(cache, index=False)
    return d


def _ycols(d):
    return [c for c in d.columns if str(c).startswith("Y_")]


@lru_cache
def edgar_global() -> pd.Series:
    """EDGAR 2025 global anthropogenic CH4 (Tg/yr), 1970-2024."""
    d = edgar_table()
    s = d[_ycols(d)].sum() / 1000
    s.index = [int(c[2:]) for c in s.index]
    return s


SECTOR_GROUPS = {  # EDGAR IPCC-2006 categories -> source process (for plots and India)
    "Enteric Fermentation": "Livestock", "Manure Management": "Livestock",
    "Rice cultivations": "Rice", "Solid Waste Disposal": "Waste", "Wastewater Treatment and Discharge": "Waste",
    "Oil and Natural Gas": "Oil & gas", "Solid Fuels": "Coal",
}


@lru_cache
def edgar_sectors(country: str | None = None) -> pd.DataFrame:
    """EDGAR CH4 by source group (Tg/yr), years as index; optionally one country by name."""
    d = edgar_table("IPCC 2006")
    if country:
        d = d[d.Name == country]
    name = d["ipcc_code_2006_for_standard_report_name"]
    grp = name.map(SECTOR_GROUPS).fillna("Other (fuel use, burning, industry)")
    g = d[_ycols(d)].groupby(grp.values).sum().T / 1000
    g.index = [int(c[2:]) for c in g.index]
    return g


@lru_cache
def edgar_country_totals() -> pd.DataFrame:
    d = edgar_table()
    t = d.set_index("Country_code_A3")[_ycols(d)] / 1000
    t.columns = [int(c[2:]) for c in t.columns]
    return t


@lru_cache
def ci_slcf() -> pd.DataFrame:
    return pd.read_csv(RAW / "climate_indicator/slcf_emissions.csv", index_col=0)


@lru_cache
def rcmip_emissions() -> pd.DataFrame:
    return pd.read_csv(RAW / "rcmip/rcmip_ch4_nox_co_voc_mcf_emissions_world.csv")


def rcmip_series(variable: str, scenario: str, kind: str = "emissions") -> pd.Series:
    d = rcmip_emissions() if kind == "emissions" else pd.read_csv(RAW / "rcmip/rcmip_ch4_mcf_concentrations_world.csv")
    r = d[(d.Variable == variable) & (d.Scenario == scenario)].iloc[0]
    s = pd.Series({int(c): float(r[c]) for c in d.columns if c.isdigit() and pd.notna(r[c])})
    full = pd.Series(np.interp(np.arange(s.index[0], s.index[-1] + 1), s.index, s.values),
                     index=np.arange(s.index[0], s.index[-1] + 1))
    return full


@lru_cache
def anthro_emissions(inventory: str = "edgar") -> pd.Series:
    """Global anthropogenic CH4 emissions 1750-2025 (Tg/yr).

    edgar: EDGAR 2025 for 1970-2024; 2025 = 2024 x (ClimateIndicator 2025/2024 ratio);
           1750-1969 = RCMIP/CEDS historical scaled to match EDGAR in 1970.
    ci:    ClimateIndicator (CEDS + biomass burning) 1970-2025, RCMIP scaled before 1970.
    """
    hist = rcmip_series("Emissions|CH4", "historical")
    if inventory == "edgar":
        main = edgar_global().copy()
        ci = ci_slcf().CH4
        main.loc[2025] = main.loc[2024] * ci.loc[2025] / ci.loc[2024]
    else:
        main = ci_slcf().CH4.dropna()
        main.index = main.index.astype(int)
    pre = hist.loc[:1969] * main.loc[1970] / hist.loc[1970]
    return pd.concat([pre, main]).sort_index()


@lru_cache
def drivers() -> pd.DataFrame:
    """Annual OH drivers 1750-2025: NOx and CO emissions normalised to REF_YEAR, and temperature anomaly.

    dT is the HadCRUT5 anomaly minus its 2001-2020 mean (the reference state of the OH model),
    smoothed with a 3-yr running mean; before 1850 it is set to the 1850-1879 mean.
    """
    s = ci_slcf().copy()
    s.index = s.index.astype(int)
    out = pd.DataFrame({"nox": s.NOx / s.NOx.loc[REF_YEAR], "co": s.CO / s.CO.loc[REF_YEAR]})
    t = hadcrut()
    t = t.rolling(3, center=True, min_periods=1).mean()
    dT = t - t.loc[2001:2020].mean()
    out["dT"] = dT.reindex(out.index)
    out.loc[:1849, "dT"] = dT.loc[1850:1879].mean()
    return out


# ---------------------------------------------------------------------------------------------
# Other observations
# ---------------------------------------------------------------------------------------------
@lru_cache
def hadcrut(include_partial: bool = False) -> pd.Series:
    d = pd.read_csv(RAW / "hadcrut/HadCRUT.5.2.0.0.analysis.summary_series.global.annual.csv")
    s = d.set_index("Time")["Anomaly (deg C)"]
    return s if include_partial else s.loc[:2025]


@lru_cache
def ci_concentrations() -> pd.DataFrame:
    return pd.read_csv(RAW / "climate_indicator/ghg_concentrations.csv").set_index("YYYY")


@lru_cache
def ci_erf() -> pd.DataFrame:
    return pd.read_csv(RAW / "climate_indicator/ERF_best.csv", index_col=0)


@lru_cache
def ci_erf_aggregates() -> pd.DataFrame:
    d = pd.read_csv(RAW / "climate_indicator/ERF_best_aggregates.csv", index_col=0)
    d.index = (d.index - 0.5).astype(int)
    return d


@lru_cache
def mcf_monthly() -> pd.DataFrame:
    lines = (RAW / "noaa_hats/HATS_global_MC.txt").read_text(encoding="latin-1").splitlines()
    i = next(k for k, s in enumerate(lines) if s.startswith("HATS_MC_YYYY"))
    from io import StringIO
    d = pd.read_csv(StringIO("\n".join(lines[i:])), sep=r"\s+", na_values="nan", dtype={"HATS_MC_Programs": str})
    return d.rename(columns={"HATS_MC_YYYY": "year", "HATS_MC_MM": "month", "HATS_NH_MC": "nh",
                             "HATS_SH_MC": "sh", "HATS_Global_MC": "glob"})


@lru_cache
def mcf_annual() -> pd.DataFrame:
    """Annual means of complete years only (12 months with a global value)."""
    d = mcf_monthly().dropna(subset=["glob"])
    g = d.groupby("year")
    a = g[["nh", "sh", "glob"]].mean()
    return a[g.size() == 12]


@lru_cache
def mcf_emissions() -> pd.Series:
    e = pd.read_csv(RAW / "climate_indicator/ghg_equivalent_emissions.csv", index_col=0)["CH3CCl3"]
    e.index = (e.index - 0.5).astype(int)
    return e


@lru_cache
def gcp_budget(sheet: str = "2010-2019 TD Budget") -> pd.DataFrame:
    f = RAW / "gcp/Global_methane_budget_2023_2000_2020_v1.xlsx"
    raw = pd.read_excel(f, sheet_name=sheet, header=None)
    h = next(i for i in range(80) if str(raw.iloc[i, 0]).strip().lower() == "category")
    d = pd.read_excel(f, sheet_name=sheet, header=h)
    d["Global"] = pd.to_numeric(d["Global"], errors="coerce")
    return d


# ---------------------------------------------------------------------------------------------
# NOAA site files -> hemispheric / global means
# ---------------------------------------------------------------------------------------------
@lru_cache
def site_table() -> pd.DataFrame:
    return pd.read_csv(RAW / "noaa_gml/sites_ch4/site_table.csv").set_index("site")


def site_annual(min_months: int = 9) -> pd.DataFrame:
    """Annual means per site (columns = sites), only years with >= min_months months."""
    folder = RAW / "noaa_gml/sites_ch4"
    cols = {}
    for p in folder.glob("*_month.txt"):
        d = _noaa(p.relative_to(RAW), ["site", "year", "month", "value"])
        g = d.groupby("year").value
        a = g.mean()[g.count() >= min_months]
        cols[d.site.iloc[0].lower()] = a
    return pd.DataFrame(cols).sort_index()


def zonal_means(annual: pd.DataFrame, exclude=("mlo", "izo"), min_nh=5, min_sh=4) -> pd.DataFrame:
    """Area-weighted NH, SH and global means by linear interpolation in sin(latitude).

    For each year the site values are interpolated in x = sin(lat) (equal-area coordinate),
    held flat beyond the most poleward site, and integrated over each hemisphere.
    """
    lat = site_table().latitude
    x = np.linspace(-1, 1, 401)
    rows = {}
    for yr, r in annual.drop(columns=[c for c in exclude if c in annual]).iterrows():
        r = r.dropna()
        xs = np.sin(np.radians(lat.loc[r.index].values))
        if (xs > 0).sum() < min_nh or (xs < 0).sum() < min_sh:
            continue
        o = np.argsort(xs)
        y = np.interp(x, xs[o], r.values[o])
        nh, sh = y[x >= 0].mean(), y[x <= 0].mean()
        rows[yr] = {"nh": nh, "sh": sh, "glob": 0.5 * (nh + sh), "n_sites": len(r)}
    return pd.DataFrame(rows).T

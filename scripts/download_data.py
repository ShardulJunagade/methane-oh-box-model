"""Download every dataset for the methane-OH box-model project into project_methane/data/raw/.

Run from the project folder:

    uv run scripts/download_data.py            # download what is missing, then check everything
    uv run scripts/download_data.py --force    # re-download everything
    uv run scripts/download_data.py --only noaa_global hats_mcf
    uv run scripts/download_data.py --list     # show the source keys

Every file is checked after it arrives: it is parsed, its time coverage and value ranges are
tested, and overlapping datasets are compared with each other. A failed check is reported and
the script exits with status 1 at the end. Results are written to
  data/MANIFEST.csv      one row per file: source URL, size, SHA-256, date, check summary
  data/CHECK_REPORT.txt  the full check output, including notes on odd things found
All sources are public. None needs a login. The GCP methane budget sits behind a licence
page on the ICOS portal; the script accepts the licence (CC BY 4.0) with the cookie that
the portal itself sets when you click "accept".
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import io
import shutil
import sys
import time
import urllib.error
import urllib.request
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RAW = DATA / "raw"
UA = "EH605-methane-course-project/1.0 (python-urllib; IIT Gandhinagar coursework)"
TODAY = dt.date.today()

GML = "https://gml.noaa.gov"
TRENDS = f"{GML}/webdata/ccgg/trends/ch4"
SITE_TXT = f"{GML}/aftp/data/trace_gases/ch4/flask/surface/txt"
CI_TAG = "v6.4.2"  # ClimateIndicator/forcing-timeseries release, 2026-06-15 (pinned for repeatability)
CI = f"https://raw.githubusercontent.com/ClimateIndicator/forcing-timeseries/{CI_TAG}"
EDGAR = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/EDGAR_2025_GHG"
RCMIP = "https://rcmip-protocols-au.s3-ap-southeast-2.amazonaws.com/v5.1.0"
HADCRUT = ("https://www.metoffice.gov.uk/hadobs/hadcrut5/data/HadCRUT.5.2.0.0/analysis/diagnostics/"
           "HadCRUT.5.2.0.0.analysis.summary_series.global.annual.csv")
GCP_ID = "-MqGCn38zUlEi4_aQm-1_w2I"  # ICOS object: Global_methane_budget_2023_2000_2020_v1.xlsx

# Background (clean-air) NOAA flask sites from the Arctic to the South Pole.
# Latitude/longitude/altitude are read from NOAA's own event files, not typed in here.
CH4_SITES = ["alt", "zep", "brw", "ice", "mhd", "azr", "izo", "mid", "kum", "mlo", "gmi", "rpb",
             "chr", "sey", "asc", "smo", "eic", "cgo", "bhd", "crz", "psa", "syo", "spo"]

RCMIP_SCENARIOS = ["historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp370-lowNTCF-gidden",
                   "ssp434", "ssp460", "ssp534-over", "ssp585"]
RCMIP_EMIS_VARS = ["Emissions|CH4", "Emissions|CH4|MAGICC AFOLU", "Emissions|CH4|MAGICC Fossil and Industrial",
                   "Emissions|NOx", "Emissions|CO", "Emissions|VOC", "Emissions|Montreal Gases|CH3CCl3"]
RCMIP_CONC_VARS = ["Atmospheric Concentrations|CH4", "Atmospheric Concentrations|Montreal Gases|CH3CCl3"]


# --------------------------------------------------------------------------------------------
# HTTP helpers
# --------------------------------------------------------------------------------------------
def _request(url: str, cookie: str | None = None, byte_range: tuple[int, int] | None = None):
    headers = {"User-Agent": UA}
    if cookie:
        headers["Cookie"] = cookie
    if byte_range:
        headers["Range"] = f"bytes={byte_range[0]}-{byte_range[1]}"
    return urllib.request.Request(url, headers=headers)


def get_bytes(url: str, cookie=None, byte_range=None, tries=3) -> bytes:
    for k in range(tries):
        try:
            with urllib.request.urlopen(_request(url, cookie, byte_range), timeout=300) as r:
                return r.read()
        except (urllib.error.URLError, TimeoutError) as e:
            if k == tries - 1:
                raise
            print(f"      retry {k + 1} after error: {e}")
            time.sleep(5 * (k + 1))
    raise RuntimeError("unreachable")


def download(url: str, dest: Path, cookie=None, tries=3) -> None:
    """Stream url to dest via a .part file; fail if the byte count differs from Content-Length."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    for k in range(tries):
        try:
            with urllib.request.urlopen(_request(url, cookie), timeout=600) as r, open(part, "wb") as f:
                expected = r.headers.get("Content-Length")
                n = 0
                while chunk := r.read(1 << 20):
                    f.write(chunk)
                    n += len(chunk)
            if expected is not None and int(expected) != n:
                raise IOError(f"truncated: got {n} of {expected} bytes")
            head = part.read_bytes()[:200].lstrip().lower()
            if head.startswith((b"<!doctype html", b"<html")) and dest.suffix not in (".html", ".htm"):
                raise IOError("server returned an HTML page instead of data")
            part.replace(dest)
            return
        except (urllib.error.URLError, TimeoutError, IOError) as e:
            part.unlink(missing_ok=True)
            if k == tries - 1:
                raise
            print(f"      retry {k + 1} after error: {e}")
            time.sleep(5 * (k + 1))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------------------------
class CheckError(Exception):
    pass


def need(cond, msg: str) -> None:
    if not cond:
        raise CheckError(msg)


@dataclass
class Result:
    summary: str
    notes: list[str] = field(default_factory=list)


def read_noaa(path: Path, names: list[str]) -> pd.DataFrame:
    return pd.read_csv(path, comment="#", sep=r"\s+", header=None, names=names)


def months_contiguous(year: pd.Series, month: pd.Series) -> bool:
    t = year.to_numpy() * 12 + month.to_numpy()
    return bool((np.diff(t) == 1).all())


def check_ch4_mm(p: Path) -> Result:
    d = read_noaa(p, ["year", "month", "decimal", "average", "average_unc", "trend", "trend_unc"])
    need(len(d) > 480, f"only {len(d)} months")
    need(months_contiguous(d.year, d.month), "gap or duplicate in the monthly series")
    need(d.average.between(1600, 2100).all(), "CH4 outside 1600-2100 ppb")
    need(d.trend.between(1600, 2100).all(), "trend outside 1600-2100 ppb")
    # -9.99 is NOAA's missing-value code; it appears in the uncertainty columns of the newest months.
    miss = d[(d.average_unc == -9.99) | (d.trend_unc == -9.99)]
    unc = d.average_unc.drop(miss.index)
    need((unc >= 0).all() and (unc < 10).all(), "odd uncertainties")
    need(miss.empty or miss.index.min() > len(d) - 13, "missing uncertainties outside the newest months")
    jump = d.average.diff().abs().max()
    need(jump < 20, f"month-to-month jump of {jump:.1f} ppb")
    need(d.year.iloc[-1] >= TODAY.year - 1, "record does not reach last year")
    notes = [f"Seasonal cycle amplitude (average - trend) is {(d.average - d.trend).abs().max():.1f} ppb at most; "
             "use `trend` (deseasonalised) for growth rates."]
    if len(miss):
        notes.append(f"The newest {len(miss)} months ({miss.year.iloc[0]}-{miss.month.iloc[0]:02d} .. "
                     f"{miss.year.iloc[-1]}-{miss.month.iloc[-1]:02d}) have uncertainty -9.99 (NOAA's missing-value "
                     "code): they are preliminary and will be revised. Treat -9.99 as NaN.")
    return Result(f"{len(d)} months {d.year.iloc[0]}-{d.month.iloc[0]:02d} .. {d.year.iloc[-1]}-{d.month.iloc[-1]:02d}; "
                  f"{d.average.min():.1f}..{d.average.max():.1f} ppb", notes)


def check_ch4_ann(p: Path) -> Result:
    d = read_noaa(p, ["year", "mean", "unc"])
    need((np.diff(d.year) == 1).all(), "missing years")
    need(d["mean"].between(1600, 2100).all(), "annual mean outside 1600-2100 ppb")
    m = read_noaa(RAW / "noaa_gml/ch4_mm_gl.txt", ["year", "month", "decimal", "average", "u", "trend", "tu"])
    full = m.groupby("year").filter(lambda g: len(g) == 12).groupby("year").average.mean()
    j = d.set_index("year")["mean"].reindex(full.index).dropna()
    diff = (j - full.loc[j.index]).abs().max()
    need(diff < 1.0, f"annual means differ from mean of monthly values by {diff:.2f} ppb")
    return Result(f"{len(d)} years {d.year.min()}..{d.year.max()}; {d['mean'].iloc[0]:.1f} -> {d['mean'].iloc[-1]:.1f} ppb; "
                  f"matches monthly file within {diff:.2f} ppb")


def check_ch4_gr(p: Path) -> Result:
    d = read_noaa(p, ["year", "ann_inc", "unc"])
    need((np.diff(d.year) == 1).all(), "missing years")
    need(d.ann_inc.between(-5, 25).all(), "growth rate outside -5..25 ppb/yr")
    top = d.nlargest(3, "ann_inc")
    low = d[(d.year >= 2000) & (d.year <= 2006)].ann_inc.mean()
    notes = [f"Stall: mean growth 2000-2006 = {low:.1f} ppb/yr. Largest growth: "
             + ", ".join(f"{int(r.year)} ({r.ann_inc:.1f})" for r in top.itertuples()) + " ppb/yr."]
    return Result(f"{len(d)} years {d.year.min()}..{d.year.max()}; {d.ann_inc.min():.1f}..{d.ann_inc.max():.1f} ppb/yr", notes)


def site_meta(species_url: str, prefix: str, code: str, suffix: str) -> dict:
    """Read latitude/longitude/altitude from the first ~40 kB of a NOAA event file (HTTP Range)."""
    txt = get_bytes(f"{species_url}/{prefix}_{code}_{suffix}_event.txt", byte_range=(0, 40000)).decode("latin-1")
    # Two header styles exist: older files have "# data_fields: site year ...", newer ones an
    # uncommented line of column names starting with "site_code".
    fields = None
    for line in txt.splitlines():
        if line.startswith("# data_fields:"):
            fields = line.split(":", 1)[1].split()
        elif line.startswith("site_code "):
            fields = line.split()
        elif fields and line and not line.startswith("#"):
            v = line.split()
            if len(v) >= len(fields) - 2:
                rec = dict(zip(fields, v))
                return {"latitude": float(rec["latitude"]), "longitude": float(rec["longitude"]),
                        "altitude_m": float(rec["altitude"])}
    raise CheckError(f"could not read coordinates for site {code}")


def check_site_folder(folder: Path) -> tuple:
    kind, lo, hi = "ch4", 1500, 2300
    rows, notes = [], []
    for p in sorted(folder.glob("*_month.txt")):
        d = read_noaa(p, ["site", "year", "month", "value"])
        code = d.site.iloc[0].lower()
        need(len(d) > 24, f"{p.name}: only {len(d)} months")
        need(d.value.between(lo, hi).all(), f"{p.name}: values outside {lo}..{hi}")
        gaps = int((d.year * 12 + d.month).diff().gt(1).sum())
        rows.append({"site": code, "first": f"{d.year.iloc[0]}-{d.month.iloc[0]:02d}",
                     "last": f"{d.year.iloc[-1]}-{d.month.iloc[-1]:02d}", "n_months": len(d), "gaps": gaps})
    t = pd.DataFrame(rows)
    need(len(t) == len(CH4_SITES), f"expected all {kind} sites, found {len(t)}")
    last = t["last"].max()
    stale = t[t["last"] < f"{int(last[:4]) - 1}-{last[5:]}"]
    if len(stale):
        notes.append(f"{kind} sites whose record ends more than a year before the newest: "
                     + ", ".join(f"{r.site} (ends {r.last})" for r in stale.itertuples()))
    gappy = t[t.gaps > 6]
    if len(gappy):
        notes.append(f"{kind} sites with more than 6 missing-month gaps: "
                     + ", ".join(f"{r.site} ({r.gaps})" for r in gappy.itertuples()))
    return Result(f"{len(t)} sites; records {t['first'].min()} .. {last}", notes), t


def check_ch4_sites(p: Path) -> Result:
    res, t = check_site_folder(p)
    meta = pd.read_csv(p / "site_table.csv")
    need(set(meta.site) >= set(t.site), "site_table.csv is missing sites")
    need(meta.latitude.between(-90, 90).all(), "bad latitude")
    nh, sh = (meta.latitude > 0).sum(), (meta.latitude < 0).sum()
    high = meta[meta.altitude_m > 1500].site.tolist()
    res.notes.append(f"{nh} NH and {sh} SH sites. Mountain/ice-sheet sites (>1500 m): {', '.join(high)}; "
                     "they sample the free troposphere rather than the marine boundary layer.")
    res.summary += f"; coordinates in site_table.csv ({nh} NH, {sh} SH)"
    return res


def check_hats_mcf(p: Path) -> Result:
    lines = p.read_text(encoding="latin-1").splitlines()
    i = next(k for k, s in enumerate(lines) if s.startswith("HATS_MC_YYYY"))
    d = pd.read_csv(io.StringIO("\n".join(lines[i:])), sep=r"\s+", na_values="nan",
                    dtype={"HATS_MC_Programs": str})
    g = d.dropna(subset=["HATS_Global_MC"])
    need(len(g) > 300, f"only {len(g)} months with a global mean")
    need(months_contiguous(g.HATS_MC_YYYY, g.HATS_MC_MM), "gap in the global monthly series")
    need(g.HATS_Global_MC.between(0, 200).all(), "global MCF outside 0..200 ppt")
    need(g.HATS_Global_MC.iloc[-1] < g.HATS_Global_MC.iloc[0] / 20, "MCF has not fallen as expected")
    yrs = g.groupby("HATS_MC_YYYY").HATS_Global_MC.mean()
    post = yrs.loc[2005:2020]
    e_fold = -1 / np.polyfit(post.index, np.log(post.values), 1)[0]
    notes = [f"Global MCF falls from {g.HATS_Global_MC.iloc[0]:.1f} ppt ({g.HATS_MC_YYYY.iloc[0]}-{g.HATS_MC_MM.iloc[0]:02d}) to "
             f"{g.HATS_Global_MC.iloc[-1]:.2f} ppt ({g.HATS_MC_YYYY.iloc[-1]}-{g.HATS_MC_MM.iloc[-1]:02d}). Quick check: "
             f"the 2005-2020 e-folding time is {e_fold:.1f} yr (total lifetime incl. ocean/stratosphere loss and "
             "small residual emissions).",
             "Rows before 1992-05 have NH data only (SH and global columns are NaN)."]
    return Result(f"{len(g)} global months {g.HATS_MC_YYYY.iloc[0]}-{g.HATS_MC_MM.iloc[0]:02d} .. "
                  f"{g.HATS_MC_YYYY.iloc[-1]}-{g.HATS_MC_MM.iloc[-1]:02d}", notes)


def edgar_totals(xlsx: Path) -> pd.DataFrame:
    """Return the country-totals sheet of an EDGAR workbook (values in Gg CH4)."""
    xl = pd.ExcelFile(xlsx)
    sheet = next(s for s in xl.sheet_names if "TOTALS BY COUNTRY" in s.upper())
    raw = pd.read_excel(xlsx, sheet_name=sheet, header=None, nrows=30)
    hdr = next(i for i in range(30) if "Y_1970" in raw.iloc[i].astype(str).tolist())
    return pd.read_excel(xlsx, sheet_name=sheet, header=hdr)


def check_edgar(p: Path) -> Result:
    xlsx = sorted(p.glob("*.xlsx"))
    need(xlsx, "no xlsx inside the EDGAR zip")
    d = edgar_totals(xlsx[0])
    ycols = [c for c in d.columns if str(c).startswith("Y_")]
    need(ycols[0] == "Y_1970" and len(ycols) >= 55, f"year columns {ycols[0]}..{ycols[-1]}")
    tot = d[ycols].sum() / 1000  # Gg -> Tg
    need(tot.between(200, 500).all(), "global total outside 200..500 Tg CH4/yr")
    name_col = next(c for c in d.columns if str(c).lower() in ("name", "country"))
    ind = d[d[name_col].astype(str).str.fullmatch("India")][ycols].sum() / 1000
    need(ind.iloc[-1] > 10, "India total looks wrong")
    xl = pd.ExcelFile(xlsx[0])
    notes = [f"Sheets in {xlsx[0].name}: {', '.join(xl.sheet_names)}",
             f"Global anthropogenic total {tot.iloc[0]:.0f} Tg (1970) -> {tot.iloc[-1]:.0f} Tg ({ycols[-1][2:]}); "
             f"India {ind.iloc[0]:.1f} -> {ind.iloc[-1]:.1f} Tg ({100 * ind.iloc[-1] / tot.iloc[-1]:.1f}% of world)."]
    return Result(f"{len(d)} countries x {len(ycols)} years ({ycols[0][2:]}..{ycols[-1][2:]}), Gg CH4", notes)


def check_owid(p: Path) -> Result:
    d = pd.read_csv(p, low_memory=False)
    for c in ("country", "year", "iso_code", "methane", "population", "gdp"):
        need(c in d.columns, f"column {c} missing")
    w = d[d.country == "World"].dropna(subset=["methane"])
    i = d[d.country == "India"].dropna(subset=["methane"])
    need(len(w) > 20 and len(i) > 20, "methane series too short")
    need(w.methane.between(500, 15000).all(), "World methane (MtCO2e) outside 500..15000")
    need(w.methane.iloc[-1] > 7000, "recent World methane (MtCO2e) too low")
    code = pd.read_csv(p.with_name("owid-co2-codebook.csv"))
    desc = code.set_index("column").loc["methane", "description"]
    notes = [f"`methane` is in million tonnes CO2-equivalent (GWP100), not Tg CH4. Codebook: {desc[:160]}...",
             f"World methane present {int(w.year.min())}..{int(w.year.max())}; India {int(i.year.min())}..{int(i.year.max())}."]
    return Result(f"{len(d)} rows, {d.country.nunique()} countries/regions, years {d.year.min()}..{d.year.max()}", notes)


def check_gcp(p: Path) -> Result:
    xl = pd.ExcelFile(p)
    need(len(xl.sheet_names) >= 8, "fewer sheets than expected")
    for want in ("CH4_observation", "TD Budget", "BU Budget", "Anthropogenic"):
        need(any(want in s for s in xl.sheet_names), f"no sheet containing '{want}'")
    return Result(f"{len(xl.sheet_names)} sheets: {'; '.join(xl.sheet_names)}",
                  ["Budget tables are decadal means (2000-2009, 2010-2019) plus 2020, in Tg CH4/yr; "
                   "each sheet has a free-text header block above the table (read with skiprows)."])


def check_ci_conc(p: Path) -> Result:
    d = pd.read_csv(p)
    need({"YYYY", "CO2", "CH4", "N2O", "CH3CCl3"} <= set(d.columns), "columns missing")
    need(d.YYYY.iloc[0] == 1750 and d.YYYY.iloc[1] == 1850, "expected 1750 then 1850 onwards")
    need((np.diff(d.YYYY.iloc[1:]) == 1).all(), "gap after 1850")
    need(d.CH4.between(600, 2100).all(), "CH4 outside 600..2100 ppb")
    noaa = read_noaa(RAW / "noaa_gml/ch4_annmean_gl.txt", ["year", "mean", "unc"]).set_index("year")["mean"]
    j = d.set_index("YYYY").CH4.reindex(noaa.index).dropna()
    diff = (j - noaa.loc[j.index])
    notes = [f"CH4 {d.CH4.iloc[0]:.1f} ppb (1750), {d.CH4.iloc[1]:.1f} (1850), {d.CH4.iloc[-1]:.1f} ({d.YYYY.iloc[-1]}). "
             "There are no rows for 1751-1849.",
             f"Difference from the NOAA global annual mean over {j.index.min()}-{j.index.max()}: "
             f"{diff.min():+.1f}..{diff.max():+.1f} ppb (this file blends NOAA and AGAGE after 1979, ice cores before)."]
    need(diff.abs().max() < 15, f"differs from NOAA by up to {diff.abs().max():.1f} ppb")
    return Result(f"{len(d)} rows 1750, 1850..{d.YYYY.iloc[-1]}; {len(d.columns) - 1} gases", notes)


def check_ci_erf(p: Path) -> Result:
    d = pd.read_csv(p, index_col=0)
    need("CH4" in d.columns, "no CH4 column")
    need(int(d.index[0]) == 1750 and (np.diff(d.index.astype(int)) == 1).all(), "years not 1750.. contiguous")
    need(d.CH4.between(-0.01, 1.0).all(), "CH4 ERF outside 0..1 W/m2")
    return Result(f"{len(d)} years {int(d.index[0])}..{int(d.index[-1])}; CH4 ERF {d.CH4.iloc[-1]:.3f} W/m2 in {int(d.index[-1])}",
                  ["CH4 ERF here includes only the direct CH4 term. CH4's indirect effects (tropospheric O3, "
                   "stratospheric H2O) are in the O3 and H2O_stratospheric columns of ERF_best_aggregates.csv."])


def check_ci_agg(p: Path) -> Result:
    d = pd.read_csv(p, index_col=0)
    for c in ("CH4", "O3", "H2O_stratospheric", "total", "anthro"):
        need(c in d.columns, f"column {c} missing")
    need(d.total.between(-7, 5).all(), "total ERF out of range")
    need(d.total.loc[1950:].between(-3, 5).all(), "total ERF after 1950 out of range")
    worst = d.total.idxmin()
    return Result(f"{len(d)} years {d.index[0]}..{d.index[-1]}; total ERF {d.total.iloc[-1]:.2f} W/m2 at the end",
                  ["Time index is mid-year (e.g. 2024.5).",
                   f"Most negative total ERF {d.total.min():.1f} W/m2 in {worst} (Tambora eruption, volcanic forcing). "
                   "Large negative spikes in volcanic years are real, not errors."])


def check_ci_slcf(p: Path) -> Result:
    d = pd.read_csv(p, index_col=0)
    for c in ("NOx", "CO", "NMVOC", "CH4"):
        need(c in d.columns, f"column {c} missing")
    need(d.CO.between(100, 1500).all() and d.NOx.between(5, 250).all(), "CO or NOx out of range")
    ch4 = d.CH4.dropna()
    need(ch4.between(50, 500).all(), "CH4 emissions out of range")
    return Result(f"{len(d)} years {d.index[0]}..{d.index[-1]}; CH4 emissions only from {ch4.index[0]}",
                  [f"Anthropogenic + biomass-burning emissions (Tg/yr; NOx as Tg NO2). In {d.index[-1]}: "
                   f"CH4 {d.CH4.iloc[-1]:.0f}, NOx {d.NOx.iloc[-1]:.0f}, CO {d.CO.iloc[-1]:.0f} Tg/yr.",
                   f"CH4 column is empty before {ch4.index[0]}."])


def check_ci_emis(p: Path) -> Result:
    d = pd.read_csv(p, index_col=0)
    need("CH3CCl3" in d.columns, "no CH3CCl3 column")
    m = d["CH3CCl3"]
    need(m.between(0, 1000).all(), "MCF emissions out of range")
    peak = m.idxmax()
    return Result(f"{len(d)} years {d.index[0]}..{d.index[-1]}; CH3CCl3 peak {m.max():.0f} Gg/yr in {peak}",
                  ["Emissions here are derived by inverting the concentration record (not inventories); "
                   "units are Gg/yr of each gas. Index is mid-year."])


def check_hadcrut(p: Path) -> Result:
    d = pd.read_csv(p)
    need(list(d.columns[:2]) == ["Time", "Anomaly (deg C)"], "unexpected columns")
    need(d.Time.iloc[0] == 1850 and (np.diff(d.Time) == 1).all(), "years not 1850.. contiguous")
    need(d["Anomaly (deg C)"].between(-1, 2).all(), "anomaly out of range")
    notes = ["Anomalies are relative to 1961-1990."]
    if d.Time.iloc[-1] >= TODAY.year:
        notes.append(f"The last row ({d.Time.iloc[-1]}) is a PARTIAL year (months to date). Drop it for annual work.")
    return Result(f"{len(d)} years {d.Time.iloc[0]}..{d.Time.iloc[-1]}", notes)


def check_rcmip(p: Path) -> Result:
    d = pd.read_csv(p)
    need(set(d.Scenario) >= set(RCMIP_SCENARIOS), "missing scenarios")
    var = "Emissions|CH4" if "emissions" in p.name else "Atmospheric Concentrations|CH4"
    need(var in set(d.Variable), f"no {var}")
    r = d[(d.Variable == var) & (d.Scenario == "historical")].iloc[0]
    v = float(r["2014"])
    lo, hi = (300, 450) if "emissions" in p.name else (1750, 1900)
    need(lo < v < hi, f"{var} 2014 = {v}")
    years = [c for c in d.columns if c.isdigit()]
    return Result(f"{len(d)} rows ({d.Variable.nunique()} variables x {d.Scenario.nunique()} scenarios), "
                  f"years {years[0]}..{years[-1]}; historical {var} 2014 = {v:.1f} {r.Unit}",
                  ["Historical rows end in 2014; SSP rows start in 2015. Many SSP years are only every 5-10 years "
                   "(interpolate)." if "emissions" in p.name else
                   "Concentrations are MAGICC7 projections used to drive CMIP6 models."])


# --------------------------------------------------------------------------------------------
# Special download steps
# --------------------------------------------------------------------------------------------
def get_sites(force: bool) -> None:
    folder, base, sites, prefix, suffix = RAW / "noaa_gml/sites_ch4", SITE_TXT, CH4_SITES, "ch4", "surface-flask_1_ccgg"
    for code in sites:
        dest = folder / f"{prefix}_{code}_{suffix}_month.txt"
        if dest.exists() and not force:
            continue
        print(f"get   {dest.relative_to(DATA)}")
        download(f"{base}/{dest.name}", dest)
    table = folder / "site_table.csv"
    if not table.exists() or force:
        print("get   site coordinates from NOAA event files (first 40 kB of each)")
        rows = [{"site": c, **site_meta(SITE_TXT, "ch4", c, "surface-flask_1_ccgg")} for c in CH4_SITES]
        pd.DataFrame(rows).sort_values("latitude", ascending=False).to_csv(table, index=False)


def get_edgar(force: bool) -> None:
    folder = RAW / "edgar"
    z = folder / "EDGAR_CH4_1970_2024.zip"
    if not z.exists() or force:
        print(f"get   {z.relative_to(DATA)}")
        download(f"{EDGAR}/{z.name}", z)
    with zipfile.ZipFile(z) as zf:
        bad = zf.testzip()
        need(bad is None, f"corrupt member {bad} in {z.name}")
        for m in zf.namelist():
            target = (folder / m).resolve()
            need(folder.resolve() in target.parents, f"unsafe path {m} in zip")
            if not target.exists() or force:
                zf.extract(m, folder)


def get_rcmip(force: bool) -> None:
    folder = RAW / "rcmip"
    jobs = [("rcmip-emissions-annual-means-v5-1-0.csv", RCMIP_EMIS_VARS, "rcmip_ch4_nox_co_voc_mcf_emissions_world.csv"),
            ("rcmip-concentrations-annual-means-v5-1-0.csv", RCMIP_CONC_VARS, "rcmip_ch4_mcf_concentrations_world.csv")]
    for src, variables, out in jobs:
        dest = folder / out
        if dest.exists() and not force:
            continue
        full = folder / ("_full_" + src)
        print(f"get   {src} ({'48' if 'emissions' in src else '21'} MB, filtered to World + CH4/NOx/CO/VOC/MCF)")
        download(f"{RCMIP}/{src}", full)
        d = pd.read_csv(full)
        keep = d[(d.Region == "World") & d.Variable.isin(variables) & d.Scenario.isin(RCMIP_SCENARIOS)]
        keep = keep.dropna(axis=1, how="all")
        keep.to_csv(dest, index=False)
        full.unlink()


# --------------------------------------------------------------------------------------------
# Source registry
# --------------------------------------------------------------------------------------------
@dataclass
class Source:
    key: str
    path: str                      # file (or folder) under data/raw/ that the check runs on
    url: str
    check: Callable[[Path], Result]
    fetch: Callable[[bool], None] | None = None   # custom download; default = plain download of url
    cookie: str | None = None


SOURCES = [
    Source("noaa_global", "noaa_gml/ch4_mm_gl.txt", f"{TRENDS}/ch4_mm_gl.txt", check_ch4_mm),
    Source("noaa_global", "noaa_gml/ch4_annmean_gl.txt", f"{TRENDS}/ch4_annmean_gl.txt", check_ch4_ann),
    Source("noaa_global", "noaa_gml/ch4_gr_gl.txt", f"{TRENDS}/ch4_gr_gl.txt", check_ch4_gr),
    Source("noaa_sites", "noaa_gml/sites_ch4", f"{SITE_TXT}/", check_ch4_sites, get_sites),
    Source("hats_mcf", "noaa_hats/HATS_global_MC.txt",
           f"{GML}/aftp/data/hats/solvents/CH3CCl3/combined/HATS_global_MC.txt", check_hats_mcf),
    Source("edgar", "edgar", f"{EDGAR}/EDGAR_CH4_1970_2024.zip", check_edgar, get_edgar),
    Source("owid", "owid/owid-co2-codebook.csv",
           "https://raw.githubusercontent.com/owid/co2-data/master/owid-co2-codebook.csv", lambda p: Result("codebook")),
    Source("owid", "owid/owid-co2-data.csv",
           "https://raw.githubusercontent.com/owid/co2-data/master/owid-co2-data.csv", check_owid),
    Source("gcp_budget", "gcp/Global_methane_budget_2023_2000_2020_v1.xlsx",
           f"https://data.icos-cp.eu/objects/{GCP_ID}", check_gcp, cookie=f"CpLicenseAcceptedFor={GCP_ID}"),
    Source("climate_indicator", "climate_indicator/ghg_concentrations.csv", f"{CI}/output/ghg_concentrations.csv", check_ci_conc),
    Source("climate_indicator", "climate_indicator/ERF_best.csv", f"{CI}/output/ERF_best.csv", check_ci_erf),
    Source("climate_indicator", "climate_indicator/ERF_best_aggregates.csv", f"{CI}/output/ERF_best_aggregates.csv", check_ci_agg),
    Source("climate_indicator", "climate_indicator/slcf_emissions.csv", f"{CI}/output/slcf_emissions.csv", check_ci_slcf),
    Source("climate_indicator", "climate_indicator/ghg_equivalent_emissions.csv",
           f"{CI}/output/ghg_equivalent_emissions.csv", check_ci_emis),
    Source("hadcrut", "hadcrut/HadCRUT.5.2.0.0.analysis.summary_series.global.annual.csv", HADCRUT, check_hadcrut),
    Source("rcmip", "rcmip/rcmip_ch4_nox_co_voc_mcf_emissions_world.csv", f"{RCMIP}/rcmip-emissions-annual-means-v5-1-0.csv",
           check_rcmip, get_rcmip),
    Source("rcmip", "rcmip/rcmip_ch4_mcf_concentrations_world.csv", f"{RCMIP}/rcmip-concentrations-annual-means-v5-1-0.csv",
           check_rcmip, lambda f: None),
]


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except AttributeError:
        pass
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true", help="re-download files that already exist")
    ap.add_argument("--only", nargs="+", metavar="KEY", help="only these source keys")
    ap.add_argument("--list", action="store_true", help="list source keys and exit")
    a = ap.parse_args()
    keys = list(dict.fromkeys(s.key for s in SOURCES))
    if a.list:
        print("\n".join(keys))
        return 0
    todo = [s for s in SOURCES if not a.only or s.key in a.only]
    if a.only and (unknown := set(a.only) - set(keys)):
        ap.error(f"unknown keys {sorted(unknown)}; choose from {keys}")

    RAW.mkdir(parents=True, exist_ok=True)
    failures, report = [], []
    for s in todo:
        dest = RAW / s.path
        try:
            if s.fetch is not None:
                s.fetch(a.force)
            elif a.force or not dest.exists():
                print(f"get   {dest.relative_to(DATA)}")
                download(s.url, dest, cookie=s.cookie)
        except Exception as e:  # noqa: BLE001 - report and continue with the other sources
            failures.append(f"{s.path}: download failed: {e}")
            print(f"FAIL  {s.path}: download failed: {e}")
            continue
        try:
            res = s.check(dest)
            print(f"ok    {s.path}: {res.summary}")
            report.append((s, res))
        except Exception as e:  # noqa: BLE001
            failures.append(f"{s.path}: check failed: {e}")
            print(f"FAIL  {s.path}: check failed: {e}")

    # Manifest of every file under raw/, with the check summary of the source it belongs to.
    summaries = {s.path: r.summary for s, r in report}
    url_of = {s.path: s.url for s in SOURCES}
    rows = []
    for f in sorted(p for p in RAW.rglob("*") if p.is_file() and not p.name.endswith(".part")):
        rel = f.relative_to(RAW).as_posix()
        owner = rel if rel in url_of else next((k for k in url_of if rel.startswith(k + "/")), "")
        url = url_of.get(owner, "")
        if url.endswith("/"):
            url += f.name if f.suffix == ".txt" else "(derived: coordinates read from NOAA event files)"
        rows.append({"file": f"raw/{rel}", "source_url": url, "bytes": f.stat().st_size, "sha256": sha256(f),
                     "file_date": dt.datetime.fromtimestamp(f.stat().st_mtime).date().isoformat(),
                     "check": summaries.get(rel, summaries.get(owner, ""))})
    with open(DATA / "MANIFEST.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    with open(DATA / "CHECK_REPORT.txt", "w", encoding="utf-8") as fh:
        fh.write(f"Check report written {dt.datetime.now():%Y-%m-%d %H:%M} by scripts/download_data.py\n\n")
        for s, r in report:
            fh.write(f"[ok] {s.path}\n     {r.summary}\n")
            for n in r.notes:
                fh.write(f"     note: {n}\n")
            fh.write("\n")
        for f in failures:
            fh.write(f"[FAIL] {f}\n")
    print(f"\n{len(rows)} files in data/raw; manifest -> data/MANIFEST.csv; notes -> data/CHECK_REPORT.txt")
    if failures:
        print(f"{len(failures)} problem(s):\n  " + "\n  ".join(failures))
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

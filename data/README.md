# Data for the methane–OH box-model project

Everything here was downloaded on **2026-11-09** by `../scripts/download_data.py`. Re-run it from the project folder to rebuild or update this folder:

```
uv run scripts/download_data.py           # fetch what is missing, then check everything
uv run scripts/download_data.py --force   # re-download everything (e.g. after NOAA's monthly update)
```

- **Checks.** After each download the script parses the file and tests its time coverage and value ranges. It also compares overlapping datasets with each other, e.g. NOAA annual means against the ClimateIndicator concentrations.
- **Outputs.** `MANIFEST.csv` lists every file with its source URL, size and SHA-256. `CHECK_REPORT.txt` has the full check output and notes. All checks passed on 2026-11-09 and again on 2026-11-11, after the isotope files were removed.
- **Access.** No manual steps are needed: every source is public and none needs a login. The one special case is the GCP budget; see [Manual fallback](#manual-fallback-if-a-scripted-download-ever-breaks).
- **Size.** 42 files, about 27 MB in total.
- **Cache.** `processed/edgar_*.csv` are CSV copies of the EDGAR sheets, written by `analysis/common.py` because reading the xlsx takes about 40 s. They are rebuilt automatically if the xlsx is newer, and can be deleted at any time.

## What each file is

### `raw/noaa_gml/`: the observed methane record (NOAA Global Monitoring Laboratory)

| File | What it is | Coverage | Units |
|---|---|---|---|
| `ch4_mm_gl.txt` | Global monthly mean CH₄ over the marine surface: `average`, deseasonalised `trend`, and their uncertainties | 1983-07 – 2026-06 (516 months) | ppb |
| `ch4_annmean_gl.txt` | Global annual mean CH₄ | 1984 – 2025 (1644.8 → 1936.0 ppb) | ppb |
| `ch4_gr_gl.txt` | Annual increase, from 1 Jan to 31 Dec | 1984 – 2025 | ppb/yr |
| `sites_ch4/ch4_<site>_surface-flask_1_ccgg_month.txt` | Monthly CH₄ at 23 clean-air flask sites from Alert (82°N) to the South Pole | 1983 – 2025-12 | ppb |
| `sites_ch4/site_table.csv` | Latitude, longitude and altitude of those sites, read by the script from NOAA's event files | — | °, m |

- **Source:** <https://gml.noaa.gov/ccgg/trends_ch4/> (global files) and <https://gml.noaa.gov/aftp/data/trace_gases/> (site files).
- **Citation:** Lan, X., Thoning, K.W. & Dlugokencky, E.J., *Trends in globally-averaged CH₄, N₂O and SF₆*, NOAA GML, doi:10.15138/P8XG-AA10.
- **Used for:** fitting every model; the stall, restart and surge; the NH–SH two-box model.

### `raw/noaa_hats/HATS_global_MC.txt`: methyl chloroform, a tracer for north–south mixing

- **What it is:** monthly hemispheric and global means of CH₃CCl₃ (MCF), plus 12 individual sites.
  - Global series: 1992-05 – 2025-03 (125.6 → 0.82 ppt).
  - NH-only rows go back to 1990-05.
- **Why it matters:** MCF was emitted almost only in the north and its emissions have now nearly stopped. Its NH–SH difference therefore measures how fast air mixes between the hemispheres (step 5).
- **Note:** using MCF's decay rate to infer OH *trends* was tried and dropped. It is outside the lectures, and with the emission data available the result was circular (see `PLAN.md`).
- **Source:** <https://gml.noaa.gov/aftp/data/hats/solvents/CH3CCl3/combined/>.
- **Citation:** Dutton, G.S., Hall, B.D., Montzka, S.A. et al. (2026), version 2026-06-11.

### `raw/edgar/`: anthropogenic CH₄ emissions by country and sector (EU JRC EDGAR 2025)

- **Files:** `EDGAR_CH4_1970_2024.zip` is the original download. The script unpacks `EDGAR_CH4_1970_2024.xlsx` and EDGAR's own `_readme.html` from it.
- **Coverage:** 1970 – 2024, 223 countries plus international aviation and shipping, in **Gg CH₄/yr**.
- **Sheets:**
  - `TOTALS BY COUNTRY`: one row per country.
  - `IPCC 2006`: country × IPCC sector, with a `fossil_bio` column.
  - `IPCC 1996`: same data, older sector codes.
  - The table header is on row 10, so read with `header=9`.
- **Global total:** 205 Tg (1970) → 317 Tg (2020) → **333 Tg (2024)**.
- **Largest sectors in 2024:**

| Sector | Tg/yr |
|---|---|
| Enteric fermentation | 114 |
| Oil and gas | 54 |
| Coal (solid fuels) | 53 |
| Landfills | 32 |
| Rice | 25.5 |
| Wastewater | 25 |
| Manure | 12 |

- **India:** 16.1 → 30.5 Tg (9.2 % of the world total).
- **Source:** <https://edgar.jrc.ec.europa.eu/dataset_ghg2025>.
- **Citation:** Crippa, M. et al., *GHG emissions of all world countries – 2025 report*, JRC.
- **Used for:** the anthropogenic emission input E(t), the Global Methane Pledge baseline (2020), and India's share.

### `raw/owid/`: Our World in Data CO₂ and GHG dataset

- **Files:** `owid-co2-data.csv` (50,411 rows, 254 countries and regions, 1750 – 2024) and `owid-co2-codebook.csv` (definitions of every column).
- **What's in it:** the `methane` column (1850 – 2024) comes from Jones et al. (2025), *National contributions to climate change*, and includes land use. It also has `population` and `gdp`, which the L18 economy links can use (methane per person, methane per unit GDP).
- **Source:** <https://github.com/owid/co2-data>.

### `raw/gcp/Global_methane_budget_2023_2000_2020_v1.xlsx`: Global Carbon Project methane budget

- **What it is:** the data supplement of Saunois et al. (2025), *Global Methane Budget 2000–2020*, ESSD 17, 1873, doi:10.5194/essd-17-1873-2025. Dataset doi: 10.18160/GKQ9-2RHT (ICOS Carbon Portal, CC BY 4.0).
- **Sheets:**
  - **Top-down (TD) and bottom-up (BU) budgets:** for 2000–2009, 2010–2019 and 2020. Values are in Tg CH₄/yr, given per source or sink category, per inversion model, for 19 regions, the globe and latitude bands.
  - **Observed global CH₄:** from four networks (NOAA, AGAGE, CSIRO, UCI).
  - **Anthropogenic emission history:** with SSP projections.
- **2010–2019 top-down means:**

| Category | Tg CH₄/yr (model range) |
|---|---|
| Total sources | 575 (553–586) |
| Anthropogenic | 369 |
| Natural | 206 |
| Wetlands | 165 |
| Total sinks | 550 |
| Chemical loss (mainly OH) | 521 |
| Soil uptake | 35 |

- **Imbalance:** sources minus sinks ≈ 25 Tg/yr, which is about the observed ~9 ppb/yr growth.
- **Used for:** the natural-emission term, the OH sink and lifetime, and checking the fitted box model against an independent budget.

### `raw/climate_indicator/`: concentrations, forcing and short-lived emissions since 1750

Source: the *Indicators of Global Climate Change* project (Forster et al.), <https://github.com/ClimateIndicator/forcing-timeseries>. The script pins release **v6.4.2 (2026-06-15)** so that it is repeatable.

| File | What it is | Coverage |
|---|---|---|
| `ghg_concentrations.csv` | Global mean concentrations of 55 gases, including CH₄ (ppb) and CH₃CCl₃ (ppt). Ice cores before 1979, NOAA/AGAGE after. | 1750, then 1850 – 2025 |
| `ERF_best.csv` | Effective radiative forcing per agent (W/m²); CH₄ = 0.570 W/m² in 2025 | 1750 – 2025 |
| `ERF_best_aggregates.csv` | Forcing grouped: CO₂, CH₄, O₃, stratospheric H₂O, aerosols, volcanic, `anthro`, `total` | 1750.5 – 2025.5 |
| `slcf_emissions.csv` | Global anthropogenic + biomass-burning emissions of NOx, CO, NMVOC, SO₂, BC, OC, NH₃, CH₄, N₂O (Tg/yr) | 1750 – 2025 |
| `ghg_equivalent_emissions.csv` | Emissions of halocarbons back-calculated from their concentrations (Gg/yr), including **CH₃CCl₃** (peak 692 Gg/yr in 1989) | 1750.5 – 2024.5 |

What these files are used for:
- **Spin-up:** running the model from pre-industrial (729 ppb in 1750), i.e. the long history before NOAA's record starts.
- **Forcing check:** testing our square-root forcing formula against the official CH₄ ERF.
- **Warming:** total forcing for the L03/L04 energy-balance model.
- **OH drivers:** NOx raises OH and CO lowers it.
- **MCF tracer:** the emissions needed for the two-hemisphere tracer fit.

### `raw/hadcrut/HadCRUT.5.2.0.0.analysis.summary_series.global.annual.csv`: global temperature

- **What it is:** annual global-mean surface temperature anomaly (°C, relative to 1961–1990) with a 95 % range, 1850 – 2026.
- **Source:** Met Office Hadley Centre, <https://www.metoffice.gov.uk/hadobs/hadcrut5/>. Citation: Morice et al. (2021), JGR, doi:10.1029/2019JD032361.
- **Used for:**
  - calibrating the L04 energy-balance model (C dT/dt = F − λT) and turning methane forcing into warming;
  - the temperature dependence of the OH rate constant and of wetland emissions.

### `raw/rcmip/`: future scenarios (RCMIP v5.1.0, the CMIP6 SSP inputs)

Both files keep only the World rows of 10 scenarios: historical, SSP1-1.9, SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP3-7.0-lowNTCF, SSP4-3.4, SSP4-6.0, SSP5-3.4-over and SSP5-8.5. The original files (48 MB and 21 MB) are filtered on download and not kept.

| File | Variables | Units |
|---|---|---|
| `rcmip_ch4_nox_co_voc_mcf_emissions_world.csv` | CH₄ (total, AFOLU, fossil and industrial), NOx, CO, VOC and CH₃CCl₃ emissions | Mt CH₄/yr, Mt NOx/yr, Mt CO/yr, Mt VOC/yr, kt/yr |
| `rcmip_ch4_mcf_concentrations_world.csv` | CH₄ and CH₃CCl₃ concentrations from MAGICC7 | ppb, ppt |

- **Source:** <https://www.rcmip.org/>.
- **Citation:** Nicholls et al. (2020), GMD, doi:10.5194/gmd-13-5175-2020. The SSP emissions are from Gidden et al. (2019).
- **Used for:** the business-as-usual and mitigation scenarios that the Global Methane Pledge is compared with, and a check of our model's projections against MAGICC.

## Things to know before using the data (odd things found)

**NOAA global record**
1. **The newest NOAA months are preliminary.** In `ch4_mm_gl.txt`, 2026-02 to 2026-06 have uncertainty **−9.99**, NOAA's missing-value code. Treat it as NaN. These months will change at the next update.
2. **"Growth rate" depends on its definition.** `ch4_gr_gl.txt` is the 1 Jan → 31 Dec increase. The year-to-year difference of annual means can differ from it by up to **8 ppb/yr**:

   | Year | 1 Jan → 31 Dec | Difference of annual means |
   |---|---|---|
   | 2004 | −4.7 | −0.2 |
   | 2014 | 12.7 | 9.2 |
   | 2025 | 4.8 | 6.6 |

   Pick one definition and state it. `trend` in the monthly file gives a smooth derivative.
3. **The record in numbers.**
   - The 2000–2006 **stall** averaged 0.5 ppb/yr.
   - The fastest growth was **2021 (17.7)**, then 2020 (14.8) and 1991 (14.0) ppb/yr.
   - Growth has since slowed: 2025 was 4.8 ppb/yr, the lowest since 2009.

**NOAA site files**

4. **Some sites stopped or have gaps.**
   - The CH₄ records at azr (ends 2022-03), chr (2019-12), eic (2019-11) and gmi (2023-03) stop early.
   - azr and gmi each have 11 multi-month gaps.
   - izo, mlo and spo are above 1500 m, so they sample the free troposphere, not the marine boundary layer.
   - The monthly site files carry no coordinates, which is why the script builds `site_table.csv`.
5. **Two event-file header styles.** NOAA uses two header formats in its event files: older ones have `# data_fields:` and newer ones an uncommented `site_code …` line. The script handles both.

**Methyl chloroform**

6. **MCF is now below 1 ppt**, close to the measurement noise.
   - Step 5 uses 1993–2015, when the NH–SH difference is still well above the noise.
   - A quick fit gives an e-folding time of **5.7 yr** for 2005–2020. This is the *total* lifetime: OH plus ocean and stratosphere losses, plus small leftover emissions. It is not the OH lifetime alone.

**Emission inventories**

7. **Emission inventories disagree by about 10–15 %.**

   | Year | EDGAR 2025 | CEDS + fires (ClimateIndicator) | RCMIP (CEDS-based) |
   |---|---|---|---|
   | 2014 | 307 | 348 | 388 |
   | 2024 | 333 | 368 | — |

   All values are Tg CH₄/yr. EDGAR leaves out large-scale wildfire and savanna burning, which are about 15–25 Tg/yr according to GCP. The GCP top-down anthropogenic total for 2010–2019 is 369. In a fitted model the natural-emission term absorbs this offset, so report which inventory was used.
8. **EDGAR quirks.**
   - Units are Gg, not Tg.
   - The country sheet includes "Int. Aviation" and "Int. Shipping" rows.
   - The zip also contains EDGAR's `_readme.html`.
9. **OWID `methane` is in million tonnes CO₂-equivalent (GWP100), not Tg CH₄.**
    - World 2024 = 9,499 Mt CO₂e.
    - There is no single conversion, because fossil and biogenic CH₄ have different GWPs (29.8 and 27.0 in AR6). Dividing by about 28 gives roughly 340 Tg.
    - Use EDGAR when you need physical units.

**ClimateIndicator and GCP**

10. **`ghg_concentrations.csv` has no rows for 1751–1849.** It jumps from 1750 (729.2 ppb) to 1850 (807.6 ppb), so interpolate if needed. After 1984 it agrees with NOAA to within ±3 ppb.
11. **Time stamps differ between files.** `ERF_best_aggregates.csv` and `ghg_equivalent_emissions.csv` use mid-year time stamps (2024.5); the others use integer years. In `slcf_emissions.csv` the **CH₄ column is empty before 1970**.
12. **Methane's forcing is split across columns.** `CH4` in the ERF files is only the direct effect. Methane's indirect effects (tropospheric ozone, stratospheric water vapour) are in the `O3` and `H2O_stratospheric` columns.
13. **Volcanic spikes in forcing are real.** Total ERF reaches −5.7 W/m² in 1816 (Tambora).
14. **GCP workbook layout.**
    - Every sheet starts with a block of free-text references. Find the header row whose first cell is `category` (row 17 in the top-down sheets), and remember that header rows differ between sheets.
    - In the observation sheet, the **AGAGE series disagrees with NOAA before about 1990**. It is 82 ppb low in 1985 and 11 ppb low in 1988, an artefact from when the network was starting up. Use AGAGE only after 1990.
    - Wetland estimates vary a lot between inversions (145–214 Tg/yr).

**HadCRUT**

15. **The last HadCRUT row (2026) is a partial year.** Drop it for annual analysis.

## Not downloaded, and why

- **Global OH itself.** It has no direct measurement; the literature OH lifetime and the GCP chemical sink stand in for it.
- **δ¹³C-CH₄ isotope data.** These were downloaded at first (NOAA/INSTAAR, 15 sites) but removed on 2026-11-11 together with the isotope analysis (see `PLAN.md`).
- **Gridded EDGAR maps and the full GCP zip (167 MB of gridded fluxes).** A box model only needs global and national totals.
- **AGAGE direct downloads.** `agage2.eas.gatech.edu` did not respond on 2026-11-09. The AGAGE global CH₄ series is inside the GCP workbook, and the ClimateIndicator concentrations already blend AGAGE in.
- **Satellite CH₄ (TROPOMI, GOSAT).** These would be needed for regional plume detection, which is outside the scope of a box model.

## Manual fallback (if a scripted download ever breaks)

1. **GCP budget (ICOS licence page).** The script accepts the ICOS licence with the cookie the portal sets. If ICOS changes this:
   1. Open <https://www.icos-cp.eu/GCP-CH4-2024> and click **Download** next to "Global CH4 Budget 2000-2020: Global emission data".
   2. Accept the CC BY 4.0 licence.
   3. Save `Global_methane_budget_2023_2000_2020_v1.xlsx` into `data/raw/gcp/`.
   4. Re-run the script; it will skip the download and only check the file.
2. **EDGAR.** If the JRC FTP path moves, open <https://edgar.jrc.ec.europa.eu/dataset_ghg2025>, download "CH4 – emissions by country and sector (1970–2024)" (`EDGAR_CH4_1970_2024.zip`), put it in `data/raw/edgar/` and re-run.
3. **Any other file.** The URL for every file is in `MANIFEST.csv`. Download the file to the path listed there, then re-run the script to check it.

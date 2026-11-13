# Project: the methane–OH box model

**Course:** EH-605 Modeling Earth System and Sustainability, IIT Gandhinagar (Prof. R N Singh)
**Student:** Shardul Junagade (23110297)
**Updated:** 2026-11-09. Replaces the 2026-11-07 version. The data are now downloaded (see `data/README.md`).

> **Status:**
> - **2026-11-09:** implemented.
> - **2026-11-11:** checked against the lecture content and trimmed. Two items from the idea list below were dropped: the methyl-chloroform OH inference and the δ¹³C isotope budget. Both are outside the lectures and gave weak results; `PLAN.md` explains why.
>
> See `PLAN.md` for the final plan, the course-coverage check and the decisions taken, and `RESULTS.md` for every result and conclusion. Run everything with `uv run scripts/run_all.py`.

## Title

**"Why Is Methane Rising Again? A Methane–OH Box Model of Recent Trends and the Global Methane Pledge"**

## The question

Atmospheric methane rose from 729 ppb in 1750 to 1936 ppb in 2025. The NOAA record (1983–2026) has four features that any model has to explain:

| Feature | When | Observed |
|---|---|---|
| Stall | 2000–2006 | growth ≈ 0.5 ppb/yr |
| Restart | 2007 | sudden jump back to growth |
| Record surge | 2020–2021 | 14.8 and 17.7 ppb/yr, the fastest since records began |
| Slowdown | 2025 | 4.8 ppb/yr |

There are two competing explanations: changes in **emissions** (wetlands, agriculture, fossil fuels) or changes in the **sink**, the hydroxyl radical OH. Methane also uses up OH. That makes the sink nonlinear and gives a positive feedback: more methane means less OH, so methane lives longer.

The project builds a hierarchy of box models, from simple to less simple, confronts each with data, and then uses them to ask what the Global Methane Pledge (−30 % anthropogenic CH₄ by 2030 from 2020 levels) would do to concentration, forcing and temperature.

## How the project uses the course

| Course idea (lecture) | Where it appears in the project |
|---|---|
| Modelling loop: object → variables → equations → confront with observations → predict (L01, L02) | Every model is fitted to NOAA data and judged on what it gets wrong |
| Greenhouse effect, radiative forcing, energy-balance model (L03) | CH₄ forcing ΔF(M) is checked against official ERF; forcing is converted to warming |
| Equilibrium, linear stability, feedback parameter, sensitivity, tipping and hysteresis (L04) | Equilibria of the CH₄–OH system, Jacobian eigenvalues, feedback factor, an OH-collapse threshold |
| Box/compartment models, residence time, analytic solutions, two- and three-box models (L11) | One-box CH₄ model; NH–SH two-box model with exchange |
| Mass-action chemical kinetics, rate constants, Arrhenius law (L15) | CH₄ + OH → products, at rate k(T)[OH][CH₄], with k = 2.45×10⁻¹² e^(−1775/T) |
| Holling type II and Michaelis–Menten saturation (L07, L13) | With OH in quasi-steady state, the CH₄ sink becomes a saturating function of CH₄: the same mathematics as nutrient uptake |
| Feedback loops between climate and the biosphere (L04, L10 Daisyworld) | Warmer climate → more wetland CH₄ → more warming: the loop gain of a climate–methane feedback |
| Economy and the environment, intensity per person and per GDP (L18) | India's methane intensity and share; emissions per capita |
| Sustainability as dynamic optimisation, GEKKO/SciPy optimal control (L19) | The cheapest abatement path that meets the pledge or a warming target |

## What we'll do

The core (1–5) is enough for a complete project. Items 6–9 are extensions; do them in order if time allows.

1. **One-box model** (L11).
   - Model: dC/dt = E(t)/β − C/τ. C is CH₄ in ppb; β ≈ 2.75–2.84 Tg per ppb, the CH₄ analogue of L11's 2.123 GtC per ppm.
   - Use EDGAR anthropogenic E(t) plus a constant natural source, and fit τ and the natural source to NOAA 1984–2025. Solve it analytically, as in L11.
   - Show where it fails: the 2000–2006 stall, the 2007 restart and the 2020 surge.
   - Spin up from 1750 using the ClimateIndicator history.

2. **CH₄–OH model from mass-action kinetics** (L15, L04).
   - Methane: dC/dt = E − k(T)·X·C
   - OH: dX/dt = S − k·X·C − X/τ_X
   - X is OH, S is the OH source (driven by NOx) and τ_X is OH's lifetime from its other sinks (mostly CO).
   - OH lives about 1 second, so it is always in quasi-steady state: X = S/(kC + 1/τ_X). Then dC/dt = E − kSC/(kC + 1/τ_X), a Holling type II sink that saturates at S.
   - Equilibria and stability: compute F′(C*). Show that the perturbation lifetime exceeds the global lifetime, with feedback factor f = τ_pert/τ ≈ 1.3–1.4 in the literature, and derive it analytically.
   - If emissions exceed the OH supply (E > S), no equilibrium exists: OH runs out and methane runs away. That is a tipping point.
   - Draw a bifurcation diagram of C* against E and phase portraits. Report how far today's world is from the threshold.

3. **Fit and compare against observations** (L02).
   - Fit with `scipy.optimize.least_squares` and report parameter uncertainties (bootstrap).
   - Compare models by AIC and by leaving out the 2016–2025 period.
   - Test the "emissions versus OH" question:
     - drive S with NOx and τ_X with CO from `slcf_emissions.csv`;
     - check whether the 2020 NOx drop during COVID lowers OH enough to explain part of the 2020 surge.

4. ~~**Methyl chloroform as an OH thermometer**~~ *(dropped on 2026-11-11: outside the lectures, and the result was circular; see `PLAN.md`)*.

5. **Forcing, warming and the Global Methane Pledge** (L03, L04, L19).
   - Forcing: ΔF from the simplified formula, checked against `ERF_best.csv`.
   - Warming: an energy-balance model C dT/dt = F − λT, calibrated with total ERF and HadCRUT5.
   - Scenarios:
     - the RCMIP SSPs (SSP1-1.9 to SSP5-8.5, plus SSP3-7.0-lowNTCF);
     - the pledge applied to EDGAR's 2020 baseline;
     - an immediate stop to emissions.
   - Report the e-folding response time, the CH₄ peak, warming avoided by 2050 and 2100, and a comparison with MAGICC's concentrations.

6. **Interhemispheric two-box model** (L11 two-box with exchange).
   - N and S hemispheres, each with its own loss rate and an exchange time τ_ex.
   - Fit to the 23 NOAA sites (`site_table.csv` gives the latitudes). Get τ_ex (about 1 yr) and the share of emissions in each hemisphere, and compare with EDGAR's country totals.

7. ~~**δ¹³C isotope budget**~~ *(dropped on 2026-11-11: outside the lectures, complex, and the result was weak; see `PLAN.md`)*.

8. **Climate–methane feedback loop** (L04, L10).
   - Make wetland emissions depend on temperature (Q10 law) and the OH rate depend on temperature (Arrhenius), then couple them with the EBM.
   - Compute the loop gain and the extra warming this feedback adds in each scenario.

9. **Optimal abatement and India** (L18, L19).
   - Choose an abatement path u(t) that minimises a quadratic abatement cost plus a warming damage, subject to the box model, using GEKKO or SciPy, as in the L19 fishery problem.
   - India:
     - sectors (EDGAR: livestock and rice dominate);
     - methane per person and per unit GDP (OWID);
     - what a 30 % cut by India would mean for global CH₄.

**Output:** a script pipeline (`scripts/`), figures, and a results write-up in `RESULTS.md`, then a report of about 10–15 pages and a presentation.

## Data

All data are downloaded and checked; see `data/README.md` for every file, its source, citation and quirks.

| Data | Use |
|---|---|
| NOAA GML global CH₄: monthly, annual, growth rate (1983–2026) | Target for every fit |
| NOAA site CH₄ at 23 sites | NH–SH model |
| NOAA HATS methyl chloroform, NH/SH monthly (1992–2025) | Tracer for the NH–SH exchange time |
| EDGAR 2025 CH₄ by country and IPCC sector (1970–2024) | Anthropogenic emissions; pledge baseline; India |
| Our World in Data CO₂/GHG dataset (methane in CO₂e, population, GDP) | India and economy links |
| GCP Global Methane Budget 2000–2020 (Saunois et al. 2025) | Natural sources, OH sink, independent budget check |
| ClimateIndicator v6.4.2: concentrations since 1750, ERF, NOx/CO/CH₄ emissions, MCF emissions | Spin-up, forcing check, OH drivers, tracer emissions |
| HadCRUT 5.2 global temperature | EBM calibration |
| RCMIP v5.1.0 SSP emissions and concentrations | Future scenarios |

**Rebuild the data:** `uv run scripts/download_data.py` (about 2 minutes; no login needed).

## Topic check

According to the 2026 project sheet, no other student has a global methane-budget topic. Two topics are nearby but different:
- Ujwal's paddy-field GHG emissions under alternate wetting and drying (field-scale rice CH₄);
- Ginisha's permafrost carbon feedback.

Mention the difference in the introduction: this project is about the global atmospheric budget and OH chemistry.

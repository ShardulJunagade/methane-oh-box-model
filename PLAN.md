# Project plan: what we set out to do with the course, and where it stands

**Project:** *Why Is Methane Rising Again? A Methane–OH Box Model of Recent Trends and the Global Methane Pledge*
**Course:** EH-605 Modeling Earth System and Sustainability (IIT Gandhinagar)
**Student:** Shardul Junagade (23110297)
**Status:**
- 2026-11-09: implemented.
- 2026-11-11: checked against the lecture content and trimmed (see "Course coverage check" below).

All 9 steps run; the results are in `RESULTS.md`.

## Guiding idea

Use the course's toolbox on one real Earth-system question, following the course loop for every step: write the model → solve it → confront it with data → judge it. The toolbox:
- box models and residence time
- feedback, stability and tipping
- mass-action kinetics
- energy balance
- dynamic optimisation

Keep each model small enough to be solved and understood by hand. Add complexity only when the data show the simpler model fails, which is the L02 modelling loop.

## The plan, step by step

| # | Step | Course link | Method | Script |
|---|---|---|---|---|
| 1 | Look at the data: the record, growth rate, inventories, budget, gradients | L01–L02 (quantify the object), L11 (sources and sinks) | Descriptive statistics and plots | `01_explore_data.py` |
| 2 | One-box model dC/dt = E/β − C/τ | L11 (residence time, analytic solution, spin-up as in the Rodhe model), L02 (fitting) | Analytic solution, least squares, misfit surface, top-down mass balance, spin-up from 1750 | `02_one_box.py` |
| 3 | CH₄–OH system | L15 (mass-action kinetics), L04 (equilibria, feedback, tipping), L08 (Jacobian, phase plane), L13 (Michaelis–Menten / Holling) | Two-variable ODE; OH in quasi-steady state; eigenvalues; pulse test; equilibrium-vs-emission curve; phase portrait | `03_oh_feedback.py` |
| 4 | Confront the model hierarchy with data; attribute the stall, restart and surge | L02 (parameter estimation), L04 (sensitivity) | Least squares with a literature prior on the lifetime, AIC, out-of-sample test, two inventories, counterfactual runs | `04_fit_compare.py` |
| 5 | Two-hemisphere atmosphere | L11 (two-box model with exchange) | A tracer (methyl chloroform) fixes the exchange time; the CH₄ gradient then gives where emissions are | `05_two_hemisphere.py` |
| 6 | Forcing and warming | L03 (greenhouse effect, layer models), L04 (EBM, feedback parameter, sensitivity) | Three forcing formulas vs official ERF; one- and two-layer EBM fitted to HadCRUT5 | `06_forcing_warming.py` |
| 7 | Futures, the Global Methane Pledge, climate–methane loop | L11 (emission scenarios), L19 (scenarios), L04 (feedback loop gain) | SSP scenarios, the pledge, Monte Carlo; loop gain | `07_scenarios_pledge.py` |
| 8 | Optimal abatement | L19 (SLSQP, GEKKO optimal control, fishery problem), L18 (externality) | GEKKO IMODE 6: cost–benefit problem and target problem | `08_optimal_abatement.py` |
| 9 | India | L18 (economy and environment, per-capita quantities, Dasgupta "impact inequality") | Sectors, per-capita and per-GDP intensity, decomposition, the effect of a 30 % cut | `09_india.py` |

Shared code:
- `analysis/common.py`: data loaders, constants, plot style
- `analysis/models.py`: the CH₄ box model, its fitting, forcing formulas, EBM
- `analysis/scenarios.py`: SSP driver paths and the pledge path

Run everything with `uv run scripts/run_all.py --fresh` (about 2 minutes).

## Course coverage check (2026-11-11)

Each concept used was searched for in the lecture notebooks (`lectures/`), their text dumps (`lecture_notes/source_text/`) and the lecture notes (`lecture_notes/chapters/`).

### Taught in the lectures and used

| Concept | Where taught | Used in step |
|---|---|---|
| Modelling loop; confronting equations with observations; inverse problems | L02 | all, especially 2 and 4 |
| Box/compartment models, mass conservation, residence time | L11 (also L12, L13, L15) | 2, 5 |
| Analytic solution of a linear box model with exponential emissions | L11 (two-box solution, F = q e^{αt}) | 2 |
| Long spin-up with an emission history; emission scenarios (SRES A1–B2) | L11 (Rodhe & Björkström model) | 2, 7 |
| Two-box model with exchange between boxes | L11 (atmosphere–ocean) | 5 |
| Mass-action rate laws, rate constants, Arrhenius dependence k = A e^{−E/RT} | L15 (Arrhenius also in L17) | 3, 4 |
| Michaelis–Menten / Monod uptake; Holling type II saturation | L13 (Monod), L08/L11–L14 (Holling) | 3 |
| Equilibrium, linear stability (∂F/∂x < 0), feedback parameter, sensitivity s = −1/λ | L04 | 3, 6, 7 |
| Tipping points, thresholds, equilibrium vs parameter (hysteresis, washout) | L04 (snowball), L13 (Deng washout threshold) | 3 |
| Jacobian and eigenvalues of a 2-variable system | L08 (Lotka–Volterra), also L13, L15 | 3 |
| Phase plane / phase portrait | L14, L18, L19 | 3 |
| Greenhouse effect; forcing term a ln(CO₂/280) in the EBM | L03, L04 | 6 |
| Energy balance model C dT/dt = F − λT; feedback parameter B = 1.3 | L04 | 6, 7 |
| Two-layer (atmosphere) models | L03 | 6 (as motivation) |
| Numerical ODE solving (Euler, odeint, RK4) | L04/L09/L11 (Euler), most lectures (odeint), L19 (RK4) | all |
| Scenario thinking (baseline, business as usual, exploratory) | L19, L11 | 7 |
| Constrained optimisation (SLSQP), optimal control with GEKKO (IMODE 6), profit integral, discounting | L19 (fishery problem) | 8 |
| Externalities; per-capita output and intensity; Dasgupta impact inequality N·y/α | L18 | 8, 9 |

### Not taught, but kept (simple to explain, good results)

| Concept | Why it is kept | Step |
|---|---|---|
| Least squares with a literature prior on the lifetime | Without it, lifetime and natural emissions cannot be separated (step 2 shows the misfit "valley"). The prior is just one extra term in the sum of squares. | 4 |
| AIC and an out-of-sample (holdout) test | One-line formulas for comparing models fairly; they give clear answers (V4 best, holdout 6.7 ppb). | 4 |
| Quasi-steady state for OH (fast–slow system) | The step that turns mass action (L15) into the Michaelis–Menten form (L13). Error vs the full system is 4×10⁻⁶ ppb. | 3 |
| Two-layer EBM with a deep ocean | The L11 two-box idea applied to heat. γ is taken from the literature, and the fit gives λ = 1.25 (close to L04's B = 1.3). | 6 |
| Monte Carlo uncertainty | Rerunning the model with random parameter draws; a clear 5–95 % range. | 7 |
| Q10 temperature law for wetlands | Only as a what-if knob in the loop-gain calculation, which gives a clear answer (no runaway). | 7 |
| Area-weighted averaging of station data (interpolation in sin(latitude)) | Needed to turn 23 stations into hemispheric means; matches NOAA to 0.3 ppb. | 5 |
| Methane forcing formulas (Myhre 1998, Etminan 2016) | The methane analogue of L04's ln(CO₂) term; compared with the course's √M formula. | 6 |
| Social cost of methane, linear marginal-abatement-cost curve, log decomposition of emissions (population × income × intensity) | Simple economic bookkeeping on top of L18/L19. | 8, 9 |

### Removed (not taught, and the result was bad or mediocre)

| Removed part | Why it was removed |
|---|---|
| Methyl-chloroform inversion for OH trends (old step 4, and the MCF-OH model runs in the fit comparison) | New atmospheric chemistry. The only available MCF emission series was itself back-calculated with a fixed lifetime, so the inversion was partly circular. The OH trend came out anywhere from −14 to +7 %/decade, and prescribing it did not improve the methane fit. |
| δ¹³C isotope budget (old step 7, and the δ¹³C station data) | New and complex: δ notation, kinetic isotope fractionation, source signatures, a fixed-effects station index. The calibrated fractionation (8.8 ‰) disagreed with the literature (6–7 ‰), and the forward model did not reproduce the observed isotope decline. |
| Fitting the wetland temperature sensitivity Q10 (two model variants) | The fit gave Q10 < 1 (0.7 ± 0.3), which is unphysical, and fixing Q10 = 2 made every fit statistic worse. |
| Two-layer EBM with the ocean coupling γ fitted freely | The record could not determine it (λ ran to its bound, giving 13 K climate sensitivity). Replaced by the literature γ. |

### Kept although the result is mediocre (because the concept is from the lectures)

- **V1′ one-box with a free lifetime and the misfit surface (L02).** It shows *why* the one-box model fails (a 6.4-yr lifetime would be needed).
- **V2 feedback-only model (L04).** It fits worse than V1. That is itself the result: the feedback cannot explain the stall.
- **V3′ with the feedback strength fitted (L04).** It is not identifiable (0.07 ± 0.14), which justifies using the literature value of 0.31.
- **Spin-up from 1750 (L11).** Neither run matches the whole history, which shows that today's non-inventory emissions exceed the pre-industrial natural ones.

## Decisions taken along the way, and why

1. **Lifetime prior.** The one-box model cannot fit the record with a realistic lifetime, and with the lifetime free it runs to 6.4 yr. So the OH lifetime is constrained by the literature value: τ_OH = 11.2 ± 1.3 yr (Prather et al. 2012), giving a total lifetime of 9.1 yr.
2. **Feedback strength a.** It is fixed at the literature value of 0.31, because concentration data cannot determine it.
3. **Ocean coupling.** γ in the two-layer EBM is the CMIP mean (0.67 W/m²/K).
4. **Tracer emission split.** Methyl chloroform's NH emission share is fixed at 0.98 (an industrial solvent), so the hemispheric exchange time can be fitted.

## Possible extensions (not done)

- A three-box CH₄–CO–OH model with CO concentration as a state variable. That would make the CO effect explicit.
- A monthly model with the seasonal cycle.
- Using the 2023–2025 slowdown as a true out-of-sample test once more data arrive.

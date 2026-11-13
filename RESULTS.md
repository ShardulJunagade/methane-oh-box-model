# Results: why is methane rising again?

**Project:** *Why Is Methane Rising Again? A Methane–OH Box Model of Recent Trends and the Global Methane Pledge*
**Course:** EH-605 Modeling Earth System and Sustainability, IIT Gandhinagar
**Student:** Shardul Junagade (23110297)
**Run date:** 2026-11-11

Every number below is computed by the scripts and stored in `results/*.json`. All 20 figures are in `figures/`.
- `PLAN.md` has the step-by-step plan and the check of which concepts come from the lectures (and what was removed because it did not).
- `data/README.md` has the data sources and data quirks.

```
uv run scripts/download_data.py     # data (once)
uv run scripts/run_all.py --fresh   # all 9 steps, ~2 minutes
```

---

## The answer in brief

1. **A constant-lifetime box model cannot explain the record.**
   - With the known 9.1-yr lifetime it misses by up to 30 ppb (RMSE 21 ppb) and predicts no stall.
   - It reaches a fair fit only if methane lived 6.4 yr, which is inconsistent with the literature value.
   - This rules out "inventory emissions + constant sink" as the explanation.
2. **The CH₄–OH feedback is real but cannot explain the stall.**
   - Methane uses up OH, so a perturbation decays in **12.7 yr**, not 9.3 yr (feedback factor **1.36**).
   - That makes concentrations respond *more* to rising emissions, so the model's fit gets worse, not better (RMSE 32 ppb).
3. **What does explain the shape of the record** (RMSE 3 ppb, growth-rate correlation 0.80) is a sink that is *pushed from outside*:
   - OH rises when NOx emissions rise and falls when CO emissions rise;
   - the CH₄ + OH reaction speeds up with temperature (Arrhenius).
   In this model:
   - **Stall (2000–2006):** rising NOx (more OH) took **−1.6 ppb/yr** off growth.
   - **Restart (2007):** NOx began to fall.
   - **2020–21 surge:** the COVID drop in NOx added **+3.4 ppb/yr**. Part of this was offset by lower CO.
   - **Missing from the model:** it still under-predicts the surge by about 2 ppb/yr, which leaves room for extra (wetland) emissions.
4. **Where the extra methane comes from.**
   - A mass balance with the known lifetime gives total emissions of 574 Tg/yr for 2010–19 (the Global Carbon Project's independent estimate is 575).
   - Emissions not in the EDGAR inventory *fell* by about 20 Tg/yr during the stall and *rose* by about 25 Tg/yr in the surge.
   - The two-hemisphere model shows these non-inventory sources are only about 60 % in the north, so they are largely tropical.
5. **No methane tipping point is in reach.**
   - The OH-only system has a runaway threshold at **1425 Tg/yr, 2.5× today's emissions**.
   - Stratospheric, soil and Cl sinks remove the singularity altogether.
   - The climate–wetland–methane loop has a gain below 0.02 even for wetland Q10 = 4.
6. **The Global Methane Pledge** (−30 % of 2020 anthropogenic emissions by 2030):
   - lowers methane in 2050 by **≈ 330 ppb** relative to SSP2-4.5 (1674 vs 2005 ppb);
   - avoids **0.055 K (0.047–0.065 K)** of warming by 2050 through methane's direct forcing alone.
   - The benefit is fast, but it stops if the cuts stop.
7. **Optimisation (GEKKO, as in L19).**
   - With a social cost of methane of $1600/t, the optimal abatement is **≈ 35–36 %**, close to the pledge's 30 %.
   - The cheapest way to reach the pledge's 2050 concentration back-loads the cuts and costs 1.5× less. The pledge buys earlier benefits.
8. **India** emits 30.5 Tg/yr (9 % of the world), two-thirds from livestock and rice.
   - Per person it emits **20 kg/yr, half the world average**.
   - Its methane per unit GDP fell 78 % from 1990 to 2022.

---

## 1. What the data contain (step 1, `01_explore_data.py`)

| Quantity | Value |
|---|---|
| Global CH₄ | 729 ppb (1750) → 1645 ppb (1984) → **1936 ppb (2025)**: +166 % since pre-industrial |
| Mean growth | 8.5 ppb/yr (1985–99), **0.4 (2000–06)**, 5.5 (2007–13), 8.9 (2014–19), **14.8 (2020–22)**, 8.3 (2023–25) |
| Record years (Jan→Dec) | 2021: 17.7, 2020: 14.8, 1991: 14.0 ppb/yr |
| Anthropogenic emissions 2024 | EDGAR 333, CEDS + fires 368 Tg/yr (inventories differ by 10–15 %) |
| EDGAR sectors 2024 | livestock 126, waste 57, oil & gas 54, coal 53, rice 25.5 Tg/yr |
| EDGAR growth | **3.9 Tg/yr² during the 2000–06 stall**, 2.8 Tg/yr² for 2007–19: the inventory grows fastest exactly when methane stalled |
| GCP top-down budget 2010–19 | sources 575 (anthropogenic 369, natural 206, of which wetlands 165); sinks 550 (chemical 521, soil 35) Tg/yr |
| Arctic minus Antarctic CH₄ | ≈ 130–150 ppb |

Figures: `01_record_growth.png`, `01_inventories.png`, `01_sectors.png`, `01_mcf_gradient.png`. The last one shows the methyl-chloroform tracer used in step 5, and the latitude gradient.

---

## 2. The one-box model (step 2; L11 residence time, L02 fitting)

Model: dC/dt = E(t)/β − C/τ, with β = 2.75 Tg per ppb. This is the methane version of L11's 2.123 GtC per ppm. E is EDGAR anthropogenic emissions plus a constant natural source.

- **Analytic check (as in L11).** For exponential emissions E = q e^{αt}, the analytic solution
  C = C₀e^{−t/τ} + (E_n/β)τ(1 − e^{−t/τ}) + (q/β)/(α + 1/τ)·(e^{αt} − e^{−t/τ})
  agrees with the numerical solver to 0.08 ppb.
- **Two fits:**

  | Fit | Lifetime | Natural source | RMSE |
  |---|---|---|---|
  | Lifetime fixed (literature) | 9.15 yr | 274 Tg/yr | **21 ppb** |
  | Lifetime free | **6.4 yr** | 500 Tg/yr (at the bound) | 12.6 ppb |

  The fixed-lifetime model has growth-rate correlation 0.37 and no stall at all. Fitted on 1984–2005 and run forward, it overshoots 2025 by **39 ppb**.
- **Misfit surface** (`02_cost_surface.png`). The error has a long valley along τ × E_nat ≈ constant. Lifetime and natural emissions trade off, so the data alone cannot fix both. This is why the later models use the literature lifetime as a prior.
- **Top-down mass balance** (L11 mass conservation used backwards): E = β(dC/dt + C/τ) with τ = 9.1 yr.
  - Total emissions for 2010–2019: **574 Tg/yr**, against GCP's independent top-down 575 Tg/yr.
  - The part not in EDGAR (natural + missing) was 289 Tg/yr (1985–99), **268 (2000–06)**, 267 (2007–19) and **291 (2020–22)**.
- **Spin-up from 1750** (as in the L11 Rodhe model):
  - With a natural source that balances 729 ppb (219 Tg/yr), the model gives 839 ppb in 1850 (observed 808) but only 1760 ppb in 2025 (observed 1936).
  - With the 1984–2025 fitted natural source (274 Tg/yr), it matches 2025 but is far too high before 1950.
  - So emissions outside the inventory are today about 55 Tg/yr above the pre-industrial natural level: missing inventory sources, or natural emissions that have risen.

Figures: `02_one_box_fit.png`, `02_cost_surface.png`, `02_topdown_emissions.png`, `02_spinup_1750.png`.

---

## 3. The CH₄–OH system (step 3; L15 mass action, L04 stability/feedback/tipping, L08 Jacobian and phase plane, L13 Michaelis–Menten)

Mass action (L15) for CH₄ + OH → products, with OH made at rate P and destroyed by CH₄ and by CO/VOCs:

    dC/dt = E/β − k X C − k_o C
    dX/dt = [P − X (a C/C_ref + 1 − a)] / τ_X        (X = OH / OH_ref, a = 0.31, τ_X ≈ 1 s)

- **Fast and slow variables.**
  - At today's state the Jacobian has eigenvalues −3.2×10⁷ yr⁻¹ (**0.98 s**, OH) and −0.079 yr⁻¹ (**12.7 yr**, CH₄).
  - OH therefore settles almost instantly, so we can set dX/dt = 0. This gives the CH₄ sink k C/(a C/C_ref + 1 − a), a **saturating Michaelis–Menten / Holling type II** form (L13).
  - The reduced one-equation model differs from the full two-variable system by 4×10⁻⁶ ppb.
- **Feedback factor (L04).** The global lifetime is τ = C/S(C) and the perturbation lifetime is τ_p = 1/S′(C). The perturbation lifetime plays the role of L04's response time.

  | Year | τ | τ_p | f = τ_p/τ |
  |---|---|---|---|
  | 2010 | 9.15 yr | 12.25 yr | 1.34 |
  | 2025 | 9.33 yr | 12.69 yr | **1.36** |

  The literature range is 1.3–1.4. OH in 2025 is 2.3 % lower than in 2010 because of methane alone.
- **Pulse test.** A 100 Tg pulse decays with an e-folding time of 12.7 yr with the feedback and 9.15 yr without. The slow eigenvalue, the analytic τ_p and the pulse experiment all agree.
- **Equilibria and tipping (L04, compare L13's washout threshold).**
  - The sink is monotonic, so there is exactly **one stable equilibrium** for every emission level. There is no bistability or hysteresis, unlike the L04 snowball Earth.
  - Holding today's 575 Tg/yr gives 1955 ppb with the feedback and 1913 without.
  - **If OH were the only sink, no equilibrium would exist above E_crit = β k C_ref/a = 1425 Tg/yr (2.5× today)**: OH would be used up and methane would run away. With 20 % less OH source (e.g. cleaner air with less NOx) the threshold drops to 1140 Tg/yr.
  - Stratospheric, soil and Cl sinks grow linearly with C and remove the asymptote. Equilibrium methane then rises steeply but stays finite.
- **Phase portrait (L08 style).** OH is slowed to τ_X = 0.5 yr for display. Trajectories first jump vertically onto the OH nullcline (fast), then slide along it to the stable node (slow).

Figures: `03_pulse_lifetime.png`, `03_equilibrium_threshold.png`, `03_phase_portrait.png`.

---

## 4. Confronting the model hierarchy with data (step 4; L02 estimation, L04 sensitivity)

Every variant is fitted to NOAA annual means for 1984–2025 by least squares. Where k_OH is free, an extra term keeps it near the literature value τ_OH = 11.2 ± 1.3 yr. "Holdout" means: fit 1984–2012, then predict 2013–2025.

| Model | Free params | τ (yr) | RMSE (ppb) | AIC | Growth r | Holdout RMSE |
|---|---|---|---|---|---|---|
| V1 one-box, lifetime fixed | 1 | 9.15 | 21.0 | 257.8 | 0.37 | 34.7 |
| V1′ one-box, lifetime free (no prior) | 2 | 6.43 | 12.6 | 216.8 | 0.48 | 17.6 |
| V2 + CH₄–OH feedback, lifetime fixed | 1 | 9.15 | 32.2 | 293.5 | 0.20 | 62.6 |
| V3 + NOx/CO drivers of OH | 4 | 7.59 | 3.1 | 102.9 | 0.79 | 14.2 |
| **V4 + temperature-dependent OH rate (best)** | 4 | **8.47** | **3.0** | **100.1** | **0.80** | 6.7 |
| V3′ feedback strength a fitted (a = 0.07 ± 0.14) | 5 | 9.26 | 3.0 | 103.5 | 0.79 | 4.5 |

**Best model V4:**
- k_OH = 0.098 ± 0.004 yr⁻¹; E_nat = 306 ± 19 Tg/yr;
- **OH ∝ E_NOx^0.33±0.02 × E_CO^−0.19±0.02**, × Arrhenius e^{1775ΔT/T²}.

Running the same fits with the CEDS + fires inventory changes RMSE by less than 0.3 ppb, so the conclusions do not depend on the inventory.

What the table says:
- **The CH₄–OH feedback alone makes things worse (V2).** It amplifies the response to steadily rising emissions, so it cannot produce a stall.
- **Externally driven OH is what reproduces the stall–restart–surge shape (V3, V4).** NOx is OH's main source and CO its main competing sink. The fitted NOx sensitivity (0.33) is probably larger than chemistry models give. It likely also absorbs other processes that happen to correlate with NOx, so it should be read as an *effective* sensitivity.
- **The feedback strength a is not identifiable from concentrations.** It stays at the literature value 0.31.

**Attribution** (V4, each driver held at its pre-event value; contributions to mean growth in ppb/yr):

| Period | Observed | Model | Anthropogenic emissions | NOx (OH source) | CO (OH sink) | Temperature |
|---|---|---|---|---|---|---|
| Stall 2000–06 | 0.4 | 1.8 | +2.9 | **−1.6** | −0.1 | −0.3 |
| Restart 2007–13 | 5.5 | 4.9 | +3.3 | +0.9 | −2.5 | 0.0 |
| Surge 2020–21 | 14.1 | 11.8 | −0.7 | **+3.4** | −2.0 | 0.0 |

Model OH was +2.8 % higher in 2000–06 than in 1994–99, −0.9 % lower in 2020 than in 2019, and −3.7 % lower in 2025 than in 2010.

- **Stall:** the model still over-predicts growth (1.8 vs 0.4 ppb/yr). Combined with the non-inventory emissions falling by about 20 Tg/yr (section 2), this points to natural or unaccounted emissions dipping during 2000–06 while inventories kept growing.
- **Surge:** the model gets about half of the jump from 8 to 14 ppb/yr, mostly through lower OH after the COVID NOx drop. The rest needs extra emissions, matching the +25 Tg/yr in the mass balance.

Figures: `04_model_hierarchy.png`, `04_attribution.png`.

---

## 5. Two-hemisphere model (step 5; the L11 two-box model with exchange)

The model is dN/dt = 2fE/β − kN − (N − S)/τ_ex, and the same for S. The gradient then obeys dG/dt = 2(2f − 1)E/β − (k + 2/τ_ex)G.

- **Exchange time from a tracer.** Methyl chloroform (MCF) was an industrial solvent emitted almost only in the north. Fitting the same two-box model to its NH and SH records gives **τ_ex = 1.03 ± 0.09 yr** (0.99–1.09 for an NH emission share of 0.95–1.0). The MCF lifetime is 5.0 yr, and the fit is within 1.2 % RMS.
- **CH₄ hemispheres.** Built from the 23 stations by area-weighted interpolation in sin(latitude). The resulting global mean matches NOAA within 0.3 ppb on average (2.6 ppb at most).
  - The NH–SH gradient is 94 ppb (1990), 78 (2000), 84 (2010), 94 (2020) and 91 (2024). It is noisy year to year because station coverage changes.
- **Share of all methane emitted in the NH:** 0.73 (1986–99), 0.72 (2000–06), 0.71 (2007–19), 0.71 (2020–24); range 0.68–0.78 over the τ_ex uncertainty.
  - EDGAR's anthropogenic emissions are 82–86 % NH.
  - So the natural + missing remainder must be only about **60 % NH**, i.e. strongly tropical (tropical wetlands are the obvious candidate).

Figures: `05_mcf_gradient.png`, `05_nh_share.png`.

---

## 6. Forcing and warming (step 6; L03 greenhouse effect, L04 energy balance)

- **CH₄ forcing in 2025** (official ClimateIndicator ERF: **0.570 W/m²**):

  | Formula | 2025 value (W/m²) | RMS error 1750–2025 |
  |---|---|---|
  | Course formula 0.036(√M − √M₀) | 0.612 (+7 %) | 0.016 |
  | Myhre 1998, with N₂O overlap | 0.536 | 0.017 |
  | **Etminan 2016 × AR6 adjustment** | **0.567** | **0.002** |

  The square-root law is the methane analogue of L04's ln(CO₂/280) term. The overlap with N₂O bands is what the course formula leaves out.
- **Energy balance model** (fitted to HadCRUT5 1850–2025 with total ERF):

  | Model | C (W yr m⁻² K⁻¹) | λ (W/m²/K) | ECS | RMS error |
  |---|---|---|---|---|
  | One-layer (L04) | 16.4 | 1.78 | 2.2 K | 0.12 K |
  | **Two-layer with deep ocean** | 16.7 | **1.25** | **3.1 K** | 0.12 K |

  The two-layer model is L11's two-box idea applied to heat, with the ocean coupling γ = 0.67 from the literature and C_d = 100. It agrees closely with the **L04 Budyko slope B = 1.3** and the AR6 best ECS of 3 K. The one-layer model over-estimates λ because it treats ocean heat uptake as feedback.
- **Warming in 2025:** the model gives 1.50 K against 1.43 K observed. **Methane's direct forcing caused 0.28 K** (0.30 K with stratospheric water vapour), about 19 % of the total.

Figure: `06_forcing_ebm.png`.

---

## 7. Futures and the Global Methane Pledge (step 7; L11/L19 scenarios, L04 feedback)

The best model (V4) is driven by SSP emission ratios for CH₄, NOx and CO, and by IPCC AR6 assessed warming. Warming differences come from the two-layer EBM.

| Scenario | CH₄ 2030 | CH₄ 2050 | CH₄ 2100 | Peak | ΔT vs SSP2-4.5 in 2050 (K) |
|---|---|---|---|---|---|
| SSP1-1.9 | 1943 | 1860 | 1651 | 1951 (2028) | −0.022 |
| SSP1-2.6 | 1954 | 1848 | 1627 | 1955 (2029) | −0.023 |
| **SSP2-4.5 (baseline)** | 1978 | **2005** | 1759 | 2019 (2042) | 0 |
| SSP3-7.0 | 1988 | 2222 | 2958 | still rising | +0.022 |
| SSP3-7.0-lowNTCF | 1917 | 1645 | 1594 | 1947 (2026) | −0.054 |
| SSP5-8.5 | 1983 | 2266 | 2304 | 2411 (2074) | +0.026 |
| **Global Methane Pledge** | 1891 | **1674** | 1639 | 1945 (2026) | **−0.054** |
| Emissions frozen at 2025 | 1974 | 2061 | 2095 | 2111 (2078) | +0.004 |
| All anthropogenic emissions stop | 1532 | 910 | 835 | — | −0.206 |

All concentrations are in ppb.

- **Validation:** our SSP2-4.5 (2005 ppb in 2050, 1759 in 2100) is close to MAGICC7's (2020, 1683).
- **Recovery time:** if emissions stopped, methane would fall with an e-folding time of **9.9 yr**, towards about 835 ppb.
- **Pledge, Monte Carlo** (200 random draws of the fitted parameters, of a ∈ [0.25, 0.37] and of λ ± 0.25):
  - ΔC in 2050 = **−333 ppb (−352 to −310)**;
  - ΔT in 2050 = **−0.055 K (−0.065 to −0.047)**;
  - ΔT in 2100 = −0.037 K.
- **Comparison with the official claim.** The pledge's ≥ 0.2 K by 2050 is larger because it counts methane's ozone and other indirect effects; ours is the direct effect only.

**Climate–wetland–methane loop (L04 feedback gain):** g = (dT/dF)(dF/dC)(dC*/dE)(dE/dT).

- Terms: 0.80 K per W/m², 3.7×10⁻⁴ W/m² per ppb, 4.6 ppb per Tg/yr.
- dE/dT combines wetlands growing by a factor Q10 per 10 K (165 Tg/yr × ln(Q10)/10) with the Arrhenius speed-up of the OH sink, which is equivalent to −11 Tg/yr per K.

| Q10 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| Loop gain g | −0.015 | 0.000 | 0.009 | 0.016 |

- Runaway (g = 1) would need dE/dT = 736 Tg/yr/K. **There is no methane runaway.**
- Still, warming driven by everything else (mainly CO₂) can raise wetland emissions a lot. Under SSP5-8.5, Q10 = 2 adds **+200 ppb** and Q10 = 3 adds **+350 ppb** by 2100. The response is large; the self-amplification is small.

Figures: `07_scenarios.png`, `07_pledge_feedback.png`.

---

## 8. Optimal abatement with GEKKO (step 8; L19 dynamic optimisation, L18 externality)

- **Control:** u(t) ∈ [0, 0.7], the abated fraction of SSP2-4.5 anthropogenic emissions, applied to the best-fit model.
- **Cost:** linear marginal cost MAC = $4000/t × u, so the yearly cost is 0.5 m u² E.
- **Discount rate:** 2 %/yr.
- **Solver:** GEKKO IMODE 6 with IPOPT, as in the L19 fishing problem.

**(A) Cost–benefit**, minimise ∑ e^{−rt}(cost + D·C), where D prices the damage so that one tonne emitted today costs the social cost of methane (SC). This mirrors the L19 profit integral.

| SC ($/t) | u(2030) | u(2050) | u(2090) | Analytic u* = SC/m | CH₄ 2050 |
|---|---|---|---|---|---|
| 800 | 0.19 | 0.18 | 0.12 | 0.20 | 1786 |
| **1600** | **0.36** | **0.35** | 0.23 | 0.40 | 1590 |
| 3200 | 0.68 | 0.65 | 0.45 | 0.70 (cap) | 1256 |

- The optimum is a near-constant abatement set where the marginal cost equals the shadow price of an emitted tonne, close to the analytic value.
- It tapers to zero at the 2100 horizon: the same end-of-horizon effect as the L19 fishery ("fish out at the end").
- With SC ≈ $1600/t (US EPA 2023), the optimum (≈ 35 %) is close to the pledge's 30 %.

**(B) Target:** the cheapest path that matches the pledge's 2050 concentration (1676 ppb) and stays at or below it afterwards.

- The optimum **delays**: u = 0.07 in 2030, rising to a spike of **0.54 in 2049**, then 0.22 in 2050.
- Methane lives about 12 yr, so cuts made long before the deadline have decayed by 2050, and discounting rewards waiting.
- **The pledge path costs 1.49× the optimum** (present value $1271 bn vs $855 bn) for the same 2050 concentration. But it keeps methane lower *before* 2050.
- **Lesson:** a single-date concentration target creates odd incentives (last-minute action). A cost–benefit objective gives steady action.

The cost numbers are illustrative: they depend on the assumed cost curve. The *shapes* of the paths are the robust result.

Figure: `08_optimal_abatement.png`.

---

## 9. India (step 9; L18 economy and environment)

- **Emissions:** EDGAR 16.1 Tg (1970) → 23.7 (2000) → 28.0 (2020) → **30.5 Tg (2024)**, 9.2 % of the world.
- **2024 sources:**

  | Source | Tg/yr | Share |
  |---|---|---|
  | Livestock | 16.1 | 53 % |
  | Waste | 4.7 | |
  | Rice | 4.0 | |
  | Coal | 3.1 | |
  | Other | 2.2 | |
  | Oil & gas | 0.4 | |

  Agriculture (livestock + rice) is **66 %**. These are the hardest sources to abate, unlike oil and gas leaks.
- **Per person and per GDP:**

  | | 1990 | 2022 |
  |---|---|---|
  | India, kg per person | 25 | 20 |
  | China, kg per person | 31 | 36 |
  | USA, kg per person | 119 | 84 |
  | Brazil, kg per person | 81 | 105 |
  | World, kg per person | 48 | 41 |
  | India, g per $ of GDP | 12.5 | 2.8 |
  | World, g per $ of GDP | 6.0 | 2.5 |

- **Decomposition, 1990→2022** (emissions = population × GDP per person × methane per GDP, the L18 "impact" identity): emissions +33 % = population +65 % × GDP per person +263 % × methane per GDP **−78 %**. India's methane intensity fell fast; its emissions grew because of population and income.
- **A 30 % Indian cut** of 2020 emissions by 2030 is 8.4 Tg/yr, 8.8 % of the global pledge cut. It lowers 2050 CH₄ by **32 ppb** (−0.005 K).

Figure: `09_india.png`.

---

## Conclusions

1. **Why is methane rising again?** Emissions keep growing, and the **OH sink is not constant**: it follows NOx and CO emissions and is weakened by methane itself (feedback factor 1.36).
   - **Stall:** inventory growth was offset by rising OH and by a dip in non-inventory emissions.
   - **Restart:** OH stopped rising while emissions kept growing.
   - **2020–21 surge:** lower OH after the COVID NOx cuts, plus about 25 Tg/yr of extra emissions, which the hemispheric balance places largely in the tropics.
2. **Box models with a few physically motivated terms explain 1984–2025 to 3 ppb.** A constant lifetime cannot, and neither can the methane–OH feedback by itself.
3. **There is no methane tipping point within reach.** But warming from CO₂ can push wetlands to add hundreds of ppb by 2100.
4. **The Global Methane Pledge works fast but modestly.** It gives about −330 ppb and −0.05 K (direct) by 2050, and its benefit lasts only while the cuts last. Economically it sits close to the optimum for a social cost of about $1600/t.
5. **For India**, methane is mostly a farming question (livestock and rice), and its emissions per person are half the world average.

## Limitations and honest caveats

- **One box, annual means.** There is no seasonal cycle, no spatial OH, and no explicit CO chemistry. The NOx and CO exponents are *effective* sensitivities fitted to one record.
- **Natural emissions are a single constant.** Real wetlands vary with ENSO and rainfall, which is why the model under-predicts the 2020–22 surge. The model cannot say which *type* of source grew; it can only say how much and roughly where (north vs south).
- **Inventories differ by 10–15 %.** Results are robust to swapping EDGAR for CEDS, but absolute natural-emission numbers are not.
- **Warming numbers are methane's direct forcing only**, from a two-layer EBM with a literature ocean coupling.
- **Abatement costs are illustrative.**

## Literature values used

| Quantity | Value | Source |
|---|---|---|
| Tg per ppb CH₄ | 2.75 | GCP convention (Saunois et al. 2020, 2025) |
| CH₄ lifetime vs tropospheric OH | 11.2 ± 1.3 yr; total 9.1 yr | Prather, Holmes & Hsu (2012) |
| Stratosphere / soil / Cl partial lifetimes | 120 / 150 / 200 yr | Prather et al. (2012) |
| d ln OH / d ln CH₄ | −0.31 | Holmes et al. (2013); IPCC AR6 WG1 Ch. 6 |
| CH₄ + OH rate constant | 2.45×10⁻¹² e^{−1775/T} | JPL kinetics evaluation |
| CH₄ forcing | Etminan et al. (2016) with −14 % AR6 adjustment | IPCC AR6 Ch. 7 |
| EBM ocean coupling | γ = 0.67 W/m²/K, C_d = 100 | Geoffroy et al. (2013) |
| Future warming by SSP | Table SPM.1 | IPCC AR6 WG1 SPM |
| Social cost of methane | ≈ $1600/t (2 % discount) | US EPA (2023) |

## Suggested report and presentation structure

1. The question: the record, the stall, restart and surge (Fig. 01).
2. The one-box model and why it fails (Figs. 02).
3. The CH₄–OH kinetics: fast–slow reduction, feedback factor, threshold, phase plane (Figs. 03).
4. Confronting the data: the model hierarchy and attribution (Figs. 04).
5. Where the emissions are: the two-hemisphere box model (Figs. 05).
6. Consequences: forcing, warming, scenarios, the pledge, loop gain (Figs. 06–07).
7. Decisions: GEKKO optimal control; India (Figs. 08–09).
8. Conclusions and limitations.

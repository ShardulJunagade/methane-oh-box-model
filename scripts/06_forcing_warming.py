"""Step 6 - From methane to forcing to warming (L03 greenhouse effect, L04 energy-balance model).

1. Methane radiative forcing: the course square-root formula vs Myhre (1998, with N2O overlap) vs
   Etminan (2016) x AR6 adjustment, all checked against the official ERF series (ClimateIndicator).
2. Energy-balance models fitted to HadCRUT5 with the total ERF:
      one-layer (L04):  C dT/dt = F - lam T
      two-layer (deep ocean, like the L11 two-box; gamma from the literature):  C dT/dt = F - lam T - g (T - Td),  Cd dTd/dt = g (T - Td)
3. Warming caused so far by methane.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from analysis import common as cm
from analysis import models as M

print("06 forcing and warming")
conc = cm.ci_concentrations()[["CH4", "N2O"]]
conc = conc.reindex(np.arange(1750, 2026)).interpolate()
erf = cm.ci_erf()
erf.index = erf.index.astype(int)
Mx, Nx = conc.CH4.values, conc.N2O.values
forms = {"Square-root formula 0.036($\\sqrt{M}-\\sqrt{M_0}$)": M.erf_ch4_simple(Mx),
         "Myhre 1998 (with N$_2$O overlap)": M.erf_ch4_myhre(Mx),
         "Etminan 2016 x AR6 adjustment": M.erf_ch4_etminan(Mx, Nx)}
ref = erf.CH4.reindex(conc.index).values
out = {"erf_2025": {"ClimateIndicator": float(ref[-1])}}
for k, v in forms.items():
    out["erf_2025"][k] = float(v[-1])
    out.setdefault("erf_rmse_1750_2025", {})[k] = float(np.sqrt(((v - ref) ** 2).mean()))
print("  ERF 2025:", {k[:12]: round(v, 3) for k, v in out["erf_2025"].items()})

# --- EBM calibration --------------------------------------------------------------------------
agg = cm.ci_erf_aggregates()
F = agg.total.reindex(np.arange(1750, 2026)).values
T_obs = cm.hadcrut()
T_obs = T_obs - T_obs.loc[1850:1900].mean()
yrs = np.arange(1750, 2026)
mask = yrs >= 1850


def two_layer(F, C, lam, g, Cd=100.0, steps=12):
    T, Td, h, out_ = 0.0, 0.0, 1 / steps, []
    for f in F:
        acc = 0.0
        for _ in range(steps):
            dT = (f - lam * T - g * (T - Td)) / C
            dTd = g * (T - Td) / Cd
            T, Td = T + h * dT, Td + h * dTd
            acc += T
        out_.append(acc / steps)
    return np.array(out_)


def rebase(T):
    return T - T[(yrs >= 1850) & (yrs <= 1900)].mean()


r1 = least_squares(lambda v: rebase(M.ebm(F, v[0], v[1]))[mask] - T_obs.values, [8, 1.3], bounds=([1, 0.3], [100, 5]))
# The deep-ocean coupling gamma is taken from the literature (CMIP5 multi-model mean, Geoffroy et
# al. 2013: 0.67 W/m2/K, with Cd = 100 W yr/m2/K); only C and lambda are fitted to the record.
G_FIX = 0.67
r2 = least_squares(lambda v: rebase(two_layer(F, v[0], v[1], G_FIX))[mask] - T_obs.values, [8, 1.2],
                   bounds=([1, 0.3], [50, 5]))
C1, lam1 = r1.x
(C2, lam2), g2 = r2.x, G_FIX
F2X = 3.93  # W/m2 for doubled CO2 (AR6)
out["ebm_one_layer"] = {"C_Wyr_m2K": C1, "lambda_Wm2K": lam1, "rmse_K": float(np.sqrt((r1.fun ** 2).mean())),
                        "ECS_K": F2X / lam1, "course_L04_B": 1.3}
out["ebm_two_layer"] = {"C_Wyr_m2K": C2, "lambda_Wm2K": lam2, "gamma_Wm2K": g2, "Cd_fixed": 100.0,
                        "rmse_K": float(np.sqrt((r2.fun ** 2).mean())), "ECS_K": F2X / lam2,
                        "TCR_like_K": F2X / (lam2 + g2), "gamma_fixed": True}
print(f"  one-layer: C={C1:.1f}, lambda={lam1:.2f} W/m2/K (L04 used B=1.3), ECS={F2X / lam1:.2f} K")
print(f"  two-layer: C={C2:.1f}, lambda={lam2:.2f}, gamma={g2:.2f}, ECS={F2X / lam2:.2f} K")

# --- warming caused by methane ---------------------------------------------------------------
F_ch4 = erf.CH4.reindex(yrs).values
F_ch4_ind = F_ch4 + agg.H2O_stratospheric.reindex(yrs).values
T_ch4 = rebase(two_layer(F_ch4, C2, lam2, g2))
T_ch4_ind = rebase(two_layer(F_ch4_ind, C2, lam2, g2))
T_all = rebase(two_layer(F, C2, lam2, g2))
out["warming_2025_K"] = {"total_model": float(T_all[-1]), "observed": float(T_obs.loc[2025]),
                         "CH4_direct": float(T_ch4[-1]), "CH4_direct_plus_stratH2O": float(T_ch4_ind[-1]),
                         "CH4_share_of_model_warming": float(T_ch4[-1] / T_all[-1])}
print(f"  warming 2025: total {T_all[-1]:.2f} K, methane (direct) {T_ch4[-1]:.2f} K")
cm.save_json(out, "06_forcing_warming")

fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.3))
for i, (k, v) in enumerate(forms.items()):
    ax[0].plot(conc.index, v, color=cm.SERIES[i], label=k)
ax[0].plot(conc.index, ref, color=cm.INK, lw=2.4, ls=":", label="Official ERF (ClimateIndicator)")
ax[0].set_ylabel("CH$_4$ forcing (W/m$^2$)")
ax[0].set_title("Methane radiative forcing")
ax[0].legend(fontsize=6.8)
ax[1].plot(T_obs.index, T_obs.values, color=cm.MUTED, lw=1.0, label="HadCRUT5")
ax[1].plot(yrs[mask], T_all[mask], color=cm.SERIES[0], label="Two-layer EBM, all forcing")
ax[1].plot(yrs[mask], rebase(M.ebm(F, C1, lam1))[mask], color=cm.SERIES[1], lw=1.1, label="One-layer EBM")
ax[1].plot(yrs[mask], T_ch4[mask], color=cm.SERIES[2], label="Two-layer EBM, CH$_4$ only")
ax[1].set_ylabel("Warming vs 1850-1900 (K)")
ax[1].set_title("Energy balance model")
ax[1].legend(fontsize=7)
cm.savefig(fig, "06_forcing_ebm")

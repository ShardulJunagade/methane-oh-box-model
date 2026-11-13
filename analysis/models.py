"""Box models used in the project.

Methane budget (one equation, OH in quasi-steady state):

    dC/dt = E(t)/beta - k_OH * X(C, t) * C - k_other * C

    X = P(t) * A(T) / (a * C/C_ref + 1 - a)          (normalised OH, X = 1 in the reference state)
    P(t) = (NOx/NOx_ref)^gN * (CO/CO_ref)^(-gC)     (OH source / competing-sink modulation)
    A(T) = exp(E_arrh * dT / T_trop^2)              (Arrhenius speed-up of CH4 + OH with warming)
    E(t) = E_anth(t) + E_nat * [(1-w) + w * Q10^(dT/10)]   (w = wetland share of natural sources)

a = 0 gives the classic one-box model with a constant lifetime; a > 0 is the CH4-OH feedback:
the OH sink k*X*C = k*C/(a C/C_ref + 1 - a) is a saturating (Holling type II / Michaelis-Menten)
function of C.  The quasi-steady OH comes from the full two-variable mass-action system
    dX/dt = (P - X (a C/C_ref + 1 - a)) / tau_X,   tau_X ~ 1 s,
see scripts/03_oh_feedback.py.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from analysis import common as cm

C_REF = float(cm.noaa_annual().loc[cm.REF_YEAR])  # ppb, reference CH4 (2010)


@dataclass
class Params:
    k_oh: float = 1 / 11.2      # 1/yr, OH loss frequency at the reference state
    e_nat: float = 200.0        # Tg/yr, natural (+ unaccounted) emissions at the reference state
    a: float = 0.0              # CH4-OH feedback strength (share of OH loss due to CH4)
    gN: float = 0.0             # d ln OH / d ln E_NOx
    gC: float = 0.0             # -d ln OH / d ln E_CO
    q10: float = 1.0            # wetland temperature sensitivity (1 = none)
    arrhenius: bool = False     # temperature dependence of k_OH

    def lifetime(self, C=C_REF):
        x = 1 / (self.a * C / C_REF + 1 - self.a)
        return 1 / (self.k_oh * x + cm.K_OTHER)


def yearly_coeffs(p: Params, years, e_anth: pd.Series, drv: pd.DataFrame):
    """Per-year emission (Tg/yr) and OH multiplier P*A (dimensionless)."""
    years = np.asarray(years)
    d = drv.reindex(years).ffill().bfill()
    dT = d.dT.to_numpy()
    w = cm.WETLAND_SHARE
    E = e_anth.reindex(years).ffill().to_numpy() + p.e_nat * ((1 - w) + w * p.q10 ** (dT / 10))
    P = d.nox.to_numpy() ** p.gN * d.co.to_numpy() ** (-p.gC)
    if p.arrhenius:
        P = P * np.exp(cm.E_ARRH * dT / cm.T_TROP ** 2)
    return E, P


def rhs(C, E, P, p: Params):
    return E / cm.BETA - p.k_oh * P * C / (p.a * C / C_REF + 1 - p.a) - cm.K_OTHER * C


def simulate(p: Params, y0: int, y1: int, C0: float, e_anth: pd.Series, drv: pd.DataFrame,
             steps_per_year: int = 12, E_override=None, P_override=None):
    """RK4 with monthly steps. Emissions and drivers are constant within each calendar year.

    Returns (annual_mean Series indexed y0..y1, Jan-1 values array of length n+1).
    """
    years = np.arange(y0, y1 + 1)
    E, P = yearly_coeffs(p, years, e_anth, drv)
    if E_override is not None:
        E = np.asarray(E_override, float)
    if P_override is not None:
        P = np.asarray(P_override, float)
    h = 1 / steps_per_year
    C = C0
    ann, jan1 = np.empty(len(years)), np.empty(len(years) + 1)
    jan1[0] = C0
    for i in range(len(years)):
        vals = [C]
        for _ in range(steps_per_year):
            k1 = rhs(C, E[i], P[i], p)
            k2 = rhs(C + 0.5 * h * k1, E[i], P[i], p)
            k3 = rhs(C + 0.5 * h * k2, E[i], P[i], p)
            k4 = rhs(C + h * k3, E[i], P[i], p)
            C = C + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            vals.append(C)
        v = np.asarray(vals)
        ann[i] = (0.5 * v[0] + v[1:-1].sum() + 0.5 * v[-1]) / steps_per_year   # trapezoid mean
        jan1[i + 1] = C
    return pd.Series(ann, index=years), jan1


def noaa_jan1(year: int) -> float:
    """Deseasonalised NOAA CH4 at 1 January of `year` (mean of Dec and Jan trend values)."""
    m = cm.noaa_monthly()
    t = m.set_index(["year", "month"]).trend
    return float((t.loc[(year - 1, 12)] + t.loc[(year, 1)]) / 2)


# ---------------------------------------------------------------------------------------------
# Fitting
# ---------------------------------------------------------------------------------------------
K_FIX = 1 / 11.2   # OH loss frequency fixed at the literature value (Prather et al. 2012; total lifetime 9.15 yr)
VARIANTS = {
    # name: (free parameter names, fixed settings).  k_oh free => fitted with the literature lifetime prior.
    "V1 one-box, lifetime fixed": (["e_nat"], {"a": 0.0, "k_oh": K_FIX}),
    "V1' one-box, lifetime free (no prior)": (["k_oh", "e_nat"], {"a": 0.0}),
    "V2 + CH4-OH feedback, lifetime fixed": (["e_nat"], {"a": 0.31, "k_oh": K_FIX}),
    "V3 + NOx/CO drivers of OH": (["k_oh", "e_nat", "gN", "gC"], {"a": 0.31}),
    "V4 + temperature-dependent OH rate": (["k_oh", "e_nat", "gN", "gC"], {"a": 0.31, "arrhenius": True}),
    "V3' feedback strength a fitted": (["k_oh", "e_nat", "gN", "gC", "a"], {}),
}
X0 = {"k_oh": 1 / 11.2, "e_nat": 200.0, "gN": 0.15, "gC": 0.05, "q10": 2.0, "a": 0.3}
BOUNDS = {"k_oh": (0.03, 0.2), "e_nat": (0, 500), "gN": (-1, 2), "gC": (-1, 2), "q10": (0.2, 10), "a": (0, 0.95)}


@dataclass
class Fit:
    name: str
    params: Params
    free: list
    values: np.ndarray
    se: np.ndarray
    years: np.ndarray
    model: pd.Series
    obs: pd.Series
    rss: float
    aic: float
    rmse: float
    growth_r: float
    growth_rmse: float
    cov: object = None

    def summary(self):
        return {"variant": self.name, "fit_years": [int(self.years[0]), int(self.years[-1])],
                "params": {k: [float(v), float(s)] for k, v, s in zip(self.free, self.values, self.se)},
                "lifetime_ref_yr": float(self.params.lifetime()),
                "rmse_ppb": self.rmse, "aic": self.aic, "growth_corr": self.growth_r,
                "growth_rmse_ppb_per_yr": self.growth_rmse}


# Prior on the OH loss frequency: tropospheric OH lifetime of CH4 = 11.2 +- 1.3 yr (Prather et al. 2012,
# from methyl chloroform). Without it the lifetime and the natural source trade off (see 02_one_box).
K_OH_PRIOR = (1 / 11.2, (1 / 11.2) * 1.3 / 11.2)
SIGMA_OBS = 3.0   # ppb, assumed model-data mismatch used to weight data against the prior


def fit_variant(name, free, fixed, y0=1984, y1=2025, inventory="edgar", x0=None,
                prior=True):
    e_anth = cm.anthro_emissions(inventory)
    drv = cm.drivers()
    obs = cm.noaa_annual().loc[y0:y1]
    C0 = noaa_jan1(y0)
    start = np.array([(x0 or X0)[k] for k in free])
    lo = [BOUNDS[k][0] for k in free]
    hi = [BOUNDS[k][1] for k in free]

    def make(v):
        p = Params(**fixed)
        for k, x in zip(free, v):
            setattr(p, k, x)
        return p

    use_prior = prior and "k_oh" in free

    def resid(v):
        ann, _ = simulate(make(v), y0, y1, C0, e_anth, drv)
        r = (ann.values - obs.values) / SIGMA_OBS
        if use_prior:
            r = np.append(r, (v[free.index("k_oh")] - K_OH_PRIOR[0]) / K_OH_PRIOR[1])
        return r

    r = least_squares(resid, start, bounds=(lo, hi), x_scale="jac")
    n, k = len(obs), len(free)
    data_res = r.fun[:n] * SIGMA_OBS
    rss = float((data_res ** 2).sum())
    try:
        # Posterior covariance with sigma_obs rescaled to the actual residual scatter
        s2 = max(rss / max(n - k, 1) / SIGMA_OBS ** 2, 1.0)
        cov = np.linalg.inv(r.jac.T @ r.jac) * s2
        se = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        cov, se = None, np.full(k, np.nan)
    p = make(r.x)
    model, _ = simulate(p, y0, y1, C0, e_anth, drv)
    gm, go = model.diff().dropna(), obs.diff().dropna()
    return Fit(name, p, free, r.x, se, obs.index.values, model, obs, rss,
               float(n * np.log(rss / n) + 2 * k), float(np.sqrt(rss / n)),
               float(np.corrcoef(gm, go)[0, 1]), float(np.sqrt(((gm - go) ** 2).mean())), cov)


# ---------------------------------------------------------------------------------------------
# Radiative forcing and energy balance
# ---------------------------------------------------------------------------------------------
def erf_ch4_simple(M, M0=cm.C_PI):
    """Course/IPCC-TAR square-root formula without the N2O overlap (W/m2)."""
    return 0.036 * (np.sqrt(M) - np.sqrt(M0))


def _overlap(M, N):
    return 0.47 * np.log(1 + 2.01e-5 * (M * N) ** 0.75 + 5.31e-15 * M * (M * N) ** 1.52)


def erf_ch4_myhre(M, M0=cm.C_PI, N0=cm.N2O_PI):
    """Myhre et al. (1998): square-root law minus the CH4-N2O band overlap."""
    return 0.036 * (np.sqrt(M) - np.sqrt(M0)) - (_overlap(M, N0) - _overlap(M0, N0))


def erf_ch4_etminan(M, N, M0=cm.C_PI, N0=cm.N2O_PI, adjustment=-0.14):
    """Etminan et al. (2016) stratospheric-temperature-adjusted forcing x (1 + AR6 adjustment)."""
    Mb, Nb = (M + M0) / 2, (N + N0) / 2
    sarf = (-1.3e-6 * Mb - 8.2e-6 * Nb + 0.043) * (np.sqrt(M) - np.sqrt(M0))
    return sarf * (1 + adjustment)


def ebm(F, heat_capacity, lam, T0=0.0):
    """Zero-dimensional energy balance C dT/dt = F - lam T (L04), exact for yearly-constant F.

    Returns annual-mean temperature for each year of F.
    """
    F = np.asarray(F, float)
    a = np.exp(-lam / heat_capacity)
    T, out = T0, np.empty_like(F)
    for i, f in enumerate(F):
        Teq = f / lam
        # mean over the year of T(t) = Teq + (T - Teq) e^{-lam t / C}
        out[i] = Teq + (T - Teq) * (1 - a) * heat_capacity / lam
        T = Teq + (T - Teq) * a
    return out

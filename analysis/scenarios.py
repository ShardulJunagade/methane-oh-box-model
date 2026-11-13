"""Future driver paths shared by steps 9 and 10.

Anthropogenic CH4, NOx and CO follow each RCMIP/SSP scenario as a ratio to its 2025 value, applied
to our 2025 values. Global temperature follows the IPCC AR6 SPM (Table SPM.1) assessed warming of each
SSP; 2001-2020 was 0.99 K above 1850-1900.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from analysis import common as cm

AR6 = {"ssp119": (1.5, 1.6, 1.4), "ssp126": (1.5, 1.7, 1.8), "ssp245": (1.5, 2.0, 2.7),
       "ssp370": (1.5, 2.1, 3.6), "ssp585": (1.6, 2.4, 4.4)}   # K above 1850-1900: 2021-40, 2041-60, 2081-2100
T_REF_PI = 0.99
FUT = np.arange(2026, 2101)


def temp_path(ssp: str) -> pd.Series:
    """Temperature anomaly relative to 2001-2020 (the OH-model reference), 2026-2100."""
    drv = cm.drivers()
    t30, t50, t90 = AR6[ssp]
    t2025 = float(drv.dT.loc[2025]) + T_REF_PI
    pts = pd.Series({2025: t2025, 2030: t30, 2050: t50, 2090: t90, 2100: t90 + (t90 - t50) / 4})
    return pd.Series(np.interp(FUT, pts.index, pts.values), index=FUT) - T_REF_PI


def scenario(ssp: str, ch4_override: pd.Series | None = None, ntcf: str | None = None):
    """(anthropogenic CH4 emissions 1750-2100, drivers DataFrame 1750-2100) for one scenario."""
    e_hist, drv = cm.anthro_emissions("edgar"), cm.drivers()
    src = ntcf or ssp

    def ratio(var):
        s = cm.rcmip_series(var, src)
        return (s / s.loc[2025]).loc[2026:2100]

    fut = e_hist.loc[2025] * ratio("Emissions|CH4") if ch4_override is None else ch4_override
    e = pd.concat([e_hist, fut])
    fd = pd.DataFrame({"nox": drv.nox.loc[2025] * ratio("Emissions|NOx"),
                       "co": drv.co.loc[2025] * ratio("Emissions|CO"), "dT": temp_path(ssp)})
    return e, pd.concat([drv, fd])


def gmp_path(cut: float = 0.30, target_year: int = 2030) -> pd.Series:
    """Global Methane Pledge: anthropogenic CH4 falls linearly from 2025 to (1-cut) x 2020 level, then held."""
    e = cm.anthro_emissions("edgar")
    goal = (1 - cut) * e.loc[2020]
    f = np.clip((FUT - 2025) / (target_year - 2025), 0, 1)
    return pd.Series(e.loc[2025] + (goal - e.loc[2025]) * f, index=FUT)


def two_layer(F, C, lam, g, Cd=100.0, steps=12):
    """Two-layer energy balance model (annual means), starting from zero anomaly."""
    T, Td, out = 0.0, 0.0, []
    for f in np.asarray(F, float):
        acc = 0.0
        for _ in range(steps):
            T, Td = T + (f - lam * T - g * (T - Td)) / C / steps, Td + g * (T - Td) / Cd / steps
            acc += T
        out.append(acc / steps)
    return np.array(out)

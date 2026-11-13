"""Step 4 - Confront the model hierarchy with the NOAA record and attribute the stall, restart and surge.

L02: parameters estimated from observations, model choice by AIC and out-of-sample prediction.
L04: counterfactual runs switch one driver off at a time (a sensitivity / feedback analysis).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analysis import common as cm
from analysis import models as M

print("04 fit and compare")
obs = cm.noaa_annual()
drv = cm.drivers()
k_fix = M.K_OH_PRIOR[0]

variants = M.VARIANTS
fits, table = {}, []
for name, (free, fixed) in variants.items():
    f = M.fit_variant(name, free, fixed, prior="no prior" not in name)
    # out-of-sample: fit 1984-2012, predict 2013-2025
    h = M.fit_variant(name, free, fixed, y1=2012, prior="no prior" not in name)
    pred, _ = M.simulate(h.params, 1984, 2025, M.noaa_jan1(1984), cm.anthro_emissions("edgar"), drv)
    f_ci = M.fit_variant(name, free, fixed, inventory="ci", prior="no prior" not in name)
    fits[name] = f
    row = f.summary()
    row.update({"n_params": len(free), "holdout_rmse_2013_2025": float(np.sqrt(((pred.loc[2013:] - obs.loc[2013:2025]) ** 2).mean())),
                "rmse_with_CEDS_inventory": f_ci.rmse, "aic_with_CEDS_inventory": f_ci.aic,
                "lifetime_with_CEDS_inventory": float(f_ci.params.lifetime())})
    table.append(row)
    print(f"  {name:45s} tau={f.params.lifetime():5.2f} RMSE={f.rmse:5.1f} AIC={f.aic:6.1f} "
          f"growth r={f.growth_r:.2f} holdout={row['holdout_rmse_2013_2025']:.1f}")

# --- choose the best model --------------------------------------------------------------
# best = lowest AIC among models with a physically constrained lifetime (V1' is the diagnostic fit)
best_name = min((r for r in table if "no prior" not in r["variant"]), key=lambda r: r["aic"])["variant"]
best = fits[best_name]
bp = best.params
cm.save_json({"variant": best_name, "k_oh": bp.k_oh, "e_nat": bp.e_nat, "a": bp.a, "gN": bp.gN, "gC": bp.gC,
              "q10": bp.q10, "arrhenius": bp.arrhenius, "lifetime_2010": float(bp.lifetime()),
              "se": dict(zip(best.free, best.se.tolist())), "free": best.free,
              "cov": np.asarray(best.cov).tolist()}, "04_best_fit")
print(f"  best by AIC: {best_name}")

# --- attribution: hold one driver at its value before the event -----------------------------
e_anth = cm.anthro_emissions("edgar")
C0 = M.noaa_jan1(1984)
full, _ = M.simulate(bp, 1984, 2025, C0, e_anth, drv)


def run_hold(what, year):
    d = drv.copy()
    e = e_anth.copy()
    if what == "anthropogenic emissions":
        e.loc[year:] = e.loc[year - 1]
    elif what == "NOx (OH source)":
        d.loc[year:, "nox"] = d.loc[year - 1, "nox"]
    elif what == "CO (OH sink)":
        d.loc[year:, "co"] = d.loc[year - 1, "co"]
    elif what == "temperature":
        d.loc[year:, "dT"] = d.loc[year - 1, "dT"]
    s, _ = M.simulate(bp, 1984, 2025, C0, e, d)
    return s


def g(s, a, b):
    return float(s.diff().loc[a:b].mean())


events = {"stall 2000-2006": (2000, 2006, 2000), "restart 2007-2013": (2007, 2013, 2007),
          "surge 2020-2021": (2020, 2021, 2020)}
attrib = {}
for ev, (a, b, y) in events.items():
    row = {"observed_growth": g(obs, a, b), "model_growth": g(full, a, b),
           "observed_growth_before": g(obs, y - 4, y - 1), "model_growth_before": g(full, y - 4, y - 1)}
    for what in ["anthropogenic emissions", "NOx (OH source)", "CO (OH sink)", "temperature"]:
        row[f"contribution: {what}"] = g(full, a, b) - g(run_hold(what, y), a, b)
    attrib[ev] = row
    print(f"  {ev}: obs {row['observed_growth']:.1f} model {row['model_growth']:.1f} | " +
          ", ".join(f"{k.split(': ')[1]} {v:+.1f}" for k, v in row.items() if k.startswith("contribution")))

# OH implied by the best model
years = np.arange(1984, 2026)
E, P = M.yearly_coeffs(bp, years, e_anth, drv)
X = P / (bp.a * full.values / M.C_REF + 1 - bp.a)
oh_model = pd.Series(X, index=years)
cm.save_json({"table": table, "attribution": attrib,
              "oh_model_rel": {int(k): float(v) for k, v in oh_model.items()},
              "oh_model_change_pct": {"2000-2006 vs 1994-1999": float(100 * (oh_model.loc[2000:2006].mean() / oh_model.loc[1994:1999].mean() - 1)),
                                      "2020 vs 2019": float(100 * (oh_model.loc[2020] / oh_model.loc[2019] - 1)),
                                      "2025 vs 2010": float(100 * (oh_model.loc[2025] / oh_model.loc[2010] - 1))}},
             "04_fit_compare")

# --- figures --------------------------------------------------------------------------------
show = ["V1 one-box, lifetime fixed", "V2 + CH4-OH feedback, lifetime fixed", "V3 + NOx/CO drivers of OH",
        "V4 + temperature-dependent OH rate"]
fig, ax = plt.subplots(2, 1, figsize=(7.2, 5.8), sharex=True)
go = obs.loc[1984:2025].diff()
ax[0].bar(go.index, go.values, color=cm.GRID, edgecolor=cm.MUTED, lw=0.4, label="Observed (NOAA)")
for i, n in enumerate(show):
    f = fits[n]
    ax[0].plot(f.model.index, f.model.diff().values, color=cm.SERIES[i], lw=1.4,
               label=f"{n.split(' ', 1)[0]} (r = {f.growth_r:.2f})")
    ax[1].plot(f.model.index, (f.model - f.obs).values, color=cm.SERIES[i], lw=1.4,
               label=f"{n} (RMSE {f.rmse:.1f})")
ax[0].set_ylabel("Growth (ppb/yr)")
ax[0].set_title("Growth rate: the model hierarchy against NOAA")
ax[0].legend(loc="upper left", ncols=3, fontsize=7.5)
ax[1].axhline(0, color=cm.MUTED, lw=0.8)
ax[1].set_ylabel("Model - observed (ppb)")
ax[1].legend(loc="lower left", fontsize=7.5)
cm.savefig(fig, "04_model_hierarchy")

fig, ax = plt.subplots(figsize=(7.2, 3.4))
labels = list(attrib)
drivers_ = ["anthropogenic emissions", "NOx (OH source)", "CO (OH sink)", "temperature"]
w = 0.18
for j, dname in enumerate(drivers_):
    vals = [attrib[e][f"contribution: {dname}"] for e in labels]
    ax.bar(np.arange(3) + (j - 1.5) * w, vals, width=w - 0.02, color=cm.SERIES[j], label=dname)
ax.axhline(0, color=cm.MUTED, lw=0.8)
ax.set_xticks(np.arange(3), [f"{e}\nobs {attrib[e]['observed_growth']:.1f}, model {attrib[e]['model_growth']:.1f} ppb/yr"
                             for e in labels], fontsize=7.5)
ax.set_ylabel("Contribution to growth (ppb/yr)")
ax.set_title(f"What changed the growth rate? ({best_name.split(' ')[0]}, drivers held at pre-event values)")
ax.legend(ncols=2, fontsize=7.5)
cm.savefig(fig, "04_attribution")

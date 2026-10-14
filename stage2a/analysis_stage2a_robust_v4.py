"""Stage 2A robustness, model v4 (lk=0.3, wM=0.5) (Python computes -> CSV; R draws).  Seeded and logged."""
import numpy as np, pandas as pd, platform, sys
from scipy.stats import spearmanr
from crossfeed_v4 import DEFAULT, rhs, run, final, y_init, p_only

HI, P0, T_END = 2.0, 0.01, 2000
rng = np.random.default_rng(20261010)

def wins(B0, F0, p, Pe, P0=P0): return final(B0, F0, P0, p, T=T_END)[7] < 0.05 * Pe
def thr(p, mode, P0=P0, lo=1e-3, hi=HI, n=14):
    Pe = p_only(p, P0)[7]
    f = {"pair": lambda x: (x, x), "B": lambda x: (x, 0), "F": lambda x: (0, x)}[mode]
    if Pe < 0.1 or not wins(*f(hi), p, Pe, P0): return np.nan
    for _ in range(n):
        mid = np.sqrt(lo * hi); hi, lo = (mid, lo) if wins(*f(mid), p, Pe, P0) else (hi, mid)
    return hi
def can(p, mode, P0=P0):                      # does this mode ever succeed up to HI?
    Pe = p_only(p, P0)[7]
    if Pe < 0.1: return np.nan
    f = {"pair": (HI, HI), "B": (HI, 0), "F": (0, HI)}[mode]
    return wins(*f, p, Pe, P0)
p0 = dict(DEFAULT, lk=0.3, wM=0.5)

# ---- 8. robustness over random parameter sets ----
logn = ["muB", "KB", "muP", "KP", "muF", "KAc", "betaB", "betaF", "aB", "qF", "gH", "kH", "gOP", "Oin", "aO", "Hin", "KH", "KO"]
N, rows = 300, []
for i in range(N):
    p = dict(p0)
    for k in logn: p[k] = p0[k] * np.exp(rng.normal(0, 0.3))
    p["phi"] = rng.uniform(0.02, 0.4); p["eps"] = rng.uniform(0.0, 0.6); p["rho"] = rng.uniform(0.6, 0.9)
    p["aP"] = rng.uniform(0, 0.5) * p["aB"]; p["Ac_bg"] = rng.uniform(0, 3.0)
    Pe = p_only(p)[7]
    if Pe < 0.1: continue
    tp = thr(p, "pair"); cB, cF = can(p, "B"), can(p, "F")
    cat = "Nobody" if np.isnan(tp) and not (cB or cF) else ("Single can" if (cB or cF) else "Pair only")
    rows.append(dict(id=i, category=cat, thr_pair=tp, **{k: p[k] for k in logn + ["phi", "eps", "rho", "aP", "Ac_bg"]}))
rb = pd.DataFrame(rows); rb.to_csv("s2a_robustness_v4.csv", index=False)
ok = rb[rb.thr_pair.notna()]
res = [(k, *spearmanr(np.log(ok[k] + 1e-9) if k in logn else ok[k], np.log(ok.thr_pair))) for k in logn + ["phi", "eps", "rho", "aP", "Ac_bg"]]
pd.DataFrame(res, columns=["param", "spearman_rho", "p"]).sort_values("spearman_rho", key=abs, ascending=False).to_csv("s2a_sensitivity_v4.csv", index=False)
print(rb.category.value_counts(normalize=True).round(3).to_string())
print("valid sets:", len(rb), "| default thresholds: pair", thr(p0, "pair"), " F", thr(p0, "F"), " B", thr(p0, "B"))

import scipy, numpy, pandas
open("environment_stage2a_v4.txt", "w").write(f"Python {platform.python_version()}\nnumpy {numpy.__version__}\nscipy {scipy.__version__}\npandas {pandas.__version__}\nseed 20261010\nT_END {T_END}  HI {HI}\n")

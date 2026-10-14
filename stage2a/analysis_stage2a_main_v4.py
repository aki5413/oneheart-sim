"""Stage 2A analyses (Python computes -> CSV; R draws).  Seeded and logged."""
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

# ---- 1. time courses ----
rows = []
for name, (b, f) in {"B alone": (0.3, 0), "F alone": (0, 0.3), "B + F pair": (0.3, 0.3)}.items():
    s = run(y_init(p0, b, f, P0), p0, T=300, npts=300)
    for i, t in enumerate(s.t):
        for k, lab in enumerate(["Sugar", "HostCarbon", "Acetate", "Butyrate", "Oxygen", "B", "F", "P", "Breakdown"]):
            rows.append((name, t, lab, max(s.y[k, i], 0)))
pd.DataFrame(rows, columns=["case", "time", "variable", "value"]).to_csv("s2a_timecourse_v4.csv", index=False)

# ---- 2. equilibria and stability (finite-difference Jacobian) ----
def jac(y, p, h=1e-6):
    n = len(y); J = np.zeros((n, n)); f0 = np.array(rhs(0, y, p))
    for j in range(n):
        yy = np.array(y, float); st = h * max(1, abs(yy[j])); yy[j] += st
        J[:, j] = (np.array(rhs(0, yy, p)) - f0) / st
    return J, f0
rows = []
for name, (b, f) in {"Pathogen-dominated": (0, 0), "Pair-dominated": (0.3, 0.3)}.items():
    y = final(b, f, P0, p0, T=8000); J, f0 = jac(y, p0)
    rows.append(dict(state=name, max_abs_rate=np.abs(f0).max(), max_real_eig=np.linalg.eigvals(J).real.max(),
                     **dict(zip(["Sugar", "HostCarbon", "Acetate", "Butyrate", "Oxygen", "B", "F", "P", "Breakdown"], np.maximum(y, 0)))))
pd.DataFrame(rows).to_csv("s2a_equilibria_v4.csv", index=False)

# ---- 3. threshold vs cross-feeding strength (acetate release aB) ----
rows = []
for aB in np.arange(0, 41, 2.0):
    p = dict(p0, aB=aB)
    for m in ["pair", "B", "F"]: rows.append((aB, m, thr(p, m)))
pd.DataFrame(rows, columns=["aB", "mode", "threshold"]).to_csv("s2a_threshold_vs_aB_v4.csv", index=False)

# ---- 4. invasion analysis vs aB (the analytic backbone) ----
def growth_rates(y, p):
    C, H, Ac, Bu, O, B, F, P, M = np.maximum(y, 0)
    cF = p["muF"] * ((1 - p["omega"]) * H / (p["KH"] + H) + p["omega"] * C / (p["KF"] + C) + p["wM"] * M / (p["KM"] + M))
    gF = cF * (p["eps"] + (1 - p["eps"]) * Ac / (p["KAc"] + Ac))
    Ofac = p["phi"] + (1 - p["phi"]) * O / (p["KO"] + O)
    gP = p["muP"] * (p["wc"] * C / (p["KP"] + C) + (1 - p["wc"]) * H / (p["KPH"] + H)) * Ofac
    return dict(lam_P=gP - p["D"], lam_F=gF - p["D"] - p["betaF"] * O, B=B, F=F)
rows = []
for aB in np.arange(0, 41, 1.0):
    p = dict(p0, aB=aB)
    yB = final(2.0, 0, 0, p, T=8000)         # B-only equilibrium (no pathogen)
    yPair = final(2.0, 2.0, 0, p, T=8000)    # pair equilibrium (no pathogen)
    gB_, gPair = growth_rates(yB, p), growth_rates(yPair, p)
    rows.append(dict(aB=aB, lamP_Bonly=gB_["lam_P"], lamF_Bonly=gB_["lam_F"],
                     lamP_pair=gPair["lam_P"], F_pair=gPair["F"], B_pair=gPair["B"]))
pd.DataFrame(rows).to_csv("s2a_invasion_v4.csv", index=False)

# ---- 5. boundary map: background acetate x F's acetate dependence ----
rows = []
for Ac_bg in [0, 0.5, 1, 2, 3, 4, 6, 8, 12]:
    for eps in [0.02, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8]:
        p = dict(p0, Ac_bg=Ac_bg, eps=eps)
        pair, F_, B_ = can(p, "pair"), can(p, "F"), can(p, "B")
        cat = "Nobody" if not pair and not F_ and not B_ else ("Single suffices" if (F_ or B_) else "Pair only")
        rows.append((Ac_bg, eps, cat, thr(p, "pair") if pair else np.nan))
pd.DataFrame(rows, columns=["Ac_bg", "eps", "category", "pair_threshold"]).to_csv("s2a_boundary_map_v4.csv", index=False)

# ---- 6. mechanism ablation (which ingredients are necessary?) ----
variants = {
  "Full model": {},
  "No acetate cross-feeding (aB=0)": dict(aB=0.0),
  "No acetate or breakdown transfer (aB=0, lk=0)": dict(aB=0.0, lk=0.0),
  "No butyrate -> O feedback": dict(gH=0.0),
  "Pathogen indifferent to O": dict(phi=1.0),
  "No oxygen stress on B, F": dict(betaB=0.0, betaF=0.0),
  "Weakly acetate-dependent F (eps=0.5)": dict(eps=0.5),
  "Competition only": dict(aB=0.0, gH=0.0, lk=0.0),
  "Niche route, weak competitors": dict(eps=1.0, aB=0.0, gH=0.0, lk=0.0, wc=0.5),
  "Niche route, strong O2-tolerant competitors": dict(eps=1.0, aB=0.0, gH=0.0, lk=0.0, wc=0.5, betaB=0.0, betaF=0.0, muB=0.5, KB=0.3, muF=0.5, KH=0.3, KF=0.3),
  "Niche route + full machinery": dict(wc=0.5),
}
rows = []
for name, ch in variants.items():
    p = dict(p0, **ch)
    rows.append((name, thr(p, "pair"), thr(p, "F"), thr(p, "B")))
pd.DataFrame(rows, columns=["variant", "pair", "F_alone", "B_alone"]).to_csv("s2a_ablation_v4.csv", index=False)

# ---- 7. basin boundary: threshold vs initial pathogen load ----
rows = []
for P0v in np.logspace(-3, -0.5, 11):
    for m in ["pair", "B", "F"]: rows.append((P0v, m, thr(p0, m, P0=P0v)))
pn = dict(p0, eps=1.0, aB=0.0, gH=0.0, lk=0.0, wc=0.5, betaB=0.0, betaF=0.0, muB=0.5, KB=0.3, muF=0.5, KH=0.3, KF=0.3)
for P0v in np.logspace(-3, -0.5, 11): rows.append((P0v, "pair, niche route", thr(pn, "pair", P0=P0v)))
pd.DataFrame(rows, columns=["P0", "mode", "threshold"]).to_csv("s2a_threshold_vs_P0_v4.csv", index=False)


# ---- 9. critical cross-feeding strength (root finding) ----
from scipy.optimize import brentq
def lamP_pair(aB):
    p = dict(p0, aB=aB); y = final(2.0, 2.0, 0, p, T=8000); return growth_rates(y, p)["lam_P"]
def lamF_Bonly(aB):
    p = dict(p0, aB=aB); y = final(2.0, 0, 0, p, T=8000); return growth_rates(y, p)["lam_F"]
aB_P = brentq(lamP_pair, 3.0, 6.0, xtol=1e-3)     # pair equilibrium becomes uninvadable by P
aB_F = brentq(lamF_Bonly, 6.0, 12.0, xtol=1e-3)   # F can self-assemble into B-only state
pd.DataFrame([("pair_state_uninvadable_by_P", aB_P), ("F_invades_B_only_state", aB_F)], columns=["criterion", "aB_critical"]).to_csv("s2a_critical_v4.csv", index=False)
print("critical aB:", round(aB_P, 3), round(aB_F, 3))

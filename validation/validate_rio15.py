"""Module-level check against Rio15 batch co-culture (B. adolescentis L2-32 + F. prausnitzii S3L/3 on FOS P95).
Model v3 run in BATCH mode (no flow, no oxygen, no pathogen). Published facts used (text of Rio15):
  - co-culture on FOS: butyrate increase ~8 mM after 24 h (S3L/3 + L2-32)
  - F monoculture on FOS grows poorly; acetate supplementation stimulates it
  - acetate falls between 8 and 24 h in co-culture and ends BELOW the B monoculture
  - starch + A2-165: butyrate only 2.9 mM (weaker)
"""
import sys; sys.path.insert(0, "../stage2a")
import numpy as np, pandas as pd
from scipy.integrate import solve_ivp
from crossfeed_v3 import DEFAULT, rhs

def batch_params(**kw):
    p = dict(DEFAULT, D=0.0, Oin=0.0, aO=0.0, gH=0.0, gOP=0.0, betaB=0.0, betaF=0.0, omega=1.0, Hin=0.0, wc=1.0,
             muB=0.28, KB=0.5, YB=0.035, muF=0.20, KF=0.5, YF=0.05, eps=0.05, aB=23.0, qF=14.0, rho=0.8)
    p.update(kw); return p
def run_batch(B0, F0, p, C0=12.3, T=24, n=97):
    y0 = [C0, 0.0, 0.0, 0.0, 0.0, B0, F0, 0.0]
    s = solve_ivp(rhs, [0, T], y0, args=(p,), method="LSODA", t_eval=np.linspace(0, T, n), rtol=1e-8, atol=1e-10)
    return s
rows = []
p = batch_params()
for name, (b, f) in {"B alone": (0.01, 0.0), "F alone": (0.0, 0.01), "B + F": (0.01, 0.01)}.items():
    s = run_batch(b, f, p)
    for i, t in enumerate(s.t):
        rows.append((name, t, s.y[2, i], s.y[3, i], s.y[5, i], s.y[6, i], s.y[0, i]))
df = pd.DataFrame(rows, columns=["case", "time", "acetate", "butyrate", "B", "F", "sugar"]); df.to_csv("rio15_batch_model.csv", index=False)
for c in ["B alone", "F alone", "B + F"]:
    d = df[df.case == c].set_index("time")
    print(f"{c:8s}  t=8h: Ac={d.loc[8,'acetate']:.2f} Bu={d.loc[8,'butyrate']:.2f} | t=24h: Ac={d.loc[24,'acetate']:.2f} Bu={d.loc[24,'butyrate']:.2f} | peak Ac={d.acetate.max():.2f} at {d.acetate.idxmax():.0f} h | sugar24={d.loc[24,'sugar']:.2f}")
# one-parameter calibration: choose qF so that co-culture butyrate at 24 h ~ 8 mM (Rio15), then re-check the qualitative acetate predictions
from scipy.optimize import brentq
def bu24(qF): return run_batch(0.01, 0.01, batch_params(qF=qF)).y[3, -1] - 8.0
try:
    q = brentq(bu24, 1.0, 120.0); print("calibrated qF for 8 mM butyrate at 24 h:", round(q, 1))
    pc = batch_params(qF=q)
    for name, (b, f) in {"B alone": (0.01, 0.0), "B + F (calibrated)": (0.01, 0.01), "F alone": (0.0, 0.01)}.items():
        s = run_batch(b, f, pc); ac = s.y[2]
        print(f"  {name:19s} acetate: peak {ac.max():.1f} mM at {s.t[ac.argmax()]:.0f} h, 8h {ac[32]:.1f}, 24h {ac[-1]:.1f} | butyrate24 {s.y[3,-1]:.2f}")
except ValueError as e:
    print("no calibration solution:", e)
# starch / A2-165-like case: weaker competitor for substrate
ps = batch_params(muF=0.12, KF=2.0, eps=0.02)
s = run_batch(0.01, 0.01, ps); print("starch-like (weak F competitor): butyrate24 =", round(s.y[3, -1], 2), "mM (Rio15: 2.9 mM)")

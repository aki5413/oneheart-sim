"""Rio15 batch check, model v4: representative parameter set (reproduces validation_summary_v4.txt).
Run from this folder:  python rio15_representative.py"""
import sys; sys.path.insert(0, "../stage2a")
from scan_rio15_v4 import bp, rb

def run(**kw):
    p = bp(lk=0.5, muF=0.4, eps=0.05, YF=0.1, omega=0.05, wM=1.0, KM=0.5, **kw)
    co = rb(0.05, 0.05, p); B = rb(0.05, 0, p); F = rb(0, 0.05, p)
    return p, co, B, F

p, co, B, F = run()
ac = co.y[2]
print(f"co-culture butyrate at 24 h: {co.y[3,-1]:.1f} mM  (Rio15: about 8)")
print(f"co-culture acetate: peak {ac.max():.1f} mM at {co.t[ac.argmax()]:.0f} h, {ac[-1]:.1f} mM at 24 h")
print(f"B alone acetate at 24 h: {B.y[2,-1]:.1f} mM;  F alone butyrate: {F.y[3,-1]:.2f} mM")
# starch-like case (Rio15: 2.9 mM): weaker F competitor for substrate (same definition as validate_rio15.py:
# muF=0.12, KF=2.0, eps=0.02), with the other parameters of the representative set
ps = bp(lk=0.5, muF=0.12, eps=0.02, KF=2.0, YF=0.1, omega=0.05, wM=1.0, KM=0.5)
st = rb(0.05, 0.05, ps)
print(f"starch-like case butyrate at 24 h: {st.y[3,-1]:.1f} mM  (Rio15: 2.9)")

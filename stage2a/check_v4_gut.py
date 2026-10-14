"""Do Stage-2A conclusions survive when B releases breakdown products usable by F (model v4)? Threshold = inoculum per species giving P<5% of P-only eq (P0=0.01)."""
import numpy as np
from crossfeed_v4 import DEFAULT, final
def thr(p,mode,lo=1e-3,hi=10.0,n=14):
    Pe=final(0,0,0.01,p,T=2000)[7]
    f={"pair":lambda x:(x,x),"B":lambda x:(x,0),"F":lambda x:(0,x)}[mode]
    w=lambda x: final(*f(x),0.01,p,T=2000)[7]<0.05*Pe
    if not w(hi): return np.nan
    for _ in range(n):
        m=np.sqrt(lo*hi); hi,lo=(m,lo) if w(m) else (hi,m)
    return hi
for lk,wM in [(0,0),(0.3,0.5),(0.3,1.0),(0.5,1.0)]:
    p=dict(DEFAULT,lk=lk,wM=wM)
    print(f"lk={lk} wM={wM}: pair={thr(p,'pair'):.3f}  B alone={thr(p,'B'):.3f}  F alone={thr(p,'F'):.3f}",flush=True)

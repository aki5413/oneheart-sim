"""Matched-total-dose design: at fixed total inoculum T, vary the B fraction f (f=0 F only, f=1 B only, 0<f<1 pair). Outcome = pathogen suppressed (P < 5% of P-only) in
(a) spatial model, (b) well-mixed model v4. Also spatial threshold of the pair, and effect of mucosal P seeding. seed-free, deterministic."""
import sys; sys.path.insert(0,"../stage2a")
import numpy as np, pandas as pd
import spatial_v4 as S, crossfeed_v4 as W
ps=S.SP; pw=dict(W.DEFAULT,lk=0.3,wM=0.5)
Pe_s=S.Ptot(S.final(0,0,0.01,ps)); Pe_w=W.final(0,0,0.01,pw,T=2000)[7]
Ts=[0.1,0.2,0.3,0.5,1,2,4]; fs=[0,0.1,0.25,0.5,0.75,0.9,1]; rows=[]
for T in Ts:
    for f in fs:
        b,fF=T*f,T*(1-f)
        rs=S.Ptot(S.final(b,fF,0.01,ps))/Pe_s; rw=W.final(b,fF,0.01,pw,T=2000)[7]/Pe_w
        rows.append((T,f,rs,rw)); print(f"T={T} Bfrac={f}: spatial P/Peq={rs:.3f}  wellmixed={rw:.3f}",flush=True)
d=pd.DataFrame(rows,columns=["total","Bfrac","P_rel_spatial","P_rel_wellmixed"]); d["win_spatial"]=d.P_rel_spatial<0.05; d["win_well"]=d.P_rel_wellmixed<0.05
d.to_csv("s2b_matched_dose.csv",index=False)
# pair threshold (per species), spatial, uniform seeding
def thr(fun,lo=1e-3,hi=10.0,n=14):
    if not fun(hi): return np.nan
    for _ in range(n):
        m=np.sqrt(lo*hi); hi,lo=(m,lo) if fun(m) else (hi,m)
    return hi
print("spatial pair threshold (per species):",round(thr(lambda x:S.Ptot(S.final(x,x,0.01,ps))<0.05*Pe_s),3))
for Db,Ds in [(0.0005,0.5),(0.05,0.5),(0.005,0.05),(0.005,5.0)]:
    q=dict(ps,Db=Db,Ds=Ds); Pe=S.Ptot(S.final(0,0,0.01,q))
    print(f"Db={Db} Ds={Ds}: pair threshold",round(thr(lambda x:S.Ptot(S.final(x,x,0.01,q))<0.05*Pe),3),"| F alone/ B alone at 4:",S.Ptot(S.final(0,4,0.01,q))/Pe<0.05,S.Ptot(S.final(4,0,0.01,q))/Pe<0.05,flush=True)

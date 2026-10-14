"""Hosts differ: each virtual host gets lognormal (CV~15%) noise on key parameters. Fraction of hosts where the pair suppresses P (P<5% of P-only eq) vs inoculum dose; two pathogen loads.
Compare qualitatively with Les22 (graded fraction colonized vs dilution; P load never varied). seed 20261010"""
import numpy as np, pandas as pd
from crossfeed_v4 import DEFAULT, final
rng=np.random.default_rng(20261010)
keys=["muB","muF","aB","qF","gH","Oin","betaF","Hin"]
doses=[0.03,0.06,0.1,0.15,0.2,0.3,0.5,1.0]; P0s=[0.01,0.1]; N=60
base=dict(DEFAULT,lk=0.3,wM=0.5)
hosts=[dict(base,**{k:base[k]*rng.lognormal(0,0.15) for k in keys}) for _ in range(N)]
Pe=[final(0,0,0.01,h,T=1500)[7] for h in hosts]
rows=[]
for P0 in P0s:
    for d in doses:
        w=sum(final(d,d,P0,h,T=1500)[7]<0.05*pe for h,pe in zip(hosts,Pe))
        rows.append((P0,d,w/N)); print(f"P0={P0} dose={d}: protected {w}/{N}",flush=True)
pd.DataFrame(rows,columns=["P0","dose","frac_protected"]).to_csv("s2a_dose_heterogeneity_v4.csv",index=False)

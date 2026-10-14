import sys; sys.path.insert(0, "."); sys.path.insert(0, "../stage2a")
import itertools, numpy as np, pandas as pd
exec(open("validate_rio15.py").read().split("rows = []")[0])
rows=[]
for muB, muF, KF, eps, YF in itertools.product([0.15,0.28],[0.2,0.4,0.6],[0.1,0.5],[0.05,0.3],[0.05,0.1]):
    p = batch_params(muB=muB, muF=muF, KF=KF, eps=eps, YF=YF)
    co = run_batch(0.01,0.01,p); b = run_batch(0.01,0,p)
    ac = co.y[2]; i8 = 32
    rows.append(dict(muB=muB,muF=muF,KF=KF,eps=eps,YF=YF,Bu24=co.y[3,-1],Ac24_co=ac[-1],Ac24_B=b.y[2,-1],
                     Ac_drop_8_24=ac[i8]-ac[-1], F_alone_Bu=run_batch(0,0.01,p).y[3,-1]))
d = pd.DataFrame(rows); d.to_csv("rio15_scan.csv", index=False)
d["ok"] = (d.Bu24.between(5,11)) & (d.Ac24_co < d.Ac24_B) & (d.Ac_drop_8_24>0) & (d.F_alone_Bu<1.5)
print(d.ok.sum(), "of", len(d), "parameter sets reproduce all 4 qualitative/quantitative facts of Rio15")
print(d.sort_values("Bu24",ascending=False).head(12).round(2).to_string(index=False))
print(d[d.ok].round(2).to_string(index=False))

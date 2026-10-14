import itertools, pandas as pd, numpy as np
from scan_rio15_v4 import bp, rb
rows=[]
for F0,lk,muF,eps,YF in itertools.product([0.0005,0.001,0.003,0.01],[0.3,0.5],[0.25,0.4,0.6],[0.05,0.3],[0.05,0.1]):
    p=bp(lk=lk,muF=muF,eps=eps,YF=YF); co=rb(.01,F0,p); B=rb(.01,0,p); F=rb(0,F0,p); a=co.y[2]
    ok=dict(bu=5<=co.y[3,-1]<=11, below=a[-1]<B.y[2,-1], drop=a[32]-a[-1]>0, Falone=F.y[3,-1]<1.5)
    rows.append(dict(F0=F0,lk=lk,muF=muF,eps=eps,YF=YF,Bu24=co.y[3,-1],Ac8=a[32],Ac24=a[-1],peak=a.max(),tpeak=co.t[a.argmax()],AcB24=B.y[2,-1],**ok,all_ok=all(ok.values())))
d=pd.DataFrame(rows); d.to_csv("rio15_scan_v4b.csv",index=False)
print(len(d),"all ok:",d.all_ok.sum()); print(d[d.all_ok].round(2).head(12).to_string())

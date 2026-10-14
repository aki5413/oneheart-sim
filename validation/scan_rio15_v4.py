"""Rio15 batch check for model v4 (B releases breakdown products M that F can use). seed-free grid scan."""
import sys; sys.path.insert(0,"../stage2a")
import numpy as np, pandas as pd, itertools
from scipy.integrate import solve_ivp
from crossfeed_v4 import DEFAULT, rhs
def bp(**kw):
    p=dict(DEFAULT,D=0.0,Oin=0.0,aO=0.0,gH=0.0,gOP=0.0,betaB=0.0,betaF=0.0,Hin=0.0,wc=1.0,
           muB=0.28,KB=0.5,YB=0.035,KF=0.5,eps=0.05,aB=23.0,qF=14.0,rho=0.8,omega=0.05,wM=1.0,KM=0.5)
    p.update(kw);return p
def rb(B0,F0,p,C0=12.3,T=24):
    s=solve_ivp(rhs,[0,T],[C0,0,0,0,0,B0,F0,0,0],args=(p,),method="LSODA",t_eval=np.linspace(0,T,97),rtol=1e-8,atol=1e-10);return s
def check(p):
    co=rb(.01,.01,p); B=rb(.01,0,p); F=rb(0,.01,p)
    Bu=co.y[3,-1]; acco=co.y[2]; 
    ok=dict(bu=5<=Bu<=11, below=co.y[2,-1]<B.y[2,-1], drop=acco[32]-acco[-1]>0, Falone=F.y[3,-1]<1.5)
    return Bu,acco[32],acco[-1],B.y[2,-1],F.y[3,-1],ok
if __name__=="__main__":
    rows=[]
    for lk,muF,eps,YF,KM,omega in itertools.product([0.1,0.3,0.5],[0.15,0.25,0.4],[0.05,0.3],[0.05,0.1],[0.5,2.0],[0.0,0.05]):
        p=bp(lk=lk,muF=muF,eps=eps,YF=YF,KM=KM,omega=omega)
        Bu,a8,a24,aB24,Fa,ok=check(p)
        rows.append(dict(lk=lk,muF=muF,eps=eps,YF=YF,KM=KM,omega=omega,Bu24=Bu,Ac8=a8,Ac24=a24,AcB24=aB24,Fa_Bu=Fa,**ok,all_ok=all(ok.values())))
    d=pd.DataFrame(rows);d.to_csv("rio15_scan_v4.csv",index=False)
    print("n sets",len(d),"all 4 ok:",d.all_ok.sum());print(d[d.all_ok].head(10).round(2).to_string())
    print(d.sort_values("Bu24",ascending=False).head(5).round(2).to_string()) if d.all_ok.sum()==0 else None

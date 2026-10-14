"""Stage 2B: 1-D radial gut model (layer 0 = lumen centre, layer n-1 = mucosa). Local kinetics = model v4 (crossfeed_v4.py).
Gradients: sugar supplied in lumen; host carbon (mucin) and oxygen supplied at mucosa; flow washout strong in lumen, weak at mucosa
(mucus-associated); solutes diffuse (Ds), bacteria diffuse slowly (Db); host oxygen use fuelled by butyrate acts mostly near the mucosa."""
import sys; sys.path.insert(0,"../stage2a")
import numpy as np
from scipy.integrate import solve_ivp
from crossfeed_v4 import DEFAULT
N=10
SP=dict(DEFAULT, lk=0.3, wM=0.5, Ds=0.5, Db=0.005, Dm=0.01, n=N, Oin_m=2.0)
def lap(X,Dc):
    d=np.zeros_like(X); d[1:-1]=X[:-2]-2*X[1:-1]+X[2:]; d[0]=X[1]-X[0]; d[-1]=X[-2]-X[-1]; return Dc*d
def rhs(t,y,p):
    n=p["n"]; Y=np.maximum(y.reshape(9,n),0.0); C,H,Ac,Bu,O,B,F,P,M=Y
    s=np.linspace(0,1,n); Dl=p["D"]*(1-s)+p["Dm"]*s
    gB=p["muB"]*C/(p["KB"]+C)
    tC=p["omega"]*C/(p["KF"]+C); tH=(1-p["omega"])*H/(p["KH"]+H); tM=p["wM"]*M/(p["KM"]+M); tot=tC+tH+tM+1e-12
    cF=p["muF"]*tot; acf=Ac/(p["KAc"]+Ac)
    gF=cF*(p["eps"]+(1-p["eps"])*acf); gFac=cF*(1-p["eps"])*acf
    Ofac=p["phi"]+(1-p["phi"])*O/(p["KO"]+O)
    gPc=p["muP"]*p["wc"]*C/(p["KP"]+C)*Ofac; gPh=p["muP"]*(1-p["wc"])*H/(p["KPH"]+H)*Ofac; gP=gPc+gPh
    upac=p["qF"]*gFac*F; Fuse=gF*F/p["YF"]
    dC=Dl*(p["Cin"]*(1-s)-C)-gB*B/p["YB"]-tC/tot*Fuse-gPc*P/p["YP"]+lap(C,p["Ds"])
    dH=Dl*(p["Hin"]*s-H)-tH/tot*Fuse-gPh*P/p["YP"]+lap(H,p["Ds"])
    dM=-Dl*M+p["lk"]*gB*B/p["YB"]-tM/tot*Fuse+lap(M,p["Ds"])
    dAc=-Dl*Ac+p["aB"]*gB*B+p["aP"]*gP*P-upac+lap(Ac,p["Ds"])
    dBu=-Dl*Bu+upac/p["rho"]+lap(Bu,p["Ds"])
    dO=(p["aO"]*s*(p["Oin_m"]-O)+s*p["gOP"]*p["kOP"]/(p["kOP"]+O)*P-p["gH"]*s*Bu*O/(p["kH"]+O)-0.2*gP*P+lap(O,p["Ds"]))
    dB=(gB-Dl-p["betaB"]*O)*B+lap(B,p["Db"]); dF=(gF-Dl-p["betaF"]*O)*F+lap(F,p["Db"]); dP=(gP-Dl)*P+lap(P,p["Db"])
    return np.concatenate([dC,dH,dAc,dBu,dO,dB,dF,dP,dM])
def y0(p,B0,F0,P0):
    n=p["n"]; s=np.linspace(0,1,n); Y=np.zeros((9,n))
    Y[0]=p["Cin"]*(1-s); Y[1]=p["Hin"]*s; Y[4]=p["Oin_m"]*s; Y[5]=B0; Y[6]=F0; Y[7]=P0
    return Y.ravel()
def final(B0,F0,P0,p,T=2000):
    s=solve_ivp(rhs,[0,T],y0(p,B0,F0,P0),args=(p,),method="LSODA",rtol=1e-7,atol=1e-9)
    return s.y[:,-1].reshape(9,p["n"])
def Ptot(Y): return Y[7].sum()
if __name__=="__main__":
    import time
    p=SP; t=time.time()
    for nm,(b,f) in {"P only":(0,0),"B .3":(.3,0),"F .3":(0,.3),"pair .3":(.3,.3),"pair 1":(1,1)}.items():
        Y=final(b,f,0.01,p); print(f"{nm:8s} P(lumen..mucosa)=",np.round(Y[7][[0,3,6,9]],3)," B=",np.round(Y[5][[0,5,9]],2)," F=",np.round(Y[6][[0,5,9]],2)," O=",np.round(Y[4][[0,5,9]],2)," Bu=",np.round(Y[3][[0,5,9]],2),f"{time.time()-t:.0f}s",flush=True)

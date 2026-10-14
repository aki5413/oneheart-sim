"""
Stage 2A model (v3): cross-feeding pair + host-mediated oxygen feedback (Gui21-style)
State: C sugar, H host-derived carbon (F only), Ac acetate, Bu butyrate, O electron acceptor (oxygen), B, F, P

Key structural choices (all documented in the literature notes):
 - F's growth needs acetate (Rio15); only the acetate-dependent part of growth makes butyrate, and butyrate
   formed is tied to acetate consumed (stoichiometry rho = acetate per butyrate, Mun11/Els14 ~0.6-0.85).
 - Butyrate does NOT inhibit P directly. It fuels host oxygen consumption (Byn17, Gui21): Bu lowers O.
 - P gains from O (respiration) and raises O via inflammation (Gui21). B and F are oxygen-sensitive (F more).
 - Acetate may also come from outside B (Ac_bg, pathogen aP): boundary conditions [Rio15].
Parameters are structural/illustrative, anchored to measured ratios only where stated in comments.
"""
import numpy as np
from scipy.integrate import solve_ivp

DEFAULT = dict(
    D=0.1, Cin=10.0, Oin=2.0, aO=0.3,                     # washout, sugar supply, O supply, O exchange
    muB=0.30, KB=0.5, YB=0.10, betaB=0.02,                 # B: muB ~ B. adolescentis on FOS 0.25-0.29 /h [Rio15]
    muF=0.31, KF=1.5, YF=0.10, betaF=0.06,                 # F: muF ~ 0.31 /h with acetate on FOS (A2-165) [Rio15]; poor sugar competitor (high K)
    eps=0.10, KAc=2.0,                                     # eps: growth without acetate relative to with (strain-dependent 0..0.5) [Rio15]
    muP=0.80, KP=0.5, YP=0.10, phi=0.15, KO=0.3,           # P: fast; phi = fermentative (O-independent) growth fraction
    aB=20.0,                                               # acetate released per unit B biomass formed (~22 mM / ~0.7 OD, Rio15 starch)
    aP=0.0, Ac_bg=0.0,                                     # acetate from P / from background community (post-antibiotic: 0)
    qF=14.0, rho=0.8,                                      # acetate consumed per F biomass (Els14 ~1.4/glucose, Y~0.1); acetate per butyrate (Mun11 0.6-0.83; Els14 0.85)
    gH=0.60, kH=0.3,                                       # host O consumption fuelled by butyrate (per mM Bu)
    gOP=1.0, kOP=1.0,                                      # pathogen-induced O release (inflammation)
    Hin=6.0, KH=1.0, omega=0.0, wc=1.0, KPH=1.0,                            # host-derived carbon for F (mucin glycans/pectin-like; Lopez-Siles 2012 via Rio15); omega = share of F growth from shared sugar
)

def rhs(t, y, p):
    C, H, Ac, Bu, O, B, F, P = np.maximum(y, 0.0)
    gB = p["muB"] * C / (p["KB"] + C)
    cF = p["muF"] * ((1 - p["omega"]) * H / (p["KH"] + H) + p["omega"] * C / (p["KF"] + C))
    acf = Ac / (p["KAc"] + Ac)
    gF = cF * (p["eps"] + (1 - p["eps"]) * acf)
    gF_ac = cF * (1 - p["eps"]) * acf                      # acetate-dependent part of growth
    Ofac = p["phi"] + (1 - p["phi"]) * O / (p["KO"] + O)
    gPc = p["muP"] * p["wc"] * C / (p["KP"] + C) * Ofac
    gPh = p["muP"] * (1 - p["wc"]) * H / (p["KPH"] + H) * Ofac
    gP = gPc + gPh
    D = p["D"]
    up_ac = p["qF"] * gF_ac * F                            # acetate uptake by F
    wC = p["omega"] * C / (p["KF"] + C) / max((1 - p["omega"]) * H / (p["KH"] + H) + p["omega"] * C / (p["KF"] + C), 1e-12)
    dC  = D * (p["Cin"] - C) - gB * B / p["YB"] - wC * gF * F / p["YF"] - gPc * P / p["YP"]
    dH  = D * (p["Hin"] - H) - (1 - wC) * gF * F / p["YF"] - gPh * P / p["YP"]
    dAc = D * (p["Ac_bg"] - Ac) + p["aB"] * gB * B + p["aP"] * gP * P - up_ac
    dBu = -D * Bu + up_ac / p["rho"]
    dO  = (p["aO"] * (p["Oin"] - O) + p["gOP"] * p["kOP"] / (p["kOP"] + O) * P
           - p["gH"] * Bu * O / (p["kH"] + O) - 0.2 * gP * P)  # last term: O used by P respiration
    dB = (gB - D - p["betaB"] * O) * B
    dF = (gF - D - p["betaF"] * O) * F
    dP = (gP - D) * P
    return [dC, dH, dAc, dBu, dO, dB, dF, dP]

def y_init(p, B0, F0, P0):
    return [p["Cin"], p["Hin"], p["Ac_bg"], 0.0, p["Oin"], B0, F0, P0]

def run(y0, p, T=3000, npts=300, dense=False):
    te = np.linspace(0, T, npts)
    return solve_ivp(rhs, [0, T], y0, args=(p,), method="LSODA", t_eval=te, rtol=1e-8, atol=1e-10)

def final(B0, F0, P0, p, T=3000):
    return run(y_init(p, B0, F0, P0), p, T=T, npts=20).y[:, -1]

def p_only(p, P0=0.01):
    return final(0, 0, P0, p)

if __name__ == "__main__":
    p = DEFAULT
    f = p_only(p); print("P-only equilibrium  C=%.2f H=%.2f Ac=%.2f Bu=%.2f O=%.2f B=%.3f F=%.3f P=%.3f" % tuple(f))
    for nm, (b, ff) in {"B alone .3": (.3, 0), "F alone .3": (0, .3), "pair .3": (.3, .3), "pair 1": (1, 1), "pair 3": (3, 3)}.items():
        g = final(b, ff, 0.01, p); print(f"{nm:11s} -> P={g[7]:.3f} B={g[5]:.3f} F={g[6]:.3f} Bu={g[3]:.2f} O={g[4]:.2f} Ac={g[2]:.2f}")

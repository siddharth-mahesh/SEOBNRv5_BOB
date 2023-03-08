import numpy as np
import Flux.auxiliaryfunctions as aux
def compute_hFlm_for_flux(m1, m2, tortoise, q, p, S1, S2, l, m, Hreal, v,chiA,chiS):
    eta = m1*m2/(m1+m2)/(m1+m2)
    
    epsilon = (l+m)%2
    
    r = np.linalg.norm(q)
    
    rholmpowl = aux.rholmpowl(m1,m2,l,m,chiA,chiS,v)
    
    Se_eff = aux.Se_eff(l,m,Hreal,v,q,p,eta)
    
    lfactorialinv = 1/np.math.factorial(l)
    
    Omega = np.power(v,3)
    
    hathatk = m*Hreal*Omega
    
    pihathatk4 = 4*np.pi*hathatk
    
    Tlmprefac = np.sqrt( pihathatk4 / ( 1 - np.exp( -pihathatk4 ) ) )
    
    Tlmprodfac = aux.Tlmprodfac(l,hathatk)
    
    T_lm = lfactorialinv*Tlmprefac*np.sqrt(Tlmprodfac)
    
    vPhi = r*Omega*np.cbrt(aux.vPhiNonKeplerian(m1,m2,tortoise,q,p,S1,S2))
    
    Vl_Phi = np.power(vPhi,l+epsilon)
    
    Ylminuse_minusm = aux.AbsSphericalHarmonicAtPiOver2(l-epsilon,-m)
    
    c_lpluse = aux.Newtonian_c(m1,m2,l,m)
    
    ne_lm = aux.Newtonian_n(m1,m2,l,m)
    
    hNe_lm = eta*ne_lm*c_lpluse*Vl_Phi*Ylminuse_minusm
    
    hF_lm = hNe_lm*Se_eff*T_lm*rholmpowl
    print("l,m = ",l,",",m)
    print("chiS, chiA = ",chiS,",",chiA)
    print("epsilon = ",(l+m)%2)
    print("r = ",np.linalg.norm(q))
    print("rholmpowl = ",aux.rholmpowl(m1,m2,l,m,chiA,chiS,v))
    print("Se_eff = ",aux.Se_eff(l,m,Hreal,v,q,p,eta))
    print("lfactorialinv = ",1/np.math.factorial(l))
    print("Omega = ",np.power(v,3))
    print("hathatk = ",m*Hreal*Omega)
    print("pihathatk4 = ",4*np.pi*hathatk)
    print("Tlmprefac = ",np.sqrt( pihathatk4 / ( 1 - np.exp( -pihathatk4 ) ) ))
    print("Tlmprodfac = ",aux.Tlmprodfac(l,hathatk))
    print("T_lm = ",lfactorialinv*Tlmprefac*np.sqrt(Tlmprodfac))
    print("vPhi = ",r*Omega*np.cbrt(aux.vPhiNonKeplerian(m1,m2,tortoise,q,p,S1,S2)))
    print("Vl_Phi = ",np.power(vPhi,l+epsilon))
    print("Ylminuse_minusm = ",aux.AbsSphericalHarmonicAtPiOver2(l-epsilon,-m))
    print("c_lpluse = ",aux.Newtonian_c(m1,m2,l,m))
    print("ne_lm = ",aux.Newtonian_n(m1,m2,l,m))
    print("hNe_lm = ",eta*ne_lm*c_lpluse*Vl_Phi*Ylminuse_minusm)
    print("hF_lm = ",hNe_lm*Se_eff*T_lm*rholmpowl)
    return hF_lm
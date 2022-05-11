import numpy as np
import Flux.auxiliaryfunctions as aux
def compute_hFlm_for_flux(m1, m2, EMgamma, tortoise, q, p, S1, S2, l, m, Hreal, v,chiA,chiS):
    eta = m1*m2/(m1+m2)/(m1+m2)
    
    epsilon = (l+m)%2
    
    rholmpowl = aux.rholmpowl(m1,m2,l,m,chiA,chiS,v,EMgamma)
    
    Se_eff = aux.Se_eff(l,m,Hreal,v,q,p,eta)
    
    lfactorialinv = 1/np.math.factorial(l)
    
    Omega = np.power(v,3)
    
    hathatk = m*Hreal*Omega
    
    pihathatk4 = 4*np.pi*hathatk
    
    Tlmprefac = np.sqrt( pihathatk4 / ( 1 - np.exp( -pihathatk4 ) ) )
    
    Tlmprodfac = aux.Tlmprodfac(l,hathatk)
    
    T_lm = lfactorialinv*Tlmprefac*np.sqrt(Tlmprodfac)
    
    vPhi = Omega*aux.vPhiNonKeplerian(m1,m2,EMgamma,tortoise,q,p,S1,S2)
    
    Vl_Phi = np.power(vPhi,l+epsilon)
    
    Ylminuse_minusm = aux.AbsSphericalHarmonicAtPiOver2(l-epsilon,-m)
    
    c_lpluse = aux.Newtonian_c(m1,m2,l,m)
    
    ne_lm = aux.Newtonian_n(m1,m2,l,m)
    
    hNe_lm = eta*ne_lm*c_lpluse*Vl_Phi*Ylminuse_minusm
    
    hF_lm = hNe_lm*Se_eff*T_lm*rholmpowl
    
    print("l = ",l," , m = ",m, " , x = %.16e"% (vPhi*vPhi), ",y  = %.16e"% Ylminuse_minusm, ",prefix = %.16e"%(eta*ne_lm*c_lpluse))
    
    return hF_lm
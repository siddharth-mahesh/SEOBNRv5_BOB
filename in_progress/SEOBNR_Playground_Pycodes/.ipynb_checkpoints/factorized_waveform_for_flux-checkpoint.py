import numpy as np
import SEOBNR_Playground_Pycodes.auxiliaryfunctions as aux
def compute_hFlm_for_flux(m1, m2, EMgamma, tortoise, q, p, S1, S2, l, m, Hreal, v,chiA,chiS):
    eta = m1*m2/(m1+m2)/(m1+m2)
    #print("eta = ",eta)
    epsilon = (l+m)%2
    #print("epsilon = ", epsilon)
    rholmpowl = aux.rholmpowl(m1,m2,l,m,chiA,chiS,v,EMgamma)
    #print("rholmpowl = %.16e"%rholmpowl)
    Se_eff = aux.Se_eff(l,m,Hreal,v,q,p,eta)
    #print("Se_eff = %.16e"% Se_eff)
    lfactorialinv = 1/np.math.factorial(l)
    #print("lfactorialinv = ", lfactorialinv)
    Omega = np.power(v,3)
    #print("Omega = ",Omega)
    hathatk = m*Hreal*Omega
    #print("hathatk = ",hathatk)
    pihathatk4 = 4*np.pi*hathatk
    #print("pihathatk4 = ",pihathatk4)
    Tlmprefac = np.sqrt( pihathatk4 / ( 1 - np.exp( -pihathatk4 ) ) )
    #print("Tlmprefac = %.16e"% (Tlmprefac*lfactorialinv) )
    Tlmprodfac = aux.Tlmprodfac(l,hathatk)
    #print("Tlmprodfac = %.16e"%Tlmprodfac)
    T_lm = lfactorialinv*Tlmprefac*np.sqrt(Tlmprodfac)
    #print("T_lm = ", T_lm)
    vPhi = np.linalg.norm(q)*Omega*np.cbrt(aux.vPhiNonKeplerian(m1,m2,EMgamma,tortoise,q,p,S1,S2))
    #print("vPhi = %.16e"%vPhi)
    Vl_Phi = np.power(vPhi,l+epsilon)
    #print("Vl_Phi = ", Vl_Phi)
    Ylminuse_minusm = aux.AbsSphericalHarmonicAtPiOver2(l-epsilon,-m)
    #print("Ylminuse_minusm = ", Ylminuse_minusm)
    c_lpluse = aux.Newtonian_c(m1,m2,l,m)
    #print("c_lpluse = ",c_lpluse)
    ne_lm = aux.Newtonian_n(m1,m2,l,m)
    #print("ne_lm = ", ne_lm)
    hNe_lm = eta*ne_lm*c_lpluse*Vl_Phi*Ylminuse_minusm
    #print("hNewton = %.16e"%hNe_lm)
    hF_lm = hNe_lm*Se_eff*T_lm*rholmpowl
    #print("hF_lm = ",hF_lm)
    return hF_lm
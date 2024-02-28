import numpy as np
from Radiation.v5HM_Flux_unoptimized import v5HM_unoptimized_flux as flux
from Derivatives.v5HM_Hamiltonian_Derivatives_unoptimized import v5HM_unoptimized_omega_circ as dHdpphi_preq0
from Derivatives.v5HM_Hamiltonian_Derivatives_unoptimized import v5HM_unoptimized_dH_dpphi as dHdpphi
from Derivatives.v5HM_Hamiltonian_Derivatives_unoptimized import v5HM_unoptimized_dH_dprstar as dHdprstar
from Derivatives.v5HM_Hamiltonian_Derivatives_unoptimized import v5HM_unoptimized_dH_dr as dHdr
from Hamiltonian.v5HM_Hamiltonian_unoptimized import v5HM_unoptimized_hamiltonian as Hamiltonian
def v5HM_unoptimized_rhs(t,y,params,verbose = False):
    r , phi , prstar , pphi = y[0] , y[1] , y[2] , y[3]
    m1 , m2 , chi1 , chi2 = params[0] , params[1] , params[2] , params[3]
    eta = m1*m2/(m1 + m2)/(m1 + m2)
    H , xi = Hamiltonian(m1,m2,r,prstar,pphi,chi1,chi2)
    Omega_circ = (1/eta)*dHdpphi_preq0(m1,m2,r,pphi,chi1,chi2)
    Omega = (1/eta)*dHdpphi(m1,m2,r,prstar,pphi,chi1,chi2)
    F_phi = flux(m1, m2, r, phi, prstar, pphi, chi1, chi2,Omega,Omega_circ,H)/Omega
    pphidot = F_phi/eta
    prstardot = (1/eta)*(-xi*dHdr(m1,m2,r,prstar,pphi,chi1,chi2) + prstar*F_phi/pphi)
    phidot = Omega
    rdot = (xi/eta)*dHdprstar(m1,m2,r,prstar,pphi,chi1,chi2)
    if not verbose:
        return np.array([rdot,phidot,prstardot,pphidot])
    else:
        return np.array([rdot,phidot,prstardot,pphidot]), F_phi, Omega, Omega_circ, xi, dHdr(m1,m2,r,prstar,pphi,chi1,chi2)/eta, eta, prstar,pphi

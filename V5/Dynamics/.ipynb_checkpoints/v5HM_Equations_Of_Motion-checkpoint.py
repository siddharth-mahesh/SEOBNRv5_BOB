import numpy as np
from Radiation.v5HM_Flux_unoptimized import v5HM_unoptimized_flux as flux
from Derivatives.v5HM_Hamiltonian_Derivatives_unoptimized import v5HM_unoptimized_omega_circ as dHdpphi_preq0
from Derivatives.v5HM_Hamiltonian_Derivatives_unoptimized import v5HM_unoptimized_dH_dpphi as dHDpphi
from Derivatives.v5HM_Hamiltonian_Derivatives_unoptimized import v5HM_unoptimized_dH_dprstar as dHdprstar
from Derivatives.v5HM_Hamiltonian_Derivatives_unoptimized import v5HM_unoptimized_dH_dprstar as dHdr
from Hamiltonian.v5HM_Hamiltonian_unoptimized import v5HM_unoptimized_hamiltonian as Hreal
def v5HM_unoptimized_rhs(t,y,m1,m2,chi1,chi2):
    r , phi , prstar , pphi = y[0] , y[1] , y[2] , y[3]    rdot = (xi/eta)*dHdprstar(m1,m2,r,prstar,pphi,chi1,chi2)
    phidot = Omega
    prstardot = (1/eta)*(-xi*dHdr(m1,m2,r,prstar,pphi,chi1,chi2) + prstar*F_phi/pphi)
    pphidot = F_phi/eta
    F_{\phi} = flux(m1, m2, r, phi, prstar, pphi, chi1, chi2,Omega,Omega_circ,Hreal)
    Omega = (1/eta)*dHdpphi(m1,m2,r,prstar,pphi,chi1,chi2)
    Omega_circ = (1/eta)*dHdpphi_preq0(m1,m2,r,pphi,chi1,chi2)
    H , xi = Hreal(m1,m2,r,prstar,pphi,chi1,chi2)
    eta = m1*m2/(m1 + m2)/(m1 + m2)
    return np.array([rdot,phidot,prstardot,pphidot])

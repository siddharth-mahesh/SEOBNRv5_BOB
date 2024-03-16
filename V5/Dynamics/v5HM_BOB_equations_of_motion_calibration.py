import numpy as np
from Radiation.v5HM_Flux_unoptimized import v5HM_unoptimized_flux as flux
from Derivatives.v5HM_BOB_unoptimized_hamiltonian_first_derivatives_calibration import v5HM_BOB_unoptimized_omega_circ_calibration as dHdpphi_preq0
from Derivatives.v5HM_BOB_unoptimized_hamiltonian_first_derivatives_calibration import v5HM_BOB_unoptimized_dH_dpphi_calibration as dHdpphi
from Derivatives.v5HM_BOB_unoptimized_hamiltonian_first_derivatives_calibration import v5HM_BOB_unoptimized_dH_dprstar_calibration as dHdprstar
from Derivatives.v5HM_BOB_unoptimized_hamiltonian_first_derivatives_calibration import v5HM_BOB_unoptimized_dH_dr_calibration as dHdr
from Hamiltonian.v5HM_BOB_unoptimized_hamiltonian_calibration import v5HM_BOB_unoptimized_hamiltonian_calibration as Hamiltonian
def v5HM_BOB_unoptimized_rhs_calibration(t,y,params,verbose = False):
    r , phi , prstar , pphi = y[0] , y[1] , y[2] , y[3]
    m1 , m2 , chi1 , chi2, a6, dSO = params[0] , params[1] , params[2] , params[3], params[4], params[5]
    eta = m1*m2/(m1 + m2)/(m1 + m2)
    H , xi = Hamiltonian(m1,m2,r,prstar,pphi,chi1,chi2,a6,dSO)
    Omega_circ = (1/eta)*dHdpphi_preq0(m1,m2,r,pphi,chi1,chi2,a6,dSO)
    Omega = (1/eta)*dHdpphi(m1,m2,r,prstar,pphi,chi1,chi2,a6,dSO)
    F_phi = flux(m1, m2, r, phi, prstar, pphi, chi1, chi2,Omega,Omega_circ,H)/Omega
    pphidot = F_phi/eta
    prstardot = (1/eta)*(-xi*dHdr(m1,m2,r,prstar,pphi,chi1,chi2,a6,dSO) + prstar*F_phi/pphi)
    phidot = Omega
    rdot = (xi/eta)*dHdprstar(m1,m2,r,prstar,pphi,chi1,chi2,a6,dSO)
    if not verbose:
        return np.array([rdot,phidot,prstardot,pphidot])
    else:
        return np.array([rdot,phidot,prstardot,pphidot]), F_phi, Omega, Omega_circ, xi, dHdr(m1,m2,r,prstar,pphi,chi1,chi2,a6,dSO)/eta, eta, prstar,pphi

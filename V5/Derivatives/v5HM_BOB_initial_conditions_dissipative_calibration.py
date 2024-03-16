import numpy as np
from Radiation.v5HM_Flux_unoptimized import v5HM_unoptimized_flux as flux
from Derivatives.v5HM_BOB_unoptimized_hamiltonian_second_derivatives_calibration import v5HM_BOB_unoptimized_dH_dr_dpphi_calibration as d2Hdrdpphi
from Derivatives.v5HM_BOB_unoptimized_hamiltonian_second_derivatives_calibration import v5HM_BOB_unoptimized_dH_dr_dr_calibration as d2Hdr2
from Derivatives.v5HM_BOB_unoptimized_hamiltonian_first_derivatives_calibration import v5HM_BOB_unoptimized_omega_circ_calibration as omega_circ
from Derivatives.v5HM_BOB_unoptimized_hamiltonian_first_derivatives_calibration import v5HM_BOB_unoptimized_dH_dpphi_calibration as omega
from Derivatives.v5HM_BOB_unoptimized_hamiltonian_first_derivatives_calibration import v5HM_BOB_unoptimized_dH_dprstar_calibration as dHdprstar
from Hamiltonian.v5HM_BOB_unoptimized_hamiltonian_calibration import v5HM_BOB_unoptimized_hamiltonian_calibration as H
def v5HM_unoptimized_calibration_IC_diss(prstar, params):
    m1, m2, r, pphi, chi1, chi2, a6, dSO = params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7]
    M = m1 + m2
    eta = m1*m2/M/M
    Hreal , xi = H(m1,m2,r,prstar,pphi,chi1,chi2,a6,dSO)
    dHdpr = xi*dHdprstar(m1,m2,r,prstar,pphi,chi1,chi2,a6,dSO)/eta
    Omega = omega(m1,m2,r,prstar,pphi,chi1,chi2,a6,dSO)/eta
    Omega_circ = omega_circ(m1,m2,r,pphi,chi1,chi2,a6,dSO)/eta
    dLdr = -d2Hdr2(m1,m2,r,prstar,pphi,chi1,chi2,a6,dSO)/d2Hdrdpphi(m1,m2,r,prstar,pphi,chi1,chi2,a6,dSO)
    dLdt = flux(m1, m2, r, 0., prstar, pphi, chi1, chi2,Omega,Omega_circ,Hreal)/eta/Omega
    rdot = dLdt/dLdr
    return rdot - dHdpr
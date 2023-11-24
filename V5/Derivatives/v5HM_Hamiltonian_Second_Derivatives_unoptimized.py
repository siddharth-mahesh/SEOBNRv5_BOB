import numpy as np
from ..Radiation.v5HM_Flux_unoptimized import v5HM_unoptimized_flux as flux
from ..Derivatives.v5HM_Hamiltonian_Second_Derivatives_unoptimized import v5HM_unoptimized_dH_dr_dpphi as d2Hdrdpphi
from ..Derivatives.v5HM_Hamiltonian_Second_Derivatives_unoptimized import v5HM_unoptimized_dH_dr_dr as d2Hdrdr2
from ..Derivatives.v5HM_Hamiltonian_Derivatives_unoptimized import v5HM_unoptimized_omega_circ as omega_circ
from ..Derivatives.v5HM_Hamiltonian_Derivatives_unoptimized import v5HM_unoptimized_dH_dpphi as omega
from ..Derivatives.v5HM_Hamiltonian_Derivatives_unoptimized import v5HM_unoptimized_dH_dprstar as dHdprstar
from ..Hamiltonian.v5HM_Hamiltonian_unoptimized import v5HM_unoptimized_hamiltonian as H
def v5HM_unoptimized_IC_diss(prstar, params):
    m1, m2, r, pphi, chi1, chi2 = params[0], params[1], params[2], params[3], params[4], params[5]
    Hreal , xi = H(m1,m2,r,0.,prstar,pphi,chi1,chi2)
    dHdpr = xi*dHdprstar(m1,m2,r,0.,prstar,pphi,chi1,chi2)
    Omega = omega(m1,m2,r,0.,prstar,pphi,chi1,chi2)
    Omega_circ = omega_circ(m1,m2,r,0.,prstar,pphi,chi1,chi2)
    dLdr = -d2Hdr2(m1,m2,r,0.,prstar,pphi,chi1,chi2)/d2Hdrdpphi(m1,m2,r,0.,prstar,pphi,chi1,chi2)
    dLdt = flux(m1, m2, r, 0., prstar, pphi, chi1, chi2,Omega,Omega_circ,Hreal)
    rdot = dLdt/dLdr
    return rdot - dHdpr
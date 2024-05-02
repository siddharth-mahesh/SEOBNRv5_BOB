import numpy as np
from Derivatives.v5HM_BOB_optimized_hamiltonian_first_derivatives import v5HM_BOB_optimized_hamiltonian_first_derivatives
from Derivatives.v5HM_BOB_optimized_hamiltonian_second_derivatives import v5HM_BOB_optimized_hamiltonian_second_derivatives
from Derivatives.v5HM_BOB_optimized_omega import v5HM_BOB_optimized_omega
from Radiation.v5HM_BOB_optimized_flux import v5HM_BOB_optimized_flux
def  v5HM_BOB_optimized_dissipative_initial_conditions(x,params):
    m1 = params[0]
    m2 = params[1]
    chi1 = params[4]
    chi2 = params[5]
    a6 = params[6]
    dSO = params[7]
    r = params[2]
    pphi = params[3]
    prstar = x
    dvalues,hamiltonian_and_xi = v5HM_BOB_optimized_hamiltonian_first_derivatives(m1,m2,r,prstar,pphi,chi1,chi2,a6,dSO)
    ddvalues = v5HM_BOB_optimized_hamiltonian_second_derivatives(m1,m2,r,prstar,pphi,chi1,chi2,a6,dSO)
    Omega_circ = v5HM_BOB_optimized_omega(m1,m2,r,pphi,chi1,chi2,a6,dSO)
    Hreal = hamiltonian_and_xi[0]
    xi = hamiltonian_and_xi[1]
    Omega = dvalues[3]
    flux = v5HM_BOB_optimized_flux(m1, m2, r, prstar, pphi, chi1, chi2, Hreal, Omega, Omega_circ)
    drdt = xi*dvalues[2]
    dLdr_inv = ddvalues[1]/ddvalues[0]
    prstareqn = drdt + dLdr_inv*flux
    return prstareqn


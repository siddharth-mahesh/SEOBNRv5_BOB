import numpy as np"
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
    omega = v5HM_BOB_optimized_omega(m1,m2,r,pphi,chi1,chi2,a6,dSO)
flux = v5HM_BOB_optimized_flux(m1,m2,r,prstar,chi1,chi2,hamiltonian_and_xi[0],dvalues[3],omega)
    drdt = xi*dvalues[2]
    dLdr_inv = ddvalues[1]/ddvalues[0]
    prstareqn = drdt + dLdr_inv*flux
    return prstareqn


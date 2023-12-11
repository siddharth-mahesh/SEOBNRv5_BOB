import numpy as np
from Derivatives.v5HM_Hamiltonian_Derivatives_unoptimized import v5HM_unoptimized_omega_circ as dHdpphi_circ
from Derivatives.v5HM_Hamiltonian_Derivatives_unoptimized import v5HM_unoptimized_dH_dr_circ as dHdr_circ
def v5HM_unoptimized_IC_cons(u , params):
    m1, m2, chi1, chi2, omega = params[0], params[1], params[2], params[3], params[4]
    r, pphi = u[0], u[1]
    eta = m1*m2/((m1 + m2)*(m1 + m2))
    return np.array([dHdr_circ(m1,m2,r,pphi,chi1,chi2), omega - dHdpphi_circ(m1,m2,r,pphi,chi1,chi2)/eta])
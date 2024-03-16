import numpy as np
from Derivatives.v5HM_BOB_unoptimized_hamiltonian_derivatives_calibration import v5HM_BOB_unoptimized_omega_circ_calibration as dHdpphi_circ
from Derivatives.v5HM_BOB_unoptimized_hamiltonian_derivatives_calibration import v5HM_BOB_unoptimized_dH_dr_circ_calibration as dHdr_circ
def v5HM_BOB_unoptimized_IC_cons_calibration(u , params):
    m1, m2, chi1, chi2, a6, dSO, omega = params[0], params[1], params[2], params[3], params[4], params[5], params[6]
    r, pphi = u[0], u[1]
    eta = m1*m2/((m1 + m2)*(m1 + m2))
    return np.array([dHdr_circ(m1,m2,r,pphi,chi1,chi2,a6,dSO), omega - dHdpphi_circ(m1,m2,r,pphi,chi1,chi2,a6,dSO)/eta])
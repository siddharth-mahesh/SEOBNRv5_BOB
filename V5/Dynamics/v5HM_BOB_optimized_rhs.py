import numpy as np
from Radiation.v5HM_BOB_optimized_flux import v5HM_BOB_optimized_flux as flux
from Derivatives.v5HM_BOB_optimized_hamiltonian_first_derivatives import v5HM_BOB_optimized_hamiltonian_first_derivatives as dH
from Derivatives.v5HM_BOB_optimized_omega import v5HM_BOB_optimized_omega as dHdpphi_preq0
def v5HM_BOB_optimized_rhs(t,y,params,verbose = False):
    r , phi , prstar , pphi = y[0] , y[1] , y[2] , y[3]
    m1 , m2 , chi1 , chi2, a6, dSO = params[0] , params[1] , params[2] , params[3], params[4], params[5]
    eta = m1*m2/(m1 + m2)/(m1 + m2)
    hamiltonian_derivs , hamiltonian_and_xi= dH(m1,m2,r,prstar,pphi,chi1,chi2,a6,dSO)
    H = hamiltonian_and_xi[0]
    xi = hamiltonian_and_xi[1]
    dHdr = hamiltonian_derivs[0]
    dHdprstar = hamiltonian_derivs[2]
    dHdpphi = hamiltonian_derivs[3]
    Omega_circ = dHdpphi_preq0(m1,m2,r,pphi,chi1,chi2,a6,dSO)
    Omega = dHdpphi
    F_phi = flux(m1, m2, r, prstar, pphi, chi1, chi2,H,Omega,Omega_circ)
    pphidot = F_phi
    prstardot = (-xi*dHdr + prstar*F_phi/pphi)
    phidot = Omega
    rdot = (xi)*dHdprstar
    if not verbose:
        return np.array([rdot,phidot,prstardot,pphidot])
    else:
        return np.array([rdot,phidot,prstardot,pphidot]), F_phi, Omega, Omega_circ, xi, dHdr(m1,m2,r,prstar,pphi,chi1,chi2,a6,dSO)/eta, eta, prstar,pphi

from __future__ import division
import numpy as np
def compute_d2Hdrdpphi(m1, m2, eta, r, theta, phi, pr, ptheta, pphi, S1r, S1theta, S1phi, S2r, S2theta, S2phi, KK, k0, k1, dSO, dSS, tortoise, EMgamma):
    u = 1/r
    sintheta = np.sin(theta)
    prT = (pr*pr+ptheta*ptheta*u*u+pphi*pphi*u*u/sintheta)
    Hreal = prT
    u_r = -1/r**2
    prT_r = 2*pphi**2*u*u_r/sintheta+2*ptheta**2*u*u_r
    prTprm_pphi = 2*pphi*u**2/sintheta
    Hrealprm_pphi = prTprm_pphi
    prT_rprm_pphi = 4*pphi*u*u_r/sintheta
    Hreal_rprm_pphi = prT_rprm_pphi
    return np.array([Hreal_rprm_pphi])
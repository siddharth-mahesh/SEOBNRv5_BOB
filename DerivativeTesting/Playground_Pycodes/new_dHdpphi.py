from __future__ import division
import numpy as np
def new_compute_dHdpphi(m1, m2, eta, r, theta, phi, pr, ptheta, pphi, S1r, S1theta, S1phi, S2r, S2theta, S2phi, KK, k0, k1, dSO, dSS, tortoise, EMgamma):
    u = 1/r
    sintheta = np.sin(theta)
    prT = (pr*pr+ptheta*ptheta*u*u+pphi*pphi*u*u/sintheta)
    prTprm_pphi = 2*pphi*u**2/sintheta
    Hrealprm_pphi = prTprm_pphi
    return np.array([Hrealprm_pphi])
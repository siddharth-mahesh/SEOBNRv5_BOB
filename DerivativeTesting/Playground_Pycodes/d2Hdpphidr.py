from __future__ import division
import numpy as np
def compute_d2Hdpphidr(m1, m2, eta, r, theta, phi, pr, ptheta, pphi, S1r, S1theta, S1phi, S2r, S2theta, S2phi, KK, k0, k1, dSO, dSS, tortoise, EMgamma):
    u = 1/r
    sintheta = np.sin(theta)
    prT = (pr*pr+ptheta*ptheta*u*u+pphi*pphi*u*u/sintheta)
    prT_pphi = 2*pphi*u**2/sintheta
    uprm_r = -1/r**2
    prTprm_r = 2*pphi**2*u*uprm_r/sintheta + 2*ptheta**2*u*uprm_r
    prT_pphiprm_r = 4*pphi*u*uprm_r/sintheta
    Hreal_pphiprm_r = prT_pphiprm_r
    return np.array([Hreal_pphiprm_r])
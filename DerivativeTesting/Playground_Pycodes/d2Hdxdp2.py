from __future__ import division
import numpy as np
def compute_d2Hdxdp2(m1, m2, eta, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z, KK, k0, k1, dSO, dSS, tortoise, EMgamma):
    Hreal = p1*p1+p2*p2+p3*p3
    Hrealprm_p2 = 2*p2
    return np.array([Hreal,Hreal_xprm_p2])
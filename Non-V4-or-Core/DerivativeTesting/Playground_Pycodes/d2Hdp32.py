from __future__ import division
import numpy as np
def compute_d2Hdp32(m1, m2, eta, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z, KK, k0, k1, dSO, dSS, tortoise, EMgamma):
    Hreal = p1*p1+p2*p2+p3*p3
    Hreal_p3 = 2*p3
    Hrealprm_p3 = 2*p3
    Hreal_p3prm_p3 = 2
    return np.array([Hreal,Hreal_p3prm_p3])
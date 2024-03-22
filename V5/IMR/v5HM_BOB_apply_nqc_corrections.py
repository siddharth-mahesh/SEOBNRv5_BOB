import numpy as np
def v5HM_BOB_apply_nqc_correction(nqc_corrections, h22_full, dynamics):
    ai = nqc_coefficients['a']

    bi = nqc_coefficients['b']

    r = dynamics[:,1]

    omega_orb = dynamics[:,6]

    pr = dynamics[:,3]

    rOmega = r * omega_orb

    h22_mode = h22_full[:,1]

    q1 = pr*pr / (rOmega*rOmega)

    q2 = q1 / r

    q3 = q2 / np.sqrt(r)

    N22amp = 1 + q1*ai[0] + q2*ai[1] + q3*ai[2]

    p1 = -pr / rOmega

    p2 = -p1 * pr * pr

    N22phase = np.exp( 1j*( p1*bi[0] + p2*bi[1] ) )

    N22 = N22amp*N22phase

    h22_inspiral_plunge = h22_mode*N22

    return np.c_[h22_ful[:,0],h22_inspiral_plunge]
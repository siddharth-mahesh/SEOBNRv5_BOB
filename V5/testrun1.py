import numpy as np
import qnm
from IMR.v5HM_BOB_generate_waveform_calibration import v5HM_BOB_generate_waveform_calibration as v5HM
from pyseobnr.generate_waveform import generate_modes_opt

M = 20
f = 20
Msol = 4.925491025543575903411922162094833998e-6
Omega_0 = M*Msol*np.pi*f
q = 1.5
qs = "15"
chi1 = 0.02
chi2 = 0.01
chi1s = "02"
chi2s = "01"

m1 = q/(1. + q)
m2 = 1./(1. + q)
ap = m1*chi1 + m2*chi2
am = m1*chi1 - m2*chi2
nu = m1*m2

dSO = (
    -7.71251231383957 * am ** 3
    - 17.2294679794015 * am ** 2 * ap
    - 238.430383378296 * am ** 2 * nu
    + 69.5461667822545 * am ** 2
    - 10.5225438990315 * am * ap ** 2
    + 362.767393298729 * am * ap * nu
    - 85.8036338010274 * am * ap
    - 1254.66845939312 * am * nu ** 2
    + 472.431937787377 * am * nu
    - 39.742317057316 * am
    - 7.58458103577458 * ap ** 3
    - 42.7601129678844 * ap ** 2 * nu
    + 18.1783435552183 * ap ** 2
    - 201.905934468847 * ap * nu ** 2
    - 90.5790079104259 * ap * nu
    + 49.6299175121658 * ap
    + 478.546231305475 * nu ** 3
    + 679.521769948995 * nu ** 2
    - 177.334831768076 * nu
    - 37.6897780220529
)

para6 = np.array(
    [4.17877875e01, -3.02193382e03, 3.34144394e04, -1.69019140e05, 3.29523262e05]
)
a6 = para6[0] + para6[1] * nu + para6[2] * nu ** 2 + para6[3] * nu ** 3 + para6[4] * nu ** 4

pardTNS = np.array(
    [1.00513217e01, -5.96231800e01, -1.05687385e03, -9.79317619e03, 5.55652392e04]
)
Deltat_NS = nu ** (-1.0 / 5 + pardTNS[0] * nu) * (
    pardTNS[1] + pardTNS[2] * nu + pardTNS[3] * nu ** 2 + pardTNS[4] * nu ** 3
)

Deltat_S = nu ** (-1.0 / 5 + 0 * nu) * (
    8.39238879807543 * am ** 2 * ap
    - 16.9056858928167 * am ** 2 * nu
    + 7.23410583477034 * am ** 2
    + 6.38975598319936 * am * ap ** 2
    + 179.569824846781 * am * ap * nu
    - 40.6063653476775 * am * ap
    + 144.253395844761 * am * nu ** 2
    - 90.1929138487509 * am * nu
    + 14.2203101910927 * am
    - 6.78913884987037 * ap ** 4
    + 5.39962303470497 * ap ** 3
    - 132.224950777226 * ap ** 2 * nu
    + 49.8016443361381 * ap ** 2
    + 384.201018794943 * ap * nu ** 2
    - 141.253181790353 * ap * nu
    + 17.5710132409988 * ap
)

Deltat_v5HM = Deltat_NS + Deltat_S

t_nocalib, h22_nocalib, h22_BOB, h22_inspiral_plunge_NQC, h22_inspiral_plunge_combined, dynamics_BOB = v5HM(M,q,chi1,chi2,f,a6,dSO,Deltat_v5HM,2.4627455127717882e-05,debug = True)

h22_save = np.c_[t_nocalib,np.real(h22_nocalib),-np.imag(h22_nocalib)]
np.savetxt(f"h22_v5HMcalib_BOB_q_{qs}_chi1_{chi1s}_chi2_{chi2s}.dat",h22_save)
np.savetxt(f"h22_v5HMcalib_BOB_inspiral_NQC_q_{qs}_chi1_{chi1s}_chi2_{chi2s}.dat",h22_inspiral_plunge_NQC)
np.savetxt(f"h22_v5HMcalib_BOB_inspiral_no_NQC_q_{qs}_chi1_{chi1s}_chi2_{chi2s}.dat",h22_inspiral_plunge_combined)
np.savetxt(f"v5HMcalib_BOB_inspiral_dynamics_q_{qs}_chi1_{chi1s}_chi2_{chi2s}.dat",dynamics_BOB)
times, modes, model = generate_modes_opt(q,chi1,chi2,Omega_0,debug = True)
v5HM_save = np.c_[times,np.real(modes['2,2']),-np.imag(modes['2,2'])]
np.savetxt(f"h22_v5HM_q_{qs}_chi1_{chi1s}_chi2_{chi2s}.dat",v5HM_save)
np.savetxt(f"v5HM_inspiral_dynamics_q_{qs}_chi1_{chi1s}_chi2_{chi2s}.dat",model.dynamics)


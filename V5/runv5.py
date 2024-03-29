import numpy as np
import qnm
from IMR.v5HM_BOB_generate_waveform_calibration import v5HM_BOB_generate_waveform_calibration as v5HM
from pyseobnr.generate_waveform import generate_modes_opt

M = 20
f = 20
Msol = 4.925491025543575903411922162094833998e-6
Omega_0 = M*Msol*np.pi*f

q = [2,1.5,2]
qs = ["20","15","20"]
chi1 = [.9,0.02,0.1]
chi2 = [.5,0.01,-0.1]
chi1s = ["9","02","1"]
chi2s = ["5","01","m1"]
pert = 1. + 4e-14


for k in range(3):
    
    a6 = 0
    dSO = 0
    Deltat = 'BOB'
    t_nocalib, h22_nocalib = v5HM(M,q[k],chi1[k],chi2[k],f,a6,dSO,Deltat,2.4627455127717882e-05)
    h22_nocalib = np.c_[t_nocalib,np.real(h22_nocalib),-np.imag(h22_nocalib)]    

    times, modes, model = generate_modes_opt(q[k],chi1[k],chi2[k],Omega_0,debug = True)
    #timespert, modespert, modelpert = generate_modes_opt(q[k]*pert,chi1[k]*pert,chi2[k]*pert,Omega_0*pert,debug = True)
    modes_22 = modes['2,2']
    h22_v5HM = np.c_[times,np.real(modes_22),-np.imag(modes_22)]
    
    pyseobnr_h22_label = f"./pyseobnr_h22_q_{qs[k]}_chi1_{chi1s[k]}_chi2_{chi2s[k]}.dat"
    #pyseobnrpert_dynamics_label = f"./pyseobnr_pertO14_dynamics_q_{qs[k]}_chi1_{chi1s[k]}_chi2_{chi2s[k]}.dat"
    #our_dynamics_label = f"./our_dynamics_q_{qs[k]}_chi1_{chi1s[k]}_chi2_{chi2s[k]}.dat"
    nocalib_h22_label = f"./nocalib_h22_q_{qs[k]}_chi1_{chi1s[k]}_chi2_{chi2s[k]}.dat"
    
    np.savetxt(pyseobnr_h22_label,h22_v5HM)
    #np.savetxt(pyseobnrpert_dynamics_label,dynamics_pyseobnr_pert)
    #np.savetxt(our_dynamics_label,dynamics_ours)
    np.savetxt(nocalib_h22_label,h22_nocalib)

"""
    # find window of interpolation pts (i.e exclude points in trusted where np.interp would end up extrapolating)
    times_trusted = dynamics_pyseobnr[:,0]
    times_pert = dynamics_pyseobnr_pert[:,0]
    times_ours = dynamics_ours[:,0]
    times_nocalib = dynamics_nocalib[:,0]

    for i in range(len(dynamics_pyseobnr_pert)-1,1,-1):
        if times_pert[i] < times_ours[-1]:
            right_pert = i
            break
    for i in range(len(dynamics_nocalib)-1,1,-1):
        if times_nocalib[i] < times_ours[-1]:
            right_nocalib = i
            break
    for i in range(len(dynamics_pyseobnr)-1,1,-1):
        if times_trusted[i] < min(times_pert[right_pert],times_nocalib[right_nocalib]):
            right_trusted = i
            break

    dynamics_trusted = dynamics_pyseobnr[:right_trusted+1]
    N = dynamics_trusted.shape[0]
    errs_trusted_pert = np.zeros([N,5])
    errs_trusted_ours = np.zeros([N,5])
    errs_trusted_nocalib = np.zeros([N,5])
    for j in range(N):
        errs_trusted_pert[j,0] = dynamics_trusted[j,0]
        errs_trusted_ours[j,0] = dynamics_trusted[j,0]
        errs_trusted_nocalib[j,0] = dynamics_trusted[j,0]
    for i in range(4):
        prims_ours_interp = np.interp(dynamics_trusted[:,0],dynamics_ours[:,0],dynamics_ours[:,i+1])
        prims_nocalib_interp = np.interp(dynamics_trusted[:,0],dynamics_nocalib[:,0],dynamics_nocalib[:,i+1])
        prims_pert_interp = np.interp(dynamics_trusted[:,0],dynamics_pyseobnr_pert[:,0],dynamics_pyseobnr_pert[:,i+1])
        for j in range(N):
            prim_trusted = dynamics_trusted[j,i+1]
            prim_pert = prims_pert_interp[j]
            prim_ours = prims_ours_interp[j]
            prim_nocalib = prims_nocalib_interp[j]

            errs_trusted_pert[j,i+1] = np.abs( (prim_trusted - prim_pert)/( prim_trusted )  )
            errs_trusted_ours[j,i+1] = np.abs( (prim_trusted - prim_ours)/( prim_trusted )  )
            errs_trusted_nocalib[j,i+1] = np.abs( (prim_trusted - prim_nocalib)/( prim_trusted )  )

    err_trustedpert_label = f"./err_trusted_pert_q_{qs[k]}_chi1_{chi1s[k]}_chi2_{chi2s[k]}.dat"
    err_trustedours_label = f"./err_trusted_ours_q_{qs[k]}_chi1_{chi1s[k]}_chi2_{chi2s[k]}.dat"
    err_trustednocalib_label = f"./err_trusted_nocalib_q_{qs[k]}_chi1_{chi1s[k]}_chi2_{chi2s[k]}.dat"
    np.savetxt(err_trustedpert_label,errs_trusted_pert)
    np.savetxt(err_trustedours_label,errs_trusted_ours)
    np.savetxt(err_trustednocalib_label,errs_trusted_nocalib)
   

from Radiation.v5HM_BOB_unoptimized_merger_ringdown import v5HM_BOB_unoptimized_merger_ringdown
from Dynamics.v5HM_BOB_initial_conditions_calibration import v5HM_BOB_unoptimized_initial_conditions_calibration

M = 20
q = 2
f = 20
Msol = 4.925491025543575903411922162094833998e-6
Omega_0 = M*Msol*np.pi*f

S1 = 0.01
S2 = -0.01
a6 = 0
dSO = 0
Deltat = -1

m1,m2,chi1,chi2,y_init,Omega_0,h,rstop,rISCO,af,Mf,h22NR,omega22NR = v5HM_BOB_unoptimized_initial_conditions_calibration(M,q,S1,S2,f,a6,dSO,Deltat)
if af > 0:
    qnm_cache = qnm.modes_cache(s = -2, l = 2, m = 2, n= 0)
    omega_complex, _, _ = qnm_cache(a = af, interp_only = True)
else:
    qnm_cache = qnm.modes_cache(s = -2, l = 2, m = -2, n= 0)
    omega_complex, _, _ = qnm_cache(a = np.abs(af), interp_only = True)

omega_complex_norm = omega_complex/Mf
omega_qnm = np.real(omega_complex_norm)
tau = -1/(np.imag(omega_complex_norm))
print(tau)
t_0 = 0
times_BOB = np.arange(-30*tau + t_0,30*tau+t_0,0.1)
omega22NR*=-1
print(omega_qnm,omega22NR)
hBOB_amp, hBOB_phase = v5HM_BOB_unoptimized_merger_ringdown(times_BOB,t_0,h22NR,omega22NR,omega_qnm,tau)
times, modes, model = generate_modes_opt(q,chi1,chi2,Omega_0,debug = True)
print(hBOB_amp)

np.savetxt("./tBOB.dat",times_BOB)
np.savetxt("./hBOBamp.dat",hBOB_amp)
np.savetxt("./hBOBphase.dat",hBOB_phase)
h22 = modes['2,2']
np.savetxt("./h22eob.dat",np.abs(h22))
np.savetxt("./teob.dat",times)
    
"""
    
    

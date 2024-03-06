import numpy as np
from Dynamics.v5HM_integrator import v5HM_integrator as v5HM
from Dynamics.v5HM_unoptimized_auxiliary_functions import get_waveforms_inspiral as wf
from pyseobnr.generate_waveform import generate_modes_opt

rng = np.random.default_rng(seed = 50)

M = 20
f = 20
Msol = 4.925491025543575903411922162094833998e-6
Omega_0 = M*Msol*np.pi*f

# Case 1: Random Mass Ratio; Random Aligned positive spins
# Case 2: Random Mass Ratio; Random Aligned negative spins
# Case 3: Random Mass Ratio; Spin1 positive, Spin2 negative
# Case 4: Random Mass Ratio; Spin1 negative, Spin2 positive

q = [2,1.5]
chi1 = [.9,0.02]
chi2 = [.5,0.01]
pert = 1. + 4e-14

for k in range(2):
    dynamics_coarse , dynamics_fine = v5HM(M,q[k],chi1[k],chi2[k],f)
    times, modes, model = generate_modes_opt(q[k],chi1[k],chi2[k],Omega_0,debug = True)
    timespert, modespert, modelpert = generate_modes_opt(q[k]*pert,chi1[k]*pert,chi2[k]*pert,Omega_0*pert,debug = True)
    dynamics_ours = np.vstack((dynamics_coarse,dynamics_fine))
    dynamics_pyseobnr = model.dynamics
    dynamics_pyseobnr_pert = modelpert.dynamics
    pyseobnr_dynamics_label = "./pyseobnr_dynamics_q"+str(q[k])+"_chi1_"+str(int(10*chi1[k]))+"_chi2"+str(int(10*chi2[k]))+".dat"
    pyseobnrpert_dynamics_label = "./pyseobnrpertO14_dynamics_q"+str(q[k])+"_chi1_"+str(int(10*chi1[k]))+"_chi2"+str(int(10*chi2[k]))+".dat"
    our_dynamics_label = "./our_dynamics_q"+str(q[k])+"_chi1_"+str(int(10*chi1[k]))+"_chi2"+str(int(10*chi2[k]))+".dat"
    np.savetxt(pyseobnr_dynamics_label,dynamics_pyseobnr)
    np.savetxt(pyseobnrpert_dynamics_label,dynamics_pyseobnr_pert)
    np.savetxt(our_dynamics_label,dynamics_ours)

    # find window of interpolation pts (i.e exclude points in trusted where np.interp would end up extrapolating)
    times_trusted = dynamics_pyseobnr[:,0]
    times_pert = dynamics_pyseobnr_pert[:,0]
    times_ours = dynamics_ours[:,0]

    for i in range(len(dynamics_pyseobnr_pert)-1,1,-1):
        if times_pert[i] < times_ours[-1]:
            right_pert = i
            break
    for i in range(len(dynamics_pyseobnr)-1,1,-1):
        if times_trusted[i] < times_pert[right_pert]:
            right_trusted = i
            break

    dynamics_trusted = dynamics_pyseobnr[:right_trusted+1]
    N = dynamics_trusted.shape[0]
    errs_trusted_pert = np.zeros([N,5])
    errs_trusted_ours = np.zeros([N,5])
    for j in range(N):
        errs_trusted_pert[j,0] = dynamics_trusted[j,0]
        errs_trusted_ours[j,0] = dynamics_trusted[j,0]
    for i in range(4):
        prims_ours_interp = np.interp(dynamics_trusted[:,0],dynamics_ours[:,0],dynamics_ours[:,i+1])
        prims_pert_interp = np.interp(dynamics_trusted[:,0],dynamics_pyseobnr_pert[:,0],dynamics_pyseobnr_pert[:,i+1])
        for j in range(N):
            prim_trusted = dynamics_trusted[j,i+1]
            prim_pert = prims_pert_interp[j]
            prim_ours = prims_ours_interp[j]

            errs_trusted_pert[j,i+1] = np.abs( (prim_trusted - prim_pert)/( prim_trusted )  )
            errs_trusted_ours[j,i+1] = np.abs( (prim_trusted - prim_ours)/( prim_trusted )  )

    err_trustedpert_label = "./err_trusted_pert_q"+str(q[k])+"_chi1_"+str(int(10*chi1[k]))+"_chi2"+str(int(10*chi2[k]))+".dat"
    err_trustedours_label = "./err_trusted_ours_q"+str(q[k])+"_chi1_"+str(int(10*chi1[k]))+"_chi2"+str(int(10*chi2[k]))+".dat"
    np.savetxt(err_trustedpert_label,errs_trusted_pert)
    np.savetxt(err_trustedours_label,errs_trusted_ours)
    
    
     
            
    
    

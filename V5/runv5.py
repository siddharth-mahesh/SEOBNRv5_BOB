import numpy as np
from Dynamics.v5HM_integrator import v5HM_integrator as v5HM
from Dynamics.v5HM_unoptimized_auxiliary_functions import get_waveforms_inspiral as wf
from pyseobnr.generate_waveform import generate_modes_opt

M = 20
q = 2
chi1 = 0.2
chi2 = 0.1
f = 20
Msol = 4.925491025543575903411922162094833998e-6
Omega_0 = M*Msol*np.pi*f
m1 = q/(1 + q)
m2 = 1/(1 + q)
dynamics_coarse , dynamics_fine = v5HM(M,q,chi1,chi2,f)
times, modes, model = generate_modes_opt(q,chi1,chi2,Omega_0,debug = True)
dynamics_ours = np.vstack((dynamics_coarse,dynamics_fine))
dynamics_pyseobnr = model.dynamics
h22_pyseobnr =np.c_[ times,modes['2,2']]
h22_ours = wf(m1,m2,dynamics_ours,chi1,chi2)


pyseobnr_dynamics_label = "./pyseobnr_dynamics_q"+str(q)+"_chi1_"+str(int(10*chi1))+"_chi2"+str(int(10*chi2))+".dat"
our_dynamics_label = "./our_dynamics_q"+str(q)+"_chi1_"+str(int(10*chi1))+"_chi2"+str(int(10*chi2))+".dat"

pyseobnr_wf_label = "./pyseobnr_wf_q"+str(q)+"_chi1_"+str(int(10*chi1))+"_chi2"+str(int(10*chi2))+".dat"
our_wf_label = "./our_wf_q"+str(q)+"_chi1_"+str(int(10*chi1))+"_chi2"+str(int(10*chi2))+".dat"

np.savetxt(pyseobnr_dynamics_label,dynamics_pyseobnr)
np.savetxt(our_dynamics_label,dynamics_ours)

np.savetxt(pyseobnr_wf_label,h22_pyseobnr)
np.savetxt(our_wf_label,h22_ours)


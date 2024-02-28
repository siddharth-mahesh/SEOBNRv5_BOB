import numpy as np
from Dynamics.v5HM_integrator import v5HM_integrator as v5HM
from pyseobnr.generate_waveform import generate_modes_opt

M = 20
q = 2
chi1 = 0.2
chi2 = 0.1
f = 20
Msol = 4.925491025543575903411922162094833998e-6
Omega_0 = M*Msol*np.pi*f

dynamics_coarse , dynamics_fine = v5HM(M,q,chi1,chi2,f)
times, modes, model = generate_modes_opt(q,chi1,chi2,Omega_0,debug = True)
dynamics_ours = np.vstack((dynamics_coarse,dynamics_fine))
dynamics_pyseobnr = model.dynamics

pyseobnr_label = "./pyseobnr_dynamics_q"+str(q)+"_chi1_"+str(int(10*chi1))+"_chi2"+str(int(10*chi2))+".dat"
our_label = "./our_dynamics_q"+str(q)+"_chi1_"+str(int(10*chi1))+"_chi2"+str(int(10*chi2))+".dat"

np.savetxt(pyseobnr_label,dynamics_pyseobnr)
np.savetxt(our_label,dynamics_ours)


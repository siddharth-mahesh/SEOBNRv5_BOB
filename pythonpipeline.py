# Import core modules
import sys,os
import numpy as np

# define relative error

def E_rel( a , b):
    diff = a - b
    avg = .5*( a + b )
    return np.abs( diff / avg )

# Import v4P Hamiltonian
from V4development.Hamiltonian.v4P_Hamiltonian import v4P_compute_Hamiltonian

# Read in randomized values and Hamiltonian and corresponding perturbed terms as numpy arrays 
fp_main = open(os.path.join("V4development","Hamiltonian","CI-InputsMain.dat"),"rb")
fp_pert = open(os.path.join("V4development","Hamiltonian","CI-InputsPert.dat"),"rb")
fp_hamiltonian = open(os.path.join("V4development","Hamiltonian","CI-Hamiltonian.dat"),"rb")
main_values = np.fromfile(fp_main,np.float64)
pert_values = np.fromfile(fp_pert,np.float64)
ham_values = np.fromfile(fp_hamiltonian,np.float64)
# For counting errors; use sys.exit(errors) to pass/fail the pipeline
mismatch_total = 0
mismatch_gt_1 = 0
mismatch_gt_10 = 0
mismatch_gt_100 = 0
mismatch_gt_1000 = 0

# Go line by line
for i in range(len(main_values)//14):
    # declare LAL's Hreal and perturbed counterpart
    # and calculate tolerance
    Hreal_LAL = ham_values[2*i+0]
    Hreal_LAL_pert = ham_values[2*i+1]
    tol = E_rel( Hreal_LAL , Hreal_LAL_pert)
    
    # declare Hamiltonian inputs
    m1 = main_values[14*i+0]
    m2 = main_values[14*i+1]
    tortoise = 2
    x = main_values[14*i+2]
    y = main_values[14*i+3]
    z = main_values[14*i+4]
    p1 = main_values[14*i+5]
    p2 = main_values[14*i+6]
    p3 = main_values[14*i+7]
    S1x = main_values[14*i+8]
    S1y = main_values[14*i+9]
    S1z = main_values[14*i+10]
    S2x = main_values[14*i+11]
    S2y = main_values[14*i+12]
    S2z = main_values[14*i+13]
    # Call Hamiltonian function
    testham = v4P_compute_Hamiltonian(m1,m2,tortoise,x,y,z,p1,p2,p3,S1x,S1y,S1z,S2x,S2y,S2z)
    #testham = compute_v4P_Hreal(m1,m2,tortoise,x,y,z,p1,p2,p3,S1x,S1y,S1z,S2x,S2y,S2z) 
    #compute relative error and count errors above tolerance
    e_rel = E_rel( testham , Hreal_LAL )
    
    #print(x,y,z,p1,p2,p3,S1x,S1y,S1z,S2x,S2y,S2z)
    if (e_rel > tol):
        mismatch_total += 1
        if (e_rel > 10*tol):
            mismatch_gt_10 += 1
            if (e_rel > 100*tol):
                mismatch_gt_100 += 1
                if (e_rel > 1000*tol):
                    mismatch_gt_1000 += 1
if (mismatch_total > 0):                   
    print("percent mismatches = %.2e"%(100.*mismatch_total/100000.))
    print("percent mismatches > 10 x tolerance = %.2e"%(100.*mismatch_gt_10/100000.))
    print("percent mismatches > 100 x tolerance = %.2e"%(100.*mismatch_gt_100/100000.))
    print("percent mismatches > 1000 x tolerance = %.2e"%(100.*mismatch_gt_1000/100000.))

sys.exit(mismatch_total)   
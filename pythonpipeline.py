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
fp_main = open(os.path.join("V4development","Hamiltonian","CI-HamiltonianTestMain.dat"),"rb")
fp_pert = open(os.path.join("V4development","Hamiltonian","CI-HamiltonianTestPert.dat"),"rb")
main_values = np.from_file(fp_main,np.float64)
pert_values = np.from_file(fp_pert,np.float64)

# For counting errors; use sys.exit(errors) to pass/fail the pipeline
errors = 0

# Go line by line
for i in range(len(main_values)):
    
    # declare LAL's Hreal and perturbed counterpart
    # and calculate tolerance
    Hreal_LAL = main_values[i,14]
    Hreal_LAL_pert = pert_values[i,14]
    tol = E_rel( Hreal_LAL , Hreal_LAL_pert)
    
    # declare Hamiltonian inputs
    m1 = main_values[i,0]
    m2 = main_values[i,1]
    x = main_values[i,2]
    y = main_values[i,3]
    z = main_values[i,4]
    p1 = main_values[i,5]
    p2 = main_values[i,6]
    p3 = main_values[i,7]
    S1x = main_values[i,8]
    S1y = main_values[i,9]
    S1z = main_values[i,10]
    S2x = main_values[i,11]
    S2y = main_values[i,12]
    S2z = main_values[i,13]
    
    # Call Hamiltonian function
    testham = v4P_compute_Hamiltonian(m1,m2,x,y,z,p1,p2,p3,S1x,S1y,S1z,S2x,S2y,S2z)
    
    #compute relative error and count errors above tolerance
    e_rel = E_rel( testham , Hreal_LAL )
    if e_rel > tol:
        errors += 1

sys.exit(errors)   
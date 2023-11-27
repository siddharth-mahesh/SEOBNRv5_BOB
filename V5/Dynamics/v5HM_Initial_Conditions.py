import numpy as np
from scipy.optimize import root, root_scalar
from Derivatives.v5HM_Initial_Conditions_Conservative import v5HM_unoptimized_IC_cons as IC_cons
from Derivatives.v5HM_Initial_Conditions_Dissipative import v5HM_unoptimized_IC_diss as IC_diss
def v5HM_initial_conditions(M,q,S1,S2,f):
    m1 = q/(1 + q) 
    m2 = 1/(1 + q) 
    chi1 = S1 
    chi2 = S2 
    Msol = 4.925491025543575903411922162094833998e-6 
    Omega_start = M*Msol*np.pi*f 
    r_guess = Omega_start**(-2/3) 
    pphi_guess = Omega_start**(-1/3) 
    sol_cons_guess = np.array([r_guess,pphi_guess]) 
    params_cons = np.array([m1,m2,chi1,chi2,Omega_start]) 
    sol_cons = root(IC_cons,sol_cons_guess,args = params_cons,tol = 6e-12) 
    r , pphi = sol_cons.x[0] , sol_cons.x[1] 
    params_diss = np.array([r,pphi,m1,m2,chi1,chi2]) 
    prstar_bracket = [-3e-2,0] 
    prstar_sol = root_scalar(IC_diss, args = params_diss, bracket = prstar_bracket, xtol = 1e-12, rtol = 1e-10) 
    prstar = prstar_sol.root 
    return np.array([r,0.,prstar,pphi])
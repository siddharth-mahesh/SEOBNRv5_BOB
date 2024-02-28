import numpy as np
from scipy.interpolate import CubicSpline
from Derivatives.v5HM_Hamiltonian_Derivatives_unoptimized import v5HM_unoptimized_omega_circ
from Dynamics.v5HM_Equations_Of_Motion import v5HM_unoptimized_rhs
from Dynamics.v5HM_Initial_Conditions import v5HM_unoptimized_initial_conditions 
from Dynamics.v5HM_unoptimized_auxiliary_functions import augment_dynamics, iterative_refinement, interpolate_dynamics 
import pygsl_lite.errno as errno
import pygsl_lite.odeiv2 as odeiv2
_control = odeiv2.pygsl_lite_odeiv2_control
class control_y_new(_control):
    def __init__(self, eps_abs, eps_rel):
        a_y = 1
        a_dydt = 1
        super().__init__(eps_abs, eps_rel, a_y, a_dydt, None)
def v5HM_integrator(M,q,chi1,chi2,f):
     
    m1,m2,chi1,chi2,y_init,Omega_0,h,rstop,risco,Deltat,af,Mf,h22NR,omega22NR = v5HM_unoptimized_initial_conditions(M,q,chi1,chi2,f) 
     
    sys = odeiv2.system(v5HM_unoptimized_rhs,None,4,[m1,m2,chi1,chi2]) 
    T = odeiv2.step_rk8pd 
    s = odeiv2.pygsl_lite_odeiv2_step(T,4) 
    atol = 1e-11 
    rtol = 1e-12 
    c = control_y_new(atol,rtol) 
    e = odeiv2.pygsl_lite_odeiv2_evolve(4) 
     
    prims = [] 
    times = [] 
    omega_previous = Omega_0 
    omega_peak = False 
    prstar_peak = False 
    times.append(0.) 
    prims.append(y_init) 
    t = times[0] 
    y = prims[0] 
    tmax = 2.0e9 
    while t < tmax: 
        status, t, h, y = e.apply(c,s,sys,t,tmax,h,y) 
        if status != errno.GSL_SUCCESS: 
                print("break status", status) 
                break 
        prims.append(y) 
        times.append(t) 
     
        r = y[0] 
        pphi = y[3] 
        if r <= 6: 
            rhs = v5HM_unoptimized_rhs(t,y,[m1,m2,chi1,chi2]) 
            drdt = rhs[0] 
            dphidt = rhs[1] 
            dprstardt = rhs[2] 
            if dphidt < omega_previous: 
                    omega_peak = True 
                    break 
            if drdt > 0: 
                    break 
            if dprstardt > 0: 
                prstar_peak = True 
                break 
            if r < rstop: 
                break 
            if r < 3: 
                phidot_circ = v5HM_unoptimized_omega_circ(m1,m2,r,pphi,chi1,chi2) 
                if phidot_circ > 1: 
                    break 
            omega_previous = dphidt 
     
    t = np.array(times) 
    dynamics = np.array(prims) 
    step_back_time = 250. 
    if omega_peak: 
        t_desired = t[-1] - step_back_time - 50 
    else: 
        t_desired = t[-1] - step_back_time 
     
    coarse_fine_separation_idx = np.argmin(np.abs(t - t_desired)) 
     
    dynamics_coarse = np.c_[t[:coarse_fine_separation_idx],dynamics[:coarse_fine_separation_idx]] 
    dynamics_fine_prelim = np.c_[t[coarse_fine_separation_idx:],dynamics[coarse_fine_separation_idx:]] 
     
    dynamics_coarse = augment_dynamics(dynamics_coarse,m1,m2,chi1,chi2) 
    dynamics_fine_prelim = augment_dynamics(dynamics_fine_prelim,m1,m2,chi1,chi2) 
     
    t_peak = None 
     
    if omega_peak: 
        interpolant = CubicSpline(dynamics_fine_prelim[:,0],dynamics_fine_prelim[:,6]) 
        t_peak = iterative_refinement(interpolant.derivative(),[dynamics_fine_prelim[0,0],dynamics_fine[-1,0]]) 
    if pr_peak: 
        interpolant = CubicSpline(dynamics_fine_prelim[:,0],dynamics_fine_prelim[:,3]) 
        t_peak = iterative_refinement(interpolant.derivative(),[dynamics_fine_prelim[-1,0] - 10, dynamics_fine[-1,0]]) 
     
    dynamics_fine = interpolate_dynamics(dynamics_fine_prelim[:,:5],omega_peak,step_back_time) 
    dynamics_fine = augment_dynamics(dynamics_fine,m1,m2,chi1,chi2) 
    return dynamics_coarse, dynamics_fine

import numpy as np
from scipy.signal import argrelmin
from scipy.interpolate import CubicSpline
from Hamiltonian.v5HM_unoptimized_hamiltonian import v5HM_unoptimized_hamiltonian as H
from Derivatives.v5HM_Hamiltonian_Derivatives_unoptimized import v5HM_unoptimized_dH_dpphi as omega
from Derivatives.v5HM_Hamiltonian_Derivatives_unoptimized import v5HM_unoptimized_omega_circ as omega_circ


def augment_dynamics(dynamics, m1, m2, chi1, chi2):
    result = []
    eta = m1*m2/(m1 + m2)/(m1 + m2)
    N = dynamics.shape[0]
    for i in range(N):
        dyn = dynamics[i]
        Hreal,xi = H(m1,m2,dyn[1],dyn[3],dyn[4],chi1,chi2)
        Omega = omega(m1,m2,dyn[1],dyn[3],dyn[4],chi1,chi2)/eta
        Omega_c = omega_circ(m1,m2,dyn[1],dyn[4],chi1,chi2)/eta
        result.append([Hreal,Omega,Omega_c])
    result = np.array(result)
    return np.c_[dynamics,result]

def iterative_refinement(var_dot, interval, peak_pr = False):
    left = interval[0]
    right = interval[1]
    dt = 0.1
    for n in range(2):
        dt /= 10
        t_array_fine = np.arange(interval[0],interval[1],dt)
        abs_var_dot = np.abs(var_dot(t_array_fine))
        minima = argrelmin(abs_var_dot, order = 3)[0]
        
        if len(minima) > 0:
            result = t_array_fine[minima[0]]
            interval = [ max(left,result - 10*dt) , min(right,result + 10*dt)]
        else:
            if peak_pr:
                return interval[-1]
            else:
                return ( interval[0] + interval[-1] ) / 2
    return result    

def interpolate_dynamics(dynamics_fine, omega_peak = None, step_back = 250):
    res = []
    n = len(dynamics_fine)
    
    if omega_peak:
        t_start_fine = max(omega_peak - step_back, dynamics_fine[0,0])
        t_new = np.arange(t_start, peak_omega, 0.1)
    else:
        t_new = np.arange(dynamics_fine[0,0], dynamics_fine[-1,0], 0.1)
    
    for i in range(1, dynamics_fine.shape[1]):
        interpolant = CubicSpline(dynamics_fine[:,0], dynamics_fine[:,1])
        res.append(interpolant(t_new))
    
    res = np.array(res)
    res_transposed = res.T
    
    return np.c_[t_new,res_transposed]

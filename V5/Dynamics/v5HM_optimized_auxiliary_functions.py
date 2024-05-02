
import numpy as np
from scipy.signal import argrelmin
from scipy.interpolate import CubicSpline
from pygsl_lite import spline
from Derivatives.v5HM_BOB_optimized_hamiltonian_first_derivatives import v5HM_BOB_optimized_hamiltonian_first_derivatives as dH
from Derivatives.v5HM_BOB_optimized_omega import v5HM_BOB_optimized_omega as omega_circ
from Radiation.v5HM_BOB_optimized_waveform import v5HM_BOB_optimized_waveform as hlm

def augment_dynamics(dynamics, m1, m2, chi1, chi2,a6,dSO):
    result = []
    N = dynamics.shape[0]
    for i in range(N):
        dyn = dynamics[i]
        dvalues, ham_and_xi = dH(m1,m2,dyn[1],dyn[3],dyn[4],chi1,chi2,a6,dSO)
        Omega_c = omega_circ(m1,m2,dyn[1],dyn[4],chi1,chi2,a6,dSO)
        result.append([ham_and_xi[0],dvalues[3],Omega_c])
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
        t_new = np.arange(t_start_fine, omega_peak, 0.1)
    else:
        t_new = np.arange(dynamics_fine[0,0], dynamics_fine[-1,0], 0.1)
    
    for i in range(1, dynamics_fine.shape[1]):
        interpolant = CubicSpline(dynamics_fine[:,0], dynamics_fine[:,i])
        res.append(interpolant(t_new))
    
    res = np.array(res)
    res_transposed = res.T
    
    return np.c_[t_new,res_transposed]

def get_waveforms_inspiral(m1,m2,augmented_dynamics,chi1,chi2):
    nu = m1*m2/(m1 + m2)/(m1 + m2)
    N = augmented_dynamics.shape[0]
    h22 = []
    for i in range(N):
        r = augmented_dynamics[i,1]
        phi = augmented_dynamics[i,2]
        prstar = augmented_dynamics[i,3]
        pphi = augmented_dynamics[i,4]
        Hreal = augmented_dynamics[i,5]
        Omega = augmented_dynamics[i,6]
        Omega_circ = augmented_dynamics[i,7]
        h22.append(hlm(m1,m2,r,phi,prstar,pphi,chi1,chi2,Hreal,Omega,Omega_circ))
    h22 = np.array(h22)
    return h22

def interpolate_modes_fast(t_new,waveform_mode, dynamics):
    n = len(waveform_mode)
    orbital_phase = dynamics[:,2]
    intrp_orbital_phase = spline.cspline(n)
    intrp_orbital_phase.init(dynamics[:,0],dynamics[:,2])
    orbital_phase_intrp = intrp_orbital_phase.eval_e_vector(t_new)
    intrp_real_mode, intrp_imag_mode = spline.cspline(n), spline.cspline(n)
    h22_remove_phase_scaling = waveform_mode*np.exp(1j*2*orbital_phase)
    h22_nophase_real = np.real(h22_remove_phase_scaling)
    h22_nophase_imag = np.imag(h22_remove_phase_scaling)
    intrp_real_mode.init(dynamics[:,0],h22_nophase_real)
    intrp_imag_mode.init(dynamics[:,0],h22_nophase_imag)
    h22_real = intrp_real_mode.eval_e_vector(t_new)
    h22_imag = intrp_imag_mode.eval_e_vector(t_new)
    h22_rescaled = (h22_real + 1j*h22_imag)*np.exp(-1j*2*orbital_phase_intrp)
    return h22_rescaled    

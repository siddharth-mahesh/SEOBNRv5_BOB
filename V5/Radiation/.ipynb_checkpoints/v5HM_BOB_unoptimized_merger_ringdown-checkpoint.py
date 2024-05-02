import numpy as np
def v5HM_BOB_unoptimized_merger_ringdown(t,t0,hNR,omegaNR,omegaQNM,tau):
    Omega0 = omegaNR/2
    OmegaQNM = omegaQNM/2
    tp = t0 - 2*tau*np.log(Omega0/OmegaQNM)
    Ap = hNR*(omegaNR**2)*np.cosh((t0 - tp)/tau)
    k = (OmegaQNM**4 - Omega0**4)/(1 - np.tanh((t0 - tp)/tau))
    kappam = Omega0*Omega0/OmegaQNM
    kappap = OmegaQNM
    Omega = ( Omega0**4 + k*( np.tanh((t - tp)/tau) - np.tanh((t0 - tp)/tau) ) )**(1/4)
    Omega0_over_kappap = Omega0/kappap
    Omega0_over_kappam = Omega0/kappam
    Omega_over_kappap = Omega/kappap
    Omega_over_kappam = Omega/kappam
    arctanhm = 0.5*kappam*tau*np.log( (1 + Omega_over_kappam)*(1 - Omega0_over_kappam)/((1 - Omega_over_kappam)*(1 + Omega0_over_kappam)) )
    arctanhp = 0.5*kappap*tau*np.log( (1 + Omega_over_kappap)*(1 - Omega0_over_kappap)/((1 - Omega_over_kappap)*(1 + Omega0_over_kappap)) )
    arctanm = kappam*tau*( np.arctan2(Omega,kappam) - np.arctan2(Omega0,kappam) )
    arctanp = kappap*tau*( np.arctan2(Omega,kappap) - np.arctan2(Omega0,kappap) )
    h = (Ap/4/(Omega**2))*(1/np.cosh((t - tp)/tau))
    Phi = arctanp + arctanhp - arctanm - arctanhm
    phi = 2*Phi
    return h,phi
import numpy as np
def v5HM_BOB_unoptimized_merger_ringdown(t,t0,hNR,omegaNR,omegaQNM,tau):
    Omega0 = omegaNR/2
    OmegaQNM = omegaQNM/2
    tp = t0 - 2*tau*np.log(Omega0/OmegaQNM)
    Ap = hNR*(omegaNR**2)*np.cosh((t0 - tp)/tau)
    k = (OmegaQNM**4 - Omega0**4)/(1 - np.tanh((t0 - tp)/tau))
    kappam = (Omega0**4 - k*( 1 + np.tanh((t0 - tp)/tau) ))**(1/4)
    kappap = (Omega0**4 + k*( 1 - np.tanh((t0 - tp)/tau) ))**(1/4)
    print(Omega0, kappap, kappam)
    Omega = ( Omega0**4 + k*( np.tanh((t - tp)/tau) - np.tanh((t0 - tp)/tau) ) )**(1/4)
    atanh_Omega0_kappam = np.divide(1,2)*np.log( (1 + Omega0/kappam)/(1 - Omega0/kappam) )
    atanh_Omega_kappam = np.divide(1,2)*np.log( (1 + Omega/kappam)/(1 - Omega/kappam) )
    atanh_Omega0_kappap = np.divide(1,2)*np.log( (1 + Omega0/kappap)/(1 - Omega0/kappap) )
    atanh_Omega_kappap = np.divide(1,2)*np.log( (1 + Omega/kappap)/(1 - Omega/kappap) )
    arctanhm = kappam*tau*( atanh_Omega_kappam - atanh_Omega0_kappam )
    arctanhp = kappap*tau*( atanh_Omega_kappap - atanh_Omega0_kappap )
    arctanm = kappam*tau*( np.arctan2(Omega,kappam) - np.arctan2(Omega0,kappam) )
    arctanp = kappap*tau*( np.arctan2(Omega,kappap) - np.arctan2(Omega0,kappap) )
    h = (Ap/4/(Omega**2))*(1/np.cosh((t - tp)/tau))
    Phi = arctanp + arctanhp - arctanm - arctanhm
    phi = 2*Phi
    return h,np.unwrap(phi)

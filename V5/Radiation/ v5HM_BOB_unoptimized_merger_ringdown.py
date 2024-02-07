import numpy as np
def v5HM_BOB_unoptimized_merger_ringdown(t,t0,hNR,omegaNR,OmegaQNM,tau,NQC = False):
    Omega0 = omegaNR/2
    tp = t0 - 2*tau*np.log(Omega0/OmegaQNM)
    Ap = hNR*(omegaNR**2)*sp.cosh((t0 - tp)/tau)
    k = (OmegaQNM**4 - Omega0**4)/(1 - sp.tanh((t0 - tp)/tau))
    kappap = (Omega0**4 + k*( 1 - sp.tanh((t0 - tp)/tau) ))**(1/4)
    Omega = ( Omega0**4 + k*( sp.tanh((t - tp)/tau) - sp.tanh((t0 - tp)/tau) ) )**(1/4)
    atanh_Omega0_kappam = np.divide(1,2)*np.log( (1 + Omega0/kappam)/(1 - Omega0/kappam) )
    atanh_Omega_kappam = np.divide(1,2)*np.log( (1 + Omega/kappam)/(1 - Omega/kappam) )
    atanh_Omega0_kappap = np.divide(1,2)*np.log( (1 + Omega0/kappap)/(1 - Omega0/kappap) )
    atanh_Omega_kappap = np.divide(1,2)*np.log( (1 + Omega/kappap)/(1 - Omega/kappap) )
    arctanhm = kappam*tau*( atanh_Omega_kappam - atanh_Omega0_kappam )
    arctanhp = kappap*tau*( atanh_Omega_kappap - atanh_Omega0_kappap )
    arctanm = kappam*tau*( sp.atan(Omega/kappam) - sp.atan(Omega0/kappam) )
    arctanp = kappap*tau*( sp.atan(Omega/kappap) - sp.atan(Omega0/kappap) )
    Phi = arctanp + arctanhp - arctanm - arctanhm
    phi = 2*Phi
    h = (Ap/4/Omega2)*sp.sech((t - tp)/tau)
    hcross = h*sp.sin(phi)
    hplus = h*sp.cos(phi)
    if not NQC:
        return hplus,hcross
    else:
        return h,phi,Omega
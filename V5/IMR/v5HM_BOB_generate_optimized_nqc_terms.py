import sympy as sp
from nrpy.c_codegen import c_codegen
def v5HM_BOB_generate_optimized_nqc_terms():
    t,t0,hNR,omegaNR,omegaQNM,tau = sp.symbols('t t0 hNR omegaNR omegaQNM tau',real = True)
    Omega0 = omegaNR/2
    OmegaQNM = omegaQNM/2
    tp = t0 - 2*tau*sp.log(Omega0/OmegaQNM)
    Ap = hNR*(omegaNR**2)*sp.cosh((t0 - tp)/tau)
    k = (OmegaQNM**4 - Omega0**4)/(1 - sp.tanh((t0 - tp)/tau))
    kappam = Omega0*Omega0/OmegaQNM
    kappap = OmegaQNM
    Omega = ( Omega0**4 + k*( sp.tanh((t - tp)/tau) - sp.tanh((t0 - tp)/tau) ) )**(1/4)
    Omega0_over_kappap = Omega0/kappap
    Omega0_over_kappam = Omega0/kappam
    Omega_over_kappap = Omega/kappap
    Omega_over_kappam = Omega/kappam
    arctanhm = 0.5*kappam*tau*sp.log( (1 + Omega_over_kappam)*(1 - Omega0_over_kappam)/((1 - Omega_over_kappam)*(1 + Omega0_over_kappam)) )
    arctanhp = 0.5*kappap*tau*sp.log( (1 + Omega_over_kappap)*(1 - Omega0_over_kappap)/((1 - Omega_over_kappap)*(1 + Omega0_over_kappap)) )
    arctanm = kappam*tau*( sp.atan2(Omega,kappam) - sp.atan2(Omega0,kappam) )
    arctanp = kappap*tau*( sp.atan2(Omega,kappap) - sp.atan2(Omega0,kappap) )
    h = (Ap/4/(Omega**2))*(1/sp.cosh((t - tp)/tau))
    Phi = arctanp + arctanhp - arctanm - arctanhm
    phi = 2*Phi
    hdot = sp.diff(h,t)
    hddot = sp.diff(hdot,t)
    omega = 2*Omega
    omegadot = sp.diff(omega,t)
    return c_codegen([h,hdot,hddot,omega,omegadot],['damp[0]','damp[1]','damp[2]','domega[0]','domega[1]'],include_braces = False,verbose = False)
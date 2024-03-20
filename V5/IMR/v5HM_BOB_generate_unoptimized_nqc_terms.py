import sympy as sp
def v5HM_BOB_generate_unoptimized_nqc_terms():
    t,t0,hNR,omegaNR,omegaQNM,tau = sp.symbols('t t0 hNR omegaNR omegaQNM tau',real = True)
    Omega0 = omegaNR/2
    OmegaQNM = omegaQNM/2
    tp = t0 - 2*tau*sp.log(Omega0/OmegaQNM)
    Ap = hNR*(omegaNR**2)*sp.cosh((t0 - tp)/tau)
    k = (OmegaQNM**4 - Omega0**4)/(1 - sp.tanh((t0 - tp)/tau))
    kappam = (Omega0**4 - k*( 1 + sp.tanh((t0 - tp)/tau) ))**(1/4)
    kappap = (Omega0**4 + k*( 1 - sp.tanh((t0 - tp)/tau) ))**(1/4)
    Omega = ( Omega0**4 + k*( sp.tanh((t - tp)/tau) - sp.tanh((t0 - tp)/tau) ) )**(1/4)
    atanh_Omega0_kappam = sp.Rational(1,2)*sp.log( (1 + Omega0/kappam)/(1 - Omega0/kappam) )
    atanh_Omega_kappam = sp.Rational(1,2)*sp.log( (1 + Omega/kappam)/(1 - Omega/kappam) )
    atanh_Omega0_kappap = sp.Rational(1,2)*sp.log( (1 + Omega0/kappap)/(1 - Omega0/kappap) )
    atanh_Omega_kappap = sp.Rational(1,2)*sp.log( (1 + Omega/kappap)/(1 - Omega/kappap) )
    arctanhm = kappam*tau*( atanh_Omega_kappam - atanh_Omega0_kappam )
    arctanhp = kappap*tau*( atanh_Omega_kappap - atanh_Omega0_kappap )
    arctanm = kappam*tau*( sp.atan2(Omega,kappam) - sp.atan2(Omega0,kappam) )
    arctanp = kappap*tau*( sp.atan2(Omega,kappap) - sp.atan2(Omega0,kappap) )
    h = (Ap/4/(Omega**2))*(1/sp.cosh((t - tp)/tau))
    Phi = arctanp + arctanhp - arctanm - arctanhm
    phi = 2*Phi
    hdot = sp.diff(h,t)
    hddot = sp.diff(hdot,t)
    omega = 2*Omega
    omegadot = sp.diff(omega,t)
    return sp.pycode(h),sp.pycode(hdot),sp.pycode(hddot),sp.pycode(omega),sp.pycode(omegadot)
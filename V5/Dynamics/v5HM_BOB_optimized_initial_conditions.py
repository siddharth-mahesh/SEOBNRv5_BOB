import numpy as np
from scipy.optimize import root, root_scalar
from Derivatives.v5HM_BOB_optimized_conservative_initial_conditions import v5HM_BOB_optimized_conservative_initial_conditions as IC_cons
from Derivatives.v5HM_BOB_optimized_dissipative_initial_conditions import v5HM_BOB_optimized_dissipative_initial_conditions as IC_diss
def v5HM_BOB_optimized_initial_conditions(M,q,S1,S2,f,a6,dSO,Deltat):
    m1 = q/(1 + q) 
    m2 = 1/(1 + q) 
    nu = m1*m2 
    chi1 = S1 
    chi2 = S2 
    Msol = 4.925491025543575903411922162094833998e-6 
    Omega_0 = M*Msol*np.pi*f 
     
    h = 2*np.pi/5/Omega_0 
    r_guess = Omega_0**(-2/3) 
    pphi_guess = Omega_0**(-1/3) 
    sol_cons_guess = np.array([r_guess,pphi_guess]) 
    params_cons = np.array([m1,m2,chi1,chi2,a6,dSO,Omega_0]) 
    sol_cons = root(IC_cons,sol_cons_guess,args = params_cons,tol = 6e-12) 
    r , pphi = sol_cons.x[0] , sol_cons.x[1] 
    params_diss = np.array([m1,m2,r,pphi,chi1,chi2,a6,dSO]) 
    prstar_bracket = [-3e-2,0] 
    prstar_sol = root_scalar(IC_diss, args = params_diss, bracket = prstar_bracket, xtol = 1e-12, rtol = 1e-10) 
    prstar = prstar_sol.root 
    atot = m1*m1*chi1 + m2*m2*chi2 
    aeff = atot + .474046*nu*(chi1 + chi2) 
    aeff_max1 = min(aeff,1) 
    Z1eff = 1 + (np.cbrt(1 - aeff_max1*aeff_max1))*(np.cbrt(1 + aeff_max1) + np.cbrt(1 - aeff_max1)) 
    Z2eff = np.sqrt(3*aeff_max1*aeff_max1 + Z1eff*Z1eff) 
    rISCOeff = 3 + Z2eff - np.sign(aeff_max1)*np.sqrt((3 - Z1eff)*(3 + Z1eff + Z2eff)) 
    LISCOeff = (2/3/np.sqrt(3))*(1 + 2*np.sqrt(3*rISCOeff - 2)) 
    EISCOeff = np.sqrt(1 - 2/(3*rISCOeff)) 
    k = np.zeros([4,5]) 
    k[0,0] = -5.97723 
    k[0,1] = 3.39221 
    k[0,2] = 4.48865 
    k[0,3] = -5.77101 
    k[0,4] = -13.0459 
    k[1,0] = 35.1278 
    k[1,1] = -72.9336 
    k[1,2] = -86.0036 
    k[1,3] = 93.7371 
    k[1,4] = 200.975 
    k[2,0] = -146.822 
    k[2,1] = 387.184 
    k[2,2] = 447.009 
    k[2,3] = -467.383 
    k[2,4] = -884.339 
    k[3,0] = 223.911 
    k[3,1] = -648.502 
    k[3,2] = -697.177 
    k[3,3] = 753.738 
    k[3,4] = 1166.89 
     
    NRfactor = 0 
    nu_i = 1 
    for i in range(len(k)): 
        aeff_j = 1 
        for j in range(len(k[i])): 
            NRfactor += k[i,j]*nu_i*aeff_j 
            aeff_j *= aeff 
        nu_i *= nu 
    ell = np.abs(LISCOeff - 2*atot*(EISCOeff - 1) + nu*NRfactor) 
    af = atot + nu*ell 
    Z1f = 1 + (np.cbrt(1 - af*af))*(np.cbrt(1 + af) + np.cbrt(1 - af)) 
    Z2f = np.sqrt(3*af*af + Z1f*Z1f) 
    rISCO = 3 + Z2f - np.sign(af)*np.sqrt((3 - Z1f)*(3 + Z1f + Z2f)) 
    rstop = -1 
    if Deltat > 0: 
        rstop = 0.98*rISCO 
    Shat = (m1*m1*chi1 + m2*m2*chi2)/(m1*m1 + m2*m2) 
    Shat2 = Shat*Shat 
    Shat3 = Shat2*Shat 
    Shat4 = Shat3*Shat 
    nu2 = nu*nu 
    nu3 = nu2*nu 
    nu4 = nu3*nu 
    sqrt1m4nu = np.sqrt(1. - 4.*nu) 
    Deltachi = chi1-chi2 
    Deltachi2 = Deltachi*Deltachi 
    a2 = 0.5609904135313374 
    a3 = -0.84667563764404 
    a4 = 3.145145224278187 
    Erad_nu_0 = a4*nu4 + a3*nu3 + a2*nu2 + (1. - 2.*np.sqrt(2.)/3.)*nu 
    b1 = -0.2091189048177395 
    b2 = -0.19709136361080587 
    b3 = -0.1588185739358418 
    b5 = 2.9852925538232014 
    f20 = 4.271313308472851 
    f30 = 31.08987570280556 
    f50 = 1.5673498395263061 
    f10 = 1.8083565298668276 
    f21 = 0. 
    d10 = -0.09803730445895877 
    d11 = -3.2283713377939134 
    d20 = 0.01118530335431078 
    d30 = -0.01978238971523653 
    d31 = -4.91667749015812 
    f11 = 15.738082204419655 
    f31 = -243.6299258830685 
    f51 = -0.5808669012986468 
    bfin1 = b1*(f10 + f11*nu + (16. - 16.*f10 - 4.*f11)*nu2) 
    bfin2 = b2*(f20 + f21*nu + (16. - 16.*f20 - 4.*f21)*nu2) 
    bfin3 = b3*(f30 + f31*nu + (16. - 16.*f30 - 4.*f31)*nu2) 
    bfin5 = b5*(f50 + f51*nu + (16. - 16.*f50 - 4.*f51)*nu2) 
    Erad_eq_Shat = (0.128*bfin3*Shat3 + 0.211*bfin2*Shat2 + 0.346*bfin1*Shat + 1.)/(1 - 0.212*bfin5*Shat) 
    Erad_nu_Shat = Erad_nu_0*Erad_eq_Shat 
    d10 = -0.09803730445895877 
    d11 = -3.2283713377939134 
    d20 = 0.01118530335431078 
    d30 = -0.01978238971523653 
    d31 = -4.91667749015812 
    A_1 = d10*sqrt1m4nu*nu2*(d11*nu+1) 
    A_2 = d20*nu3 
    A_3 = d30*sqrt1m4nu*nu*(d31*nu + 1) 
    DeltaErad_nu_Shat_Deltachi = A_1*Deltachi + A_2*Deltachi2 + A_3*Deltachi*Shat 
    Mf = 1 - (Erad_nu_Shat + DeltaErad_nu_Shat_Deltachi) 
    h22NR = nu*np.abs(71.97969776036882194603 * nu4 
                - 13.35761402231352157344 * nu3 * Shat 
                - 46.87585958426210908101 * nu3 
                + 0.61988944517825661507 * nu2 * Shat2 
                + 7.19426416189229733789 * nu2 * Shat 
                + 12.44040490932310127903 * nu2 
                + 0.43014673069078152023 * nu * Shat3 
                - 1.74313546783413597652 * nu * Shat 
                - 0.86828935763242798274 * nu 
                - 0.08493901280736430859 * Shat3 
                - 0.02082621429567295401 * Shat2 
                + 0.18693991146784910695 * Shat 
                + 1.46709663479911811557) 
    omega22NR = (5.89352329617707670906 * nu4 
                + 3.75145580491965446868 * nu3 * Shat 
                - 3.34930536209472551334 * nu3 
                - 0.97140932625194231775 * nu2 * Shat2 
                - 1.69734302394369973577 * nu2 * Shat 
                + 0.28539204856044564362 * nu2 
                + 0.2419483723662931296 * nu * Shat3 
                + 0.51801427018052081941 * nu * Shat2 
                + 0.25096450064948544467 * nu * Shat 
                - 0.31709602351033533418 * nu 
                - 0.01525897668158244028 * Shat4 
                - 0.06692658483513916345 * Shat3 
                - 0.08715176045684569495 * Shat2 
                - 0.09133931944098934441 * Shat 
                - 0.2685414392185025978 
    ) 
    return m1,m2,chi1,chi2,np.array([r,0.,prstar,pphi]),Omega_0,h,rstop,rISCO,af,Mf,h22NR,omega22NR
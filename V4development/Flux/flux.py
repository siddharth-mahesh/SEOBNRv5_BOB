import numpy as np
import Flux.factorized_waveform_for_flux as wv

## Function to compute factorized flux
# Where does omega come from more generally? 
# Spin aligned EOB version is a parameter we may not need as all precessing 
# approximants build on v2 spin-aligned background

def compute_flux(m1,m2,EMGamma,tortoise,q,p,S1,S2,omega,Hreal,NQC=1.,lMax = 8):
    r = np.linalg.norm(q)
    eta = m1*m2/(m1+m2)/(m1+m2)
    flux = 0.
    omegasq = omega*omega
    v = np.cbrt(omega)
    rcrossp = np.cross(q,p)
    rcrosspmag = np.linalg.norm(rcrossp)
    m1hat = m1/(m1+m2)
    m2hat = m2/(m1+m2)
    S1_over_m12 = S1/m1hat/m1hat
    S2_over_m22 = S2/m2hat/m2hat
    Lhat = rcrossp/rcrosspmag
    s1dotL = np.dot(S1_over_m12,Lhat)
    s2dotL = np.dot(S2_over_m22,Lhat)
    
    #print("s1dotL = ",s1dotL)
    #print("s2dotL = ",s2dotL)
    
    chiS = 0.5*(s1dotL + s2dotL)
    chiA = 0.5*(s1dotL - s2dotL)
    #print("inputs:m1,m2,EMGamma,tortoise,q,p,S1,S2,Hreal,v,chiA,chiS")
    #print(m1,m2,EMGamma,tortoise,q,p,S1,S2,Hreal,v,chiA,chiS)
    for l in range(2,lMax+1):
        for m in range(1,l+1):
            #print(l,",",m)
            #print(m1,m2,EMGamma,tortoise,q,p,S1,S2,l,m,Hreal,v,chiA,chiS)
            hLM = wv.compute_hFlm_for_flux(m1,m2,EMGamma,tortoise,q,p,S1,S2,l,m,Hreal,v,chiA,chiS)
            ## NQC terms here
            if l == 2 and m ==2:
                hLM *= NQC
            flux += m*m*omegasq*hLM*hLM
            #print(l,m,hLM)
            #print("fluxlm = %.16e"%flux)
    if flux > 5 or omegasq > 1:
        flux = 0.
        
    return flux/np.pi/8.


    
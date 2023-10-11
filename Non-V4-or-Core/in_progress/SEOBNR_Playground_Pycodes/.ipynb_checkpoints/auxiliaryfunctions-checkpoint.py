import numpy as np
from SEOBNR_Playground_Pycodes.new_dHdp1 import new_compute_dHdp1
from SEOBNR_Playground_Pycodes.new_dHdp2 import new_compute_dHdp2
from SEOBNR_Playground_Pycodes.new_dHdp3 import new_compute_dHdp3
import SEOBNR_Playground_Pycodes.factorized_modes as fm

## Function to compute a double factorial; see https://en.wikipedia.org/wiki/Double_factorial
def doublefactorial(n):
     if n <= 0:
        return 1
     else:
        return n * doublefactorial(n-2)

## Compute the Newtonian prefix m
def Newtonian_n(m1,m2,l,m):
    epsilon = (l+m)%2
    doubfact = doublefactorial(2*l+1)
    n = np.power(complex(0,m), l)*8*np.divide(np.pi,doubfact)
    if epsilon==0:
        n *= np.sqrt(np.divide((l+1)*(l+2),l*(l-1)))
    elif epsilon==1:
        num = (2*l + 1)*(l + 2)*(l*l - m*m)
        div = (2*l - 1)*(l + 1)*(l)*(l - 1)
        n *= -2.j*np.sqrt( np.divide( num , div ) )
    else:
        print("Epsilon must be 0 or 1")
        exit()
    #print("n = %.16e"%np.abs(n))
    return np.abs(n)

## Compute the Newtonian prefix c
#     Eq. 7 of Damour, Iyer and Nagar 2008.
#     For odd m, c is proportional to dM = m1-m2. In the equal-mass case, c = dM = 0.
#     In the equal-mass unequal-spin case, however, when spins are different, the odd m term is generally not zero.
#     In this case, c can be written as c0 * dM, while spins terms in PN expansion may take the form chiA/dM.
#     Although the dM's cancel analytically, we can not implement c and chiA/dM with the possibility of dM -> 0.
#     Therefore, for this case, we give numerical values of c0 for relevant modes, and c0 is calculated as
#     c / dM in the limit of dM -> 0. Consistently, for this case, we implement chiA instead of chiA/dM
#     in LALSimIMRSpinEOBFactorizedWaveform.c.
    
def Newtonian_c(m1,m2,l,m):
    Mtot = m1 + m2
    m1hat = np.divide(m1,Mtot)
    m2hat = np.divide(m2,Mtot)
    epsilon = (l+m)%2
    if (m%2)==0:
        sign = 1
    else:
        sign = -1
    lpepm1 = l + epsilon - 1
    if (m1!=m2) or sign==1:
        c = np.power(m2hat,lpepm1) + sign*np.power(m1hat,lpepm1)
    else:
        if l==2 or l==3:
            c = -1.
        elif l==4 or l==5:
            c = -0.5
        else:
            c = 0.
    #print("c = %.16e"%c)
    return c

## Compute the Associate Legendre Polynomial for input value 0
def AssociatedLegendre(l,m):
    if l==1:
        if m==1:
            return -1.
        else:
            print("You used a bad (l,m)")
    if l==2:
        if m==2:
            return 3.
        elif m==1:
            return 0.
        else:
            print("You used a bad (l,m)")
    if l==3:
        if m==3:
            return -15.
        elif m==2:
            return 0.
        elif m==1:
            return 1.5
        else:
            print("You used a bad (l,m)")
    if l==4:
        if m==4:
            return 105.
        elif m==3:
            return 0.
        elif m==2:
            return -7.5
        elif m==1:
            return 0.
        else:
            print("You used a bad (l,m)")
    if l==5:
        if m==5:
            return -945.
        elif m==4:
            return 0.
        elif m==3:
            return 52.5
        elif m==2:
            return 0.
        elif m==1:
            return -1.875
        else:
            print("You used a bad (l,m)")
    if l==6:
        if m==6:
            return 10395.
        elif m==5:
            return 0.
        elif m==4:
            return -472.5
        elif m==3:
            return 0.
        elif m==2:
            return 13.125
        elif m==1:
            return 0.
        else:
            print("You used a bad (l,m)")
    if l==7:
        if m==7:
            return -135135.
        elif m==6:
            return 0.
        elif m==5:
            return 5197.5
        elif m==4:
            return 0.
        elif m==3:
            return -118.125
        elif m==2:
            return 0.
        elif m==1:
            return 2.1875
        else:
            print("You used a bad (l,m)")
    if l==8:
        if m==8:
            return 2027025.
        elif m==7:
            return 0.
        elif m==6:
            return -67567.5
        elif m==5:
            return 0.
        elif m==4:
            return 1299.375
        elif m==3:
            return 0.
        elif m==2:
            return -19.6875
        elif m==1:
            return 0.
        else:
            print("You used a bad (l,m)")

## Compute the absolute value of Scalar Spherical Harmonics at the equatorial plane
def AbsSphericalHarmonicAtPiOver2(l,m):
    absM = abs( m )
    legendre = AssociatedLegendre(l,absM)
    #print("legendre = %.16e"%legendre)
    #Since we assume negative m as input, we compute the prefactor accordingly
    result = legendre* np.sqrt( (2*l + 1)*np.math.factorial(l - absM) / ( 4*np.pi*np.math.factorial(l+absM) ) )
    if (m < 0 and absM % 2 == 1):
        result *= -1
    #print("ylm = %.16e"%result)
    return result

## Compute the Source Term Se_eff
def Se_eff(l,m,Hreal,v,q,p,eta):
    epsilon = (l + m)%2
    if epsilon == 0:
        return (Hreal*Hreal - 1)/(2*eta) + 1
    rcrossp = np.cross(q,p)
    return v*np.linalg.norm(rcrossp)

## Compute the product factor in T_lm
def Tlmprodfac(l,hathatk):
    result = 1
    for s in range(1,l+1):
        result *= (s*s + 4*hathatk*hathatk)
    #print(result)
    return result

## Compute the exact Circular Frequency at the given phase space points
# Question: why is this necessary if we are already passing omega through the flux function?
# if we are not computing waveform then what is being passed as v?

def CalcOmega(m1,m2,EMGamma,tortoise,q,p,S1,S2):
    #print("q = ", q)
    #print("p = ", p)
    #print("S1 = ", S1)
    #print("S2 = ", S2)
    eta = m1*m2/(m1 + m2)/(m1 + m2)
    xdot = new_compute_dHdp1(m1,m2,EMGamma,tortoise,q[0],q[1],q[2],p[0],p[1],p[2],S1[0],S1[1],S1[2],S2[0],S2[1],S2[2])[0]
    ydot = new_compute_dHdp2(m1,m2,EMGamma,tortoise,q[0],q[1],q[2],p[0],p[1],p[2],S1[0],S1[1],S1[2],S2[0],S2[1],S2[2])[0]
    zdot = new_compute_dHdp3(m1,m2,EMGamma,tortoise,q[0],q[1],q[2],p[0],p[1],p[2],S1[0],S1[1],S1[2],S2[0],S2[1],S2[2])[0]
    qdot = np.array([xdot,ydot,zdot])/eta
    #print("qdot = ", qdot)
    
    Lnhat = np.cross(q,qdot)
    Lnhat /= np.linalg.norm(Lnhat)
    #print("Lnhat = ", Lnhat)
    xhat = np.array([1.,0.,0.])
    yhat = np.array([0.,1.,0.])
    zhat = np.array([0.,0.,1.])
    
    if np.dot(Lnhat,xhat) < 0.9:
        R1 = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1]])
        xprimehat = Lnhat
    else:
        invsqrt2 = 1/np.sqrt(2)
        R1 = np.array([[invsqrt2,-invsqrt2,0.],[invsqrt2,invsqrt2,0.],[0.,0.,1.]])
        xprimehat = np.matmul(R1,Lnhat)
    
    yprimehat = np.cross(xprimehat,xhat)
    yprimehat /= np.linalg.norm(yprimehat)
    
    zprimehat = np.cross(xprimehat,yprimehat)
    zprimehat /= np.linalg.norm(zprimehat)
    
    R2 = np.array([[xprimehat[0],xprimehat[1],xprimehat[2]],[yprimehat[0],yprimehat[1],yprimehat[2]],[zprimehat[0],zprimehat[1],zprimehat[2]]])
    #print("R2 = ", R2)
    
    qtemp = np.matmul(R1,q)
    ptemp = np.matmul(R1,p)
    S1temp = np.matmul(R1,S1)
    S2temp = np.matmul(R1,S2)
    Lnhattemp = np.matmul(R1,Lnhat)
    
    qprime = np.matmul(R2,qtemp)
    pprime = np.matmul(R2,ptemp)
    S1prime = np.matmul(R2,S1temp)
    S2prime = np.matmul(R2,S2temp)
    LNhatprime = np.matmul(R2,Lnhattemp)
    #print("qprime = " , np.matmul(R2,qtemp))
    #print("pprime = " , np.matmul(R2,ptemp))
    #print("S1prime = " , np.matmul(R2,S1temp))
    #print("S2prime = " , np.matmul(R2,S2temp))
    #print("LNhatprime = " , np.matmul(R2,Lnhattemp))
    
    ## No need to compute spherical momenta as we have analytical derivatives
    
    rpass = np.linalg.norm(qprime)
    thpass = np.arccos(qprime[0]/rpass)
    phipass = np.arctan2(-qprime[1],qprime[2])
    #print("rpass = " , np.linalg.norm(qprime))
    #print("thpass = " , np.arccos(qprime[0]/rpass))
    #print("phipass = " , np.arctan2(-qprime[1],qprime[2]))
    
    dHdpz = new_compute_dHdp3(m1,m2,EMGamma,tortoise,qprime[0],qprime[1],qprime[2],pprime[0],pprime[1],pprime[2],S1prime[0],S1prime[1],S1prime[2],S2prime[0],S2prime[1],S2prime[2])
    dHdpy = new_compute_dHdp2(m1,m2,EMGamma,tortoise,qprime[0],qprime[1],qprime[2],pprime[0],pprime[1],pprime[2],S1prime[0],S1prime[1],S1prime[2],S2prime[0],S2prime[1],S2prime[2])
    
    omega = -(np.cos(phipass)*dHdpy + np.sin(phipass)*dHdpz)/rpass/np.sin(thpass)/eta
    #print("omega = ", omega)
    
    return omega
    
## Compute the non-Keplerian coefficient for vPhi
def vPhiNonKeplerian(m1,m2,EMGamma,tortoise,q,p,S1,S2):
    omega_circular = CalcOmega(m1,m2,EMGamma,tortoise,q,p,S1,S2)
    r = np.linalg.norm(q)
    r3 = r*r*r
    #print("vPhiKepler = %.16e" % (1/(omega_circular*omega_circular*r3)) )
    return 1/(omega_circular*omega_circular*r3)
    
## Function to compute the Post-Newtonian Waveform modes

def rholmpowl(m1,m2,l,m,chiA,chiS,v,EMgamma):
    eta = m1*m2/(m1+m2)/(m1+m2)
    tplspin = (1. - 2.*eta)*chiS + chiA*(m1-m2)/(m1+m2)
    #print("validate_a = %.16e"%tplspin)
    #print(m1,m2,tplspin,eta,chiA,chiS)
    test = fm.compute_modes(m1,m2,tplspin,eta,chiA,chiS)
    vsq = v*v
    eulerlog = EMgamma + np.log(2.*m*v)
    auxflm = 0
    if l==2:
        if m==2:
            rholm = 1 + vsq*(test['rho22v2'] + v*(test['rho22v3'] + v*(test['rho22v4'] + v*(test['rho22v5'] + v*(test['rho22v6']
                             + test['rho22v6l']*eulerlog + v*(test['rho22v7'] + v*(test['rho22v8'] + test['rho22v8l']*eulerlog
                             + (test['rho22v10'] + test['rho22v10l']*eulerlog)*vsq)))))))
        elif m==1:
            rholm = 1. + v*(test['rho21v1'] + v*(test['rho21v2'] + v*(test['rho21v3'] + v*(test['rho21v4'] + v*(test['rho21v5']
                            + v*(test['rho21v6'] + test['rho21v6l']*eulerlog + v*(test['rho21v7'] + test['rho21v7l']*eulerlog
                            + v*(test['rho21v8'] + test['rho21v8l']*eulerlog + (test['rho21v10'] + test['rho21v10l']*eulerlog)*vsq))))))))
            auxflm = v*test['f21v1'] + vsq*v*test['f21v3']
        else:
            print("You used a bad (l,m)")
    elif l==3:
        if m==3:
            rholm = 1. + vsq*(test['rho33v2'] + v*(test['rho33v3'] + v*(test['rho33v4'] + v*(test['rho33v5'] + v*(test['rho33v6']
                            + test['rho33v6l']*eulerlog + v*(test['rho33v7'] + (test['rho33v8'] + test['rho33v8l']*eulerlog)*v))))))
            auxflm = v*vsq*test['f33v3'];
        elif m==2:
            rholm = 1. + v*(test['rho32v1'] + v*(test['rho32v2'] + v*(test['rho32v3'] + v*(test['rho32v4'] + v*(test['rho32v5']
                            + v*(test['rho32v6'] + test['rho32v6l']*eulerlog + (test['rho32v8'] + test['rho32v8l']*eulerlog)*vsq))))))
        elif m==1:
            rholm = 1. + vsq*(test['rho31v2'] + v*(test['rho31v3'] + v*(test['rho31v4'] + v*(test['rho31v5'] + v*(test['rho31v6']
                            + test['rho31v6l']*eulerlog + v*(test['rho31v7'] + (test['rho31v8'] + test['rho31v8l']*eulerlog)*v))))))
            auxflm = v*vsq*test['f31v3']
        else:
            print("You used a bad (l,m)")
    elif l==4:
        if m==4:
            rholm = 1. + vsq*(test['rho44v2'] + v*(test['rho44v3'] + v*(test['rho44v4'] + v*(test['rho44v5'] + (test['rho44v6']
                            + test['rho44v6l']*eulerlog)*v))))
        elif m==3:
            rholm = 1. + v*(test['rho43v1'] + v*(test['rho43v2'] + vsq*(test['rho43v4'] + v*(test['rho43v5'] + (test['rho43v6']
                            + test['rho43v6l']*eulerlog)*v))))
            auxflm = v*test['f43v1']
        elif m==2:
            rholm = 1. + vsq*(test['rho42v2'] + v*(test['rho42v3'] + v*(test['rho42v4'] + v*(test['rho42v5'] + (test['rho42v6']
                            + test['rho42v6l']*eulerlog)*v))))
        elif m==1:
            rholm = 1. + v*(test['rho41v1'] + v*(test['rho41v2'] + vsq*(test['rho41v4'] + v*(test['rho41v5'] + (test['rho41v6']
                            + test['rho41v6l']*eulerlog)*v))))
            auxflm = v*test['f41v1']
        else:
            print("You used a bad (l,m)")
    elif l==5:
        if m==5:
            rholm = 1. + vsq*(test['rho55v2'] + v*(test['rho55v3'] + v*(test['rho55v4'] + v*(test['rho55v5'] + test['rho55v6']*v))))
        elif m==4:
            rholm = 1. + vsq*(test['rho54v2'] + v*(test['rho54v3'] + test['rho54v4']*v))
        elif m==3:
            rholm = 1. + vsq*(test['rho53v2'] + v*(test['rho53v3'] + v*(test['rho53v4'] + test['rho53v5']*v)))
        elif m==2:
            rholm = 1. + vsq*(test['rho52v2'] + v*(test['rho52v3'] + test['rho52v4']*v))
        elif m==1:
            rholm = 1. + vsq*(test['rho51v2'] + v*(test['rho51v3'] + v*(test['rho51v4'] + test['rho51v5']*v)))
        else:
            print("You used a bad (l,m)")
    elif l==6:
        if m==6:
            rholm = 1. + vsq*(test['rho66v2'] + v*(test['rho66v3'] + test['rho66v4']*v))
        elif m==5:
            rholm = 1. + vsq*(test['rho65v2'] + test['rho65v3']*v)
        elif m==4:
            rholm = 1. + vsq*(test['rho64v2'] + v*(test['rho64v3'] + test['rho64v4']*v))
        elif m==3:
            rholm = 1. + vsq*(test['rho63v2'] + test['rho63v3']*v)
        elif m==2:
            rholm = 1. + vsq*(test['rho62v2'] + v*(test['rho62v3'] + test['rho62v4']*v))
        elif m==1:
            rholm = 1. + vsq*(test['rho61v2'] + test['rho61v3']*v)
        else:
            print("You used a bad (l,m)")
    elif l==7:
        if m==7:
            rholm = 1. + vsq*(test['rho77v2'] + test['rho77v3']*v)
        elif m==6:
            rholm = 1. + test['rho76v2']*vsq
        elif m==5:
            rholm = 1. + vsq*(test['rho75v2'] + test['rho75v3']*v)
        elif m==4:
            rholm = 1. + test['rho74v2']*vsq
        elif m==3:
            rholm = 1. + vsq*(test['rho73v2'] + test['rho73v3']*v)
        elif m==2:
            rholm = 1. + test['rho72v2']*vsq
        elif m==1:
            rholm = 1. + vsq*(test['rho71v2'] + test['rho71v3']*v)
        else:
            print("You used a bad (l,m)")
    elif l==8:
        if m==8:
            rholm = 1. + test['rho88v2']*vsq
        elif m==7:
            rholm = 1. + test['rho87v2']*vsq
        elif m==6:
            rholm = 1. + test['rho86v2']*vsq
        elif m==5:
            rholm = 1. + test['rho85v2']*vsq
        elif m==4:
            rholm = 1. + test['rho84v2']*vsq
        elif m==3:
            rholm = 1. + test['rho83v2']*vsq
        elif m==2:
            rholm = 1. + test['rho82v2']*vsq
        elif m==1:
            rholm = 1. + test['rho81v2']*vsq
        else:
            print("You used a bad (l,m)")
    else:
        print("You used a bad (l,m)")
    rholmPowl = np.power(rholm,l)
    # Since in the clm computations, we had made adjustments to the
    # equal-mass unequal spin cases, below we make the corresponding adjustment
    # to rholmPowl in order to ensure that we get the correct answer in the limt dM -> 0
    if eta==0.25 and (m % 2):
        rholmPowl = auxflm
    else:
        rholmPowl += auxflm
    #print(rholmPowl)
    return rholmPowl

    
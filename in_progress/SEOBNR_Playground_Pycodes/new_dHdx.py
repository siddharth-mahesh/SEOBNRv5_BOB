from __future__ import division
import numpy as np
def new_compute_dHdx(m1, m2, eta, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z, KK, k0, k1, dSO, dSS, tortoise, EMgamma):
    M = m1+m2
    mu = m1*m2/M
    eta = mu/M
    r = np.sqrt(x*x+y*y+z*z)
    u = 1/r
    sigmastar3 = m2/m1*S1z+m1/m2*S2z
    sigmastar2 = m2/m1*S1y+m1/m2*S2y
    sigmastar1 = m2/m1*S1x+m1/m2*S2x
    sigma3 = S1z+S2z
    sigma2 = S1y+S2y
    sigma1 = S1x+S2x
    Skerr3 = sigma3
    Skerr2 = sigma2
    Skerr1 = sigma1
    Skerrmag = np.sqrt(Skerr1*Skerr1+Skerr2*Skerr2+Skerr3*Skerr3)
    Skerrhat3 = Skerr3/Skerrmag
    Skerrhat2 = Skerr2/Skerrmag
    Skerrhat1 = Skerr1/Skerrmag
    a = Skerrmag
    n3 = z/r
    n2 = y/r
    n1 = x/r
    e33 = Skerrhat3
    e32 = Skerrhat2
    e31 = Skerrhat1
    xi3 = e31*n2-e32*n1
    xi2 = -e31*n3+e33*n1
    xi1 = e32*n3-e33*n2
    v3 = n1*xi2-n2*xi1
    v2 = n3*xi1-n1*xi3
    v1 = n2*xi3-n3*xi2
    costheta = e31*n1+e32*n2+e33*n3
    sin2theta = 1-costheta*costheta
    xisq = sin2theta
    w2 = a*a+r*r
    Sigma = r*r+a*a*costheta*costheta
    Dinv = 1+np.log(1+6*eta*u*u+2*(26-3*eta)*eta*u*u*u)
    omegatilde = 2*a*r
    K = 1.712-1.803949138004582*eta-39.77229225266885*eta*eta+103.16588921239249*eta*eta*eta
    etaKminus1 = eta*K-1
    Delta0 = K*(eta*K-2)
    Delta1 = -2*etaKminus1*(K+Delta0)
    Delta2 = np.true_divide(1,2)*Delta1*(Delta1-4*etaKminus1)-a*a*etaKminus1*etaKminus1*Delta0
    Delta3 = -np.true_divide(1,3)*Delta1*Delta1*Delta1+etaKminus1*Delta1*Delta1+Delta2*Delta1-2*etaKminus1*(Delta2-etaKminus1)-a*a*etaKminus1*etaKminus1*Delta1
    Delta4 = np.true_divide(1,12)*(6*a*a*(Delta1*Delta1-2*Delta2)*etaKminus1*etaKminus1+3*Delta1*Delta1*Delta1*Delta1-8*etaKminus1*Delta1*Delta1*Delta1-12*Delta2*Delta1*Delta1+12*(2*etaKminus1*Delta2+Delta3)*Delta1+12*(np.true_divide(94,3)-np.true_divide(41,32)*np.pi*np.pi)*etaKminus1*etaKminus1+6*(Delta2*Delta2-4*Delta3*etaKminus1))
    Delta5 = etaKminus1*etaKminus1*((np.true_divide(-4237,60)+np.true_divide(128,5)*EMgamma+np.true_divide(2275,512)*np.pi*np.pi-np.true_divide(1,3)*a*a*(Delta1*Delta1*Delta1-3*Delta1*Delta2+3*Delta3)-(Delta1*Delta1*Delta1*Delta1*Delta1-5*Delta1*Delta1*Delta1*Delta2+5*Delta1*Delta2*Delta2+5*Delta1*Delta1*Delta3-5*Delta2*Delta3-5*Delta1*Delta4)/(5*etaKminus1*etaKminus1)+(Delta1*Delta1*Delta1*Delta1-4*Delta1*Delta1*Delta2+2*Delta2*Delta2+4*Delta1*Delta3-4*Delta4)/(2*etaKminus1)+np.true_divide(256,5)*np.log(2)))
    Delta5l = etaKminus1*etaKminus1*np.true_divide(64,5)
    logarg = u*(Delta1+u*(Delta2+u*(Delta3+u*(Delta4+u*(Delta5+Delta5l*np.log(u))))))
    Deltaucalib = 1+eta*(Delta0+np.log(1+logarg))
    Deltaucalibprime = -eta*u*u*(Delta1+u*(2*Delta2+u*(3*Delta3+u*(4*Delta4+u*(5*(Delta5+Delta5l*np.log(u)))))))/(1+logarg)
    Deltaubar = a*a*u*u+(2*u+1/etaKminus1)/etaKminus1
    Deltaubarprime = -2*a*a*u*u*u-2*u*u/(etaKminus1)
    Deltau = Deltaubar*Deltaucalib
    Deltauprime = Deltaubarprime*Deltaucalib+Deltaubar*Deltaucalibprime
    Deltatprime = 2*r*Deltau+r*r*Deltauprime
    Deltat = r*r*Deltau
    Deltar = Deltat*Dinv
    Lambt = w2*w2-a*a*Deltat*sin2theta
    csi = np.sqrt(Deltar*Deltat)/w2
    csi1 = 1+(1-abs(1-tortoise))*(csi-1)
    csi2 = 1+(np.true_divide(1,2)-np.true_divide(1,2)*np.sign(np.true_divide(3,2)-tortoise))*(csi-1)
    prT = csi2*(p1*n1+p2*n2+p3*n3)
    phat3 = p3+prT*(1-1/csi1)*n3
    phat2 = p2+prT*(1-1/csi1)*n2
    phat1 = p1+prT*(1-1/csi1)*n1
    pdotxir = (phat1*xi1+phat2*xi2+phat3*xi3)*r
    pdotn = phat1*n1+phat2*n2+phat3*n3
    pdotvr = (phat1*v1+phat2*v2+phat3*v3)*r
    pphi = pdotxir
    Qcoeff2 = 1/(Sigma*sin2theta)
    Qcoeff1 = Sigma/(Lambt*sin2theta)
    DrSipn2 = Deltar*pdotn*pdotn/Sigma
    Q = 1+DrSipn2+Qcoeff1*pdotxir*pdotxir+Qcoeff2*pdotvr*pdotvr
    Qminus1 = Q-1
    Jtilde = np.sqrt(Deltar)
    exp2mu = Sigma
    expmu = np.sqrt(exp2mu)
    Brtilde = (np.sqrt(Deltar)*Deltatprime-2*Deltat)/(2*np.sqrt(Deltar*Deltat))
    Btilde = np.sqrt(Deltat)
    exp2nu = Deltat*Sigma/Lambt
    expnu = np.sqrt(exp2nu)
    omega = omegatilde/Lambt
    omegatildeprime = 2*a
    Lambtprime = 4*(a*a+r*r)*r-2*a*a*Deltatprime*sin2theta
    mucostheta = a*a*costheta/Sigma
    nucostheta = a*a*w2*costheta*(w2-Deltat)/(Lambt*Sigma)
    omegacostheta = -2*a*a*costheta*Deltat*omegatilde/(Lambt*Lambt)
    mur = r/Sigma-1/np.sqrt(Deltar)
    nur = r/Sigma+w2*(w2*Deltatprime-4*r*Deltat)/(2*Lambt*Deltat)
    omegar = (Lambt*omegatildeprime-Lambtprime*omegatilde)/(Lambt*Lambt)
    dSO = -74.71-156.*eta+627.5*eta*eta
    sigmacoeffTerm3 = eta*dSO*u*u*u
    sigmacoeffTerm2 = eta/(144*r*r)*(-896+r*(-436*Qminus1-96*DrSipn2+r*(-45*Qminus1*Qminus1+36*Qminus1*DrSipn2))+eta*(-336+r*(204*Qminus1-882*DrSipn2+r*(810*DrSipn2*DrSipn2-234*Qminus1*DrSipn2))))
    sigmacoeffTerm1 = eta/12*(-8/r+3*Qminus1-36*DrSipn2)
    sigmacoeff = sigmacoeffTerm1+sigmacoeffTerm2+sigmacoeffTerm3
    sigmastarcoeffTerm2 = eta/(72*r*r)*(706+r*(-206*Qminus1+282*DrSipn2+r*Qminus1*(96*DrSipn2-23*Qminus1))+eta*(-54+r*(120*Qminus1-324*DrSipn2+r*(360*DrSipn2*DrSipn2+Qminus1*(-126*DrSipn2-3*Qminus1)))))
    sigmastarcoeffTerm1 = eta/12*(14/r+4*Qminus1-30*DrSipn2)
    sigmastarcoeff = sigmastarcoeffTerm1+sigmastarcoeffTerm2
    Deltasigmastar3 = sigmastar3*sigmastarcoeff+sigma3*sigmacoeff
    Deltasigmastar2 = sigmastar2*sigmastarcoeff+sigma2*sigmacoeff
    Deltasigmastar1 = sigmastar1*sigmastarcoeff+sigma1*sigmacoeff
    Sstar3 = sigmastar3+Deltasigmastar3
    Sstar2 = sigmastar2+Deltasigmastar2
    Sstar1 = sigmastar1+Deltasigmastar1
    S3 = Sstar3
    S2 = Sstar2
    S1 = Sstar1
    Sstardotn = Sstar1*n1+Sstar2*n2+Sstar3*n3
    SdotSkerrhat = S1*Skerrhat1+S2*Skerrhat2+S3*Skerrhat3
    Sdotn = S1*n1+S2*n2+S3*n3
    Sdotv = S1*v1+S2*v2+S3*v3
    Sdotxi = S1*xi1+S2*xi2+S3*xi3
    HdsumTerm2 = 3*Sstardotn*Sstardotn
    HdsumTerm1 = Sstar1*Sstar1+Sstar2*Sstar2+Sstar3*Sstar3
    Hdsum = HdsumTerm1-HdsumTerm2
    Hdcoeff = np.true_divide(1,2)/(r*r*r)
    Q4 = 2*prT*prT*prT*prT*u*u*(4-3*eta)*eta
    gammappsum = Deltar/Sigma*pdotn*pdotn+1/Sigma*pdotvr*pdotvr/sin2theta+Sigma/Lambt/sin2theta*pdotxir*pdotxir
    Hnsradicand = 1+gammappsum+Q4
    alpha = np.sqrt(Deltat*Sigma/Lambt)
    betapsum = omegatilde*pphi/Lambt
    HssTerm3 = expmu*expnu*pdotxir*(Jtilde*pdotn*Sdotxi*Btilde-expmu*expnu*pdotxir*Sdotn)+(pdotvr*(Sdotn*pdotvr-Jtilde*pdotn*Sdotv)-exp2mu*(np.sqrt(Q)+Q)*Sdotn*xisq)*Btilde*Btilde
    HssTerm3coeff = omegacostheta/(2*exp2mu*expmu*expnu*Btilde*(Q+np.sqrt(Q)))
    HssTerm2 = expmu*pdotxir*(expmu*exp2nu*pdotxir*Sdotv-expnu*pdotvr*Sdotxi*Btilde)+xisq*Btilde*Btilde*(exp2mu*(np.sqrt(Q)+Q)*Sdotv+Jtilde*pdotn*(pdotvr*Sdotn-Jtilde*pdotn*Sdotv))
    HssTerm2coeff = Jtilde*omegar/(2*exp2mu*expmu*expnu*Btilde*(Q+np.sqrt(Q))*xisq)
    HssTerm1 = omega*SdotSkerrhat
    Hss = HssTerm1+HssTerm2coeff*HssTerm2+HssTerm3coeff*HssTerm3
    HsoTerm2c = Jtilde*Brtilde*expmu*expnu*pdotxir*(np.sqrt(Q)+1)*Sdotv
    HsoTerm2b = expmu*expnu*pdotxir*(2*np.sqrt(Q)+1)*(Jtilde*nur*Sdotv-nucostheta*Sdotn*xisq)*Btilde
    HsoTerm2a = Sdotxi*Jtilde*(mur*pdotvr*(np.sqrt(Q)+1)-mucostheta*pdotn*xisq-np.sqrt(Q)*(nur*pdotvr+(mucostheta-nucostheta)*pdotn*xisq))*Btilde*Btilde
    HsoTerm2 = HsoTerm2a+HsoTerm2b-HsoTerm2c
    HsoTerm2coeff = expnu/(exp2mu*Btilde*Btilde*(Q+np.sqrt(Q))*xisq)
    HsoTerm1 = exp2nu*(expmu*expnu-Btilde)*pdotxir*SdotSkerrhat/(expmu*Btilde*Btilde*np.sqrt(Q)*xisq)
    Hso = HsoTerm1+HsoTerm2coeff*HsoTerm2
    Hd = Hdcoeff*Hdsum
    Hns = betapsum+alpha*np.sqrt(Hnsradicand)
    Hs = Hso+Hss
    dSS = 8.127-154.2*eta+830.8*eta*eta
    Heff = Hs+Hns-Hd+dSS*eta*u*u*u*u*(S1x*S1x+S1y*S1y+S1z*S1z+S2x*S2x+S2y*S2y+S2z*S2z)
    rprm_x = x/np.sqrt(x**2 + y**2 + z**2)
    uprm_x = -rprm_x/r**2
    n3prm_x = -rprm_x*z/r**2
    n2prm_x = -rprm_x*y/r**2
    n1prm_x = 1/r - rprm_x*x/r**2
    xi3prm_x = e31*n2prm_x - e32*n1prm_x
    xi2prm_x = -e31*n3prm_x + e33*n1prm_x
    xi1prm_x = e32*n3prm_x - e33*n2prm_x
    v3prm_x = n1*xi2prm_x + n1prm_x*xi2 - n2*xi1prm_x - n2prm_x*xi1
    v2prm_x = -n1*xi3prm_x - n1prm_x*xi3 + n3*xi1prm_x + n3prm_x*xi1
    v1prm_x = n2*xi3prm_x + n2prm_x*xi3 - n3*xi2prm_x - n3prm_x*xi2
    costhetaprm_x = e31*n1prm_x + e32*n2prm_x + e33*n3prm_x
    sin2thetaprm_x = -2*costheta*costhetaprm_x
    xisqprm_x = sin2thetaprm_x
    w2prm_x = 2*r*rprm_x
    Sigmaprm_x = 2*a**2*costheta*costhetaprm_x + 2*r*rprm_x
    Dinvprm_x = (3*eta*u**2*uprm_x*(-6*eta + 52) + 12*eta*u*uprm_x)/(eta*u**3*(-6*eta + 52) + 6*eta*u**2 + 1)
    omegatildeprm_x = 2*a*rprm_x
    logargprm_x = u*(u*(u*(u*(Delta5l*uprm_x + uprm_x*(Delta5 + Delta5l*np.log(u))) + uprm_x*(Delta4 + u*(Delta5 + Delta5l*np.log(u)))) + uprm_x*(Delta3 + u*(Delta4 + u*(Delta5 + Delta5l*np.log(u))))) + uprm_x*(Delta2 + u*(Delta3 + u*(Delta4 + u*(Delta5 + Delta5l*np.log(u)))))) + uprm_x*(Delta1 + u*(Delta2 + u*(Delta3 + u*(Delta4 + u*(Delta5 + Delta5l*np.log(u))))))
    Deltaucalibprm_x = eta*logargprm_x/(logarg + 1)
    Deltaucalibprimeprm_x = eta*logargprm_x*u**2*(Delta1 + u*(2*Delta2 + u*(3*Delta3 + u*(4*Delta4 + u*(5*Delta5 + 5*Delta5l*np.log(u))))))/(logarg + 1)**2 - eta*u**2*(u*(u*(u*(5*Delta5l*uprm_x + uprm_x*(5*Delta5 + 5*Delta5l*np.log(u))) + uprm_x*(4*Delta4 + u*(5*Delta5 + 5*Delta5l*np.log(u)))) + uprm_x*(3*Delta3 + u*(4*Delta4 + u*(5*Delta5 + 5*Delta5l*np.log(u))))) + uprm_x*(2*Delta2 + u*(3*Delta3 + u*(4*Delta4 + u*(5*Delta5 + 5*Delta5l*np.log(u))))))/(logarg + 1) - 2*eta*u*uprm_x*(Delta1 + u*(2*Delta2 + u*(3*Delta3 + u*(4*Delta4 + u*(5*Delta5 + 5*Delta5l*np.log(u))))))/(logarg + 1)
    Deltaubarprm_x = 2*a**2*u*uprm_x + 2*uprm_x/etaKminus1
    Deltaubarprimeprm_x = -6*a**2*u**2*uprm_x - 4*u*uprm_x/etaKminus1
    Deltauprm_x = Deltaubar*Deltaucalibprm_x + Deltaubarprm_x*Deltaucalib
    Deltauprimeprm_x = Deltaubar*Deltaucalibprimeprm_x + Deltaubarprime*Deltaucalibprm_x + Deltaubarprimeprm_x*Deltaucalib + Deltaubarprm_x*Deltaucalibprime
    Deltatprimeprm_x = 2*Deltau*rprm_x + 2*Deltauprime*r*rprm_x + Deltauprimeprm_x*r**2 + 2*Deltauprm_x*r
    Deltatprm_x = 2*Deltau*r*rprm_x + Deltauprm_x*r**2
    Deltarprm_x = Deltat*Dinvprm_x + Deltatprm_x*Dinv
    Lambtprm_x = -Deltat*a**2*sin2thetaprm_x - Deltatprm_x*a**2*sin2theta + 2*w2*w2prm_x
    csiprm_x = -w2prm_x*np.sqrt(Deltar*Deltat)/w2**2 + np.sqrt(Deltar*Deltat)*(Deltar*Deltatprm_x/2 + Deltarprm_x*Deltat/2)/(Deltar*Deltat*w2)
    csi1prm_x = csiprm_x*(-abs(-tortoise + 1) + 1)
    csi2prm_x = csiprm_x*(-np.sign(-tortoise + 3/2)/2 + 1/2)
    prTprm_x = csi2*(n1prm_x*p1 + n2prm_x*p2 + n3prm_x*p3) + csi2prm_x*(n1*p1 + n2*p2 + n3*p3)
    phat3prm_x = n3*prTprm_x*(1 - 1/csi1) + n3prm_x*prT*(1 - 1/csi1) + csi1prm_x*n3*prT/csi1**2
    phat2prm_x = n2*prTprm_x*(1 - 1/csi1) + n2prm_x*prT*(1 - 1/csi1) + csi1prm_x*n2*prT/csi1**2
    phat1prm_x = n1*prTprm_x*(1 - 1/csi1) + n1prm_x*prT*(1 - 1/csi1) + csi1prm_x*n1*prT/csi1**2
    pdotxirprm_x = r*(phat1*xi1prm_x + phat1prm_x*xi1 + phat2*xi2prm_x + phat2prm_x*xi2 + phat3*xi3prm_x + phat3prm_x*xi3) + rprm_x*(phat1*xi1 + phat2*xi2 + phat3*xi3)
    pdotnprm_x = n1*phat1prm_x + n1prm_x*phat1 + n2*phat2prm_x + n2prm_x*phat2 + n3*phat3prm_x + n3prm_x*phat3
    pdotvrprm_x = r*(phat1*v1prm_x + phat1prm_x*v1 + phat2*v2prm_x + phat2prm_x*v2 + phat3*v3prm_x + phat3prm_x*v3) + rprm_x*(phat1*v1 + phat2*v2 + phat3*v3)
    pphiprm_x = pdotxirprm_x
    Qcoeff2prm_x = -sin2thetaprm_x/(Sigma*sin2theta**2) - Sigmaprm_x/(Sigma**2*sin2theta)
    Qcoeff1prm_x = -Sigma*sin2thetaprm_x/(Lambt*sin2theta**2) + Sigmaprm_x/(Lambt*sin2theta) - Lambtprm_x*Sigma/(Lambt**2*sin2theta)
    DrSipn2prm_x = 2*Deltar*pdotn*pdotnprm_x/Sigma - Deltar*Sigmaprm_x*pdotn**2/Sigma**2 + Deltarprm_x*pdotn**2/Sigma
    Qprm_x = DrSipn2prm_x + 2*Qcoeff1*pdotxir*pdotxirprm_x + Qcoeff1prm_x*pdotxir**2 + 2*Qcoeff2*pdotvr*pdotvrprm_x + Qcoeff2prm_x*pdotvr**2
    Qminus1prm_x = Qprm_x
    Jtildeprm_x = Deltarprm_x/(2*np.sqrt(Deltar))
    exp2muprm_x = Sigmaprm_x
    expmuprm_x = exp2muprm_x/(2*np.sqrt(exp2mu))
    Brtildeprm_x = (np.sqrt(Deltar)*Deltatprimeprm_x - 2*Deltatprm_x + Deltarprm_x*Deltatprime/(2*np.sqrt(Deltar)))/(2*np.sqrt(Deltar*Deltat)) + (np.sqrt(Deltar)*Deltatprime - 2*Deltat)*(-Deltar*Deltatprm_x/2 - Deltarprm_x*Deltat/2)/(2*Deltar*Deltat*np.sqrt(Deltar*Deltat))
    Btildeprm_x = Deltatprm_x/(2*np.sqrt(Deltat))
    exp2nuprm_x = Deltat*Sigmaprm_x/Lambt - Deltat*Lambtprm_x*Sigma/Lambt**2 + Deltatprm_x*Sigma/Lambt
    expnuprm_x = exp2nuprm_x/(2*np.sqrt(exp2nu))
    omegaprm_x = omegatildeprm_x/Lambt - Lambtprm_x*omegatilde/Lambt**2
    Lambtprimeprm_x = -2*Deltatprime*a**2*sin2thetaprm_x - 2*Deltatprimeprm_x*a**2*sin2theta + 8*r**2*rprm_x + rprm_x*(4*a**2 + 4*r**2)
    mucosthetaprm_x = a**2*costhetaprm_x/Sigma - Sigmaprm_x*a**2*costheta/Sigma**2
    nucosthetaprm_x = a**2*costheta*w2*(-Deltatprm_x + w2prm_x)/(Lambt*Sigma) + a**2*costheta*w2prm_x*(-Deltat + w2)/(Lambt*Sigma) + a**2*costhetaprm_x*w2*(-Deltat + w2)/(Lambt*Sigma) - Sigmaprm_x*a**2*costheta*w2*(-Deltat + w2)/(Lambt*Sigma**2) - Lambtprm_x*a**2*costheta*w2*(-Deltat + w2)/(Lambt**2*Sigma)
    omegacosthetaprm_x = -2*Deltat*a**2*costheta*omegatildeprm_x/Lambt**2 - 2*Deltat*a**2*costhetaprm_x*omegatilde/Lambt**2 + 4*Deltat*Lambtprm_x*a**2*costheta*omegatilde/Lambt**3 - 2*Deltatprm_x*a**2*costheta*omegatilde/Lambt**2
    murprm_x = rprm_x/Sigma - Sigmaprm_x*r/Sigma**2 + Deltarprm_x/(2*Deltar**(3/2))
    nurprm_x = rprm_x/Sigma - Sigmaprm_x*r/Sigma**2 + w2*(-4*Deltat*rprm_x + Deltatprime*w2prm_x + Deltatprimeprm_x*w2 - 4*Deltatprm_x*r)/(2*Deltat*Lambt) + w2prm_x*(-4*Deltat*r + Deltatprime*w2)/(2*Deltat*Lambt) - Lambtprm_x*w2*(-4*Deltat*r + Deltatprime*w2)/(2*Deltat*Lambt**2) - Deltatprm_x*w2*(-4*Deltat*r + Deltatprime*w2)/(2*Deltat**2*Lambt)
    omegarprm_x = (-Lambtprime*omegatildeprm_x - Lambtprimeprm_x*omegatilde + Lambtprm_x*omegatildeprime)/Lambt**2 + Lambtprm_x*(-2*Lambt*omegatildeprime + 2*Lambtprime*omegatilde)/Lambt**3
    sigmacoeffTerm3prm_x = 3*dSO*eta*u**2*uprm_x
    sigmacoeffTerm2prm_x = eta*(eta*(r*(-882*DrSipn2prm_x + 204*Qminus1prm_x + r*(1620*DrSipn2*DrSipn2prm_x - 234*DrSipn2*Qminus1prm_x - 234*DrSipn2prm_x*Qminus1) + rprm_x*(810*DrSipn2**2 - 234*DrSipn2*Qminus1)) + rprm_x*(-882*DrSipn2 + 204*Qminus1 + r*(810*DrSipn2**2 - 234*DrSipn2*Qminus1))) + r*(-96*DrSipn2prm_x - 436*Qminus1prm_x + r*(36*DrSipn2*Qminus1prm_x + 36*DrSipn2prm_x*Qminus1 - 90*Qminus1*Qminus1prm_x) + rprm_x*(36*DrSipn2*Qminus1 - 45*Qminus1**2)) + rprm_x*(-96*DrSipn2 - 436*Qminus1 + r*(36*DrSipn2*Qminus1 - 45*Qminus1**2)))/(144*r**2) - eta*rprm_x*(eta*(r*(-882*DrSipn2 + 204*Qminus1 + r*(810*DrSipn2**2 - 234*DrSipn2*Qminus1)) - 336) + r*(-96*DrSipn2 - 436*Qminus1 + r*(36*DrSipn2*Qminus1 - 45*Qminus1**2)) - 896)/(72*r**3)
    sigmacoeffTerm1prm_x = eta*(-36*DrSipn2prm_x + 3*Qminus1prm_x + 8*rprm_x/r**2)/12
    sigmacoeffprm_x = sigmacoeffTerm1prm_x + sigmacoeffTerm2prm_x + sigmacoeffTerm3prm_x
    sigmastarcoeffTerm2prm_x = eta*(eta*(r*(-324*DrSipn2prm_x + 120*Qminus1prm_x + r*(720*DrSipn2*DrSipn2prm_x + Qminus1*(-126*DrSipn2prm_x - 3*Qminus1prm_x) + Qminus1prm_x*(-126*DrSipn2 - 3*Qminus1)) + rprm_x*(360*DrSipn2**2 + Qminus1*(-126*DrSipn2 - 3*Qminus1))) + rprm_x*(-324*DrSipn2 + 120*Qminus1 + r*(360*DrSipn2**2 + Qminus1*(-126*DrSipn2 - 3*Qminus1)))) + r*(282*DrSipn2prm_x + Qminus1*r*(96*DrSipn2prm_x - 23*Qminus1prm_x) + Qminus1*rprm_x*(96*DrSipn2 - 23*Qminus1) + Qminus1prm_x*r*(96*DrSipn2 - 23*Qminus1) - 206*Qminus1prm_x) + rprm_x*(282*DrSipn2 + Qminus1*r*(96*DrSipn2 - 23*Qminus1) - 206*Qminus1))/(72*r**2) - eta*rprm_x*(eta*(r*(-324*DrSipn2 + 120*Qminus1 + r*(360*DrSipn2**2 + Qminus1*(-126*DrSipn2 - 3*Qminus1))) - 54) + r*(282*DrSipn2 + Qminus1*r*(96*DrSipn2 - 23*Qminus1) - 206*Qminus1) + 706)/(36*r**3)
    sigmastarcoeffTerm1prm_x = eta*(-30*DrSipn2prm_x + 4*Qminus1prm_x - 14*rprm_x/r**2)/12
    sigmastarcoeffprm_x = sigmastarcoeffTerm1prm_x + sigmastarcoeffTerm2prm_x
    Deltasigmastar3prm_x = sigma3*sigmacoeffprm_x + sigmastar3*sigmastarcoeffprm_x
    Deltasigmastar2prm_x = sigma2*sigmacoeffprm_x + sigmastar2*sigmastarcoeffprm_x
    Deltasigmastar1prm_x = sigma1*sigmacoeffprm_x + sigmastar1*sigmastarcoeffprm_x
    Sstar3prm_x = Deltasigmastar3prm_x
    Sstar2prm_x = Deltasigmastar2prm_x
    Sstar1prm_x = Deltasigmastar1prm_x
    S3prm_x = Sstar3prm_x
    S2prm_x = Sstar2prm_x
    S1prm_x = Sstar1prm_x
    Sstardotnprm_x = Sstar1*n1prm_x + Sstar1prm_x*n1 + Sstar2*n2prm_x + Sstar2prm_x*n2 + Sstar3*n3prm_x + Sstar3prm_x*n3
    SdotSkerrhatprm_x = S1prm_x*Skerrhat1 + S2prm_x*Skerrhat2 + S3prm_x*Skerrhat3
    Sdotnprm_x = S1*n1prm_x + S1prm_x*n1 + S2*n2prm_x + S2prm_x*n2 + S3*n3prm_x + S3prm_x*n3
    Sdotvprm_x = S1*v1prm_x + S1prm_x*v1 + S2*v2prm_x + S2prm_x*v2 + S3*v3prm_x + S3prm_x*v3
    Sdotxiprm_x = S1*xi1prm_x + S1prm_x*xi1 + S2*xi2prm_x + S2prm_x*xi2 + S3*xi3prm_x + S3prm_x*xi3
    HdsumTerm2prm_x = 6*Sstardotn*Sstardotnprm_x
    HdsumTerm1prm_x = 2*Sstar1*Sstar1prm_x + 2*Sstar2*Sstar2prm_x + 2*Sstar3*Sstar3prm_x
    Hdsumprm_x = HdsumTerm1prm_x - HdsumTerm2prm_x
    Hdcoeffprm_x = -3*rprm_x/(2*r**4)
    Q4prm_x = 4*eta*prT**4*u*uprm_x*(-3*eta + 4) + 8*eta*prT**3*prTprm_x*u**2*(-3*eta + 4)
    gammappsumprm_x = 2*Deltar*pdotn*pdotnprm_x/Sigma - Deltar*Sigmaprm_x*pdotn**2/Sigma**2 + Deltarprm_x*pdotn**2/Sigma - pdotvr**2*sin2thetaprm_x/(Sigma*sin2theta**2) + 2*pdotvr*pdotvrprm_x/(Sigma*sin2theta) - Sigmaprm_x*pdotvr**2/(Sigma**2*sin2theta) - Sigma*pdotxir**2*sin2thetaprm_x/(Lambt*sin2theta**2) + 2*Sigma*pdotxir*pdotxirprm_x/(Lambt*sin2theta) + Sigmaprm_x*pdotxir**2/(Lambt*sin2theta) - Lambtprm_x*Sigma*pdotxir**2/(Lambt**2*sin2theta)
    Hnsradicandprm_x = Q4prm_x + gammappsumprm_x
    alphaprm_x = Lambt*np.sqrt(Deltat*Sigma/Lambt)*(Deltat*Sigmaprm_x/(2*Lambt) - Deltat*Lambtprm_x*Sigma/(2*Lambt**2) + Deltatprm_x*Sigma/(2*Lambt))/(Deltat*Sigma)
    betapsumprm_x = omegatilde*pphiprm_x/Lambt + omegatildeprm_x*pphi/Lambt - Lambtprm_x*omegatilde*pphi/Lambt**2
    HssTerm3prm_x = Btilde**2*(Sdotn*exp2mu*xisq*(-Qprm_x - Qprm_x/(2*np.sqrt(Q))) + Sdotn*exp2mu*xisqprm_x*(-np.sqrt(Q) - Q) + Sdotn*exp2muprm_x*xisq*(-np.sqrt(Q) - Q) + Sdotnprm_x*exp2mu*xisq*(-np.sqrt(Q) - Q) + pdotvr*(-Jtilde*Sdotv*pdotnprm_x - Jtilde*Sdotvprm_x*pdotn - Jtildeprm_x*Sdotv*pdotn + Sdotn*pdotvrprm_x + Sdotnprm_x*pdotvr) + pdotvrprm_x*(-Jtilde*Sdotv*pdotn + Sdotn*pdotvr)) + Btilde*Btildeprm_x*(-2*Sdotn*exp2mu*xisq*(np.sqrt(Q) + Q) + 2*pdotvr*(-Jtilde*Sdotv*pdotn + Sdotn*pdotvr)) + expmu*expnu*pdotxir*(Btilde*Jtilde*Sdotxi*pdotnprm_x + Btilde*Jtilde*Sdotxiprm_x*pdotn + Btilde*Jtildeprm_x*Sdotxi*pdotn + Btildeprm_x*Jtilde*Sdotxi*pdotn - Sdotn*expmu*expnu*pdotxirprm_x - Sdotn*expmu*expnuprm_x*pdotxir - Sdotn*expmuprm_x*expnu*pdotxir - Sdotnprm_x*expmu*expnu*pdotxir) + expmu*expnu*pdotxirprm_x*(Btilde*Jtilde*Sdotxi*pdotn - Sdotn*expmu*expnu*pdotxir) + expmu*expnuprm_x*pdotxir*(Btilde*Jtilde*Sdotxi*pdotn - Sdotn*expmu*expnu*pdotxir) + expmuprm_x*expnu*pdotxir*(Btilde*Jtilde*Sdotxi*pdotn - Sdotn*expmu*expnu*pdotxir)
    HssTerm3coeffprm_x = omegacostheta*(-Qprm_x - Qprm_x/(2*np.sqrt(Q)))/(2*Btilde*exp2mu*expmu*expnu*(np.sqrt(Q) + Q)**2) + omegacosthetaprm_x/(Btilde*exp2mu*expmu*expnu*(2*np.sqrt(Q) + 2*Q)) - expnuprm_x*omegacostheta/(Btilde*exp2mu*expmu*expnu**2*(2*np.sqrt(Q) + 2*Q)) - expmuprm_x*omegacostheta/(Btilde*exp2mu*expmu**2*expnu*(2*np.sqrt(Q) + 2*Q)) - exp2muprm_x*omegacostheta/(Btilde*exp2mu**2*expmu*expnu*(2*np.sqrt(Q) + 2*Q)) - Btildeprm_x*omegacostheta/(Btilde**2*exp2mu*expmu*expnu*(2*np.sqrt(Q) + 2*Q))
    HssTerm2prm_x = Btilde**2*xisq*(Jtilde*pdotn*(-Jtilde*Sdotv*pdotnprm_x - Jtilde*Sdotvprm_x*pdotn - Jtildeprm_x*Sdotv*pdotn + Sdotn*pdotvrprm_x + Sdotnprm_x*pdotvr) + Jtilde*pdotnprm_x*(-Jtilde*Sdotv*pdotn + Sdotn*pdotvr) + Jtildeprm_x*pdotn*(-Jtilde*Sdotv*pdotn + Sdotn*pdotvr) + Sdotv*exp2mu*(Qprm_x + Qprm_x/(2*np.sqrt(Q))) + Sdotv*exp2muprm_x*(np.sqrt(Q) + Q) + Sdotvprm_x*exp2mu*(np.sqrt(Q) + Q)) + Btilde**2*xisqprm_x*(Jtilde*pdotn*(-Jtilde*Sdotv*pdotn + Sdotn*pdotvr) + Sdotv*exp2mu*(np.sqrt(Q) + Q)) + Btilde*Btildeprm_x*xisq*(2*Jtilde*pdotn*(-Jtilde*Sdotv*pdotn + Sdotn*pdotvr) + 2*Sdotv*exp2mu*(np.sqrt(Q) + Q)) + expmu*pdotxir*(-Btilde*Sdotxi*expnu*pdotvrprm_x - Btilde*Sdotxi*expnuprm_x*pdotvr - Btilde*Sdotxiprm_x*expnu*pdotvr - Btildeprm_x*Sdotxi*expnu*pdotvr + Sdotv*exp2nu*expmu*pdotxirprm_x + Sdotv*exp2nu*expmuprm_x*pdotxir + Sdotv*exp2nuprm_x*expmu*pdotxir + Sdotvprm_x*exp2nu*expmu*pdotxir) + expmu*pdotxirprm_x*(-Btilde*Sdotxi*expnu*pdotvr + Sdotv*exp2nu*expmu*pdotxir) + expmuprm_x*pdotxir*(-Btilde*Sdotxi*expnu*pdotvr + Sdotv*exp2nu*expmu*pdotxir)
    HssTerm2coeffprm_x = Jtilde*omegar*(-Qprm_x - Qprm_x/(2*np.sqrt(Q)))/(2*Btilde*exp2mu*expmu*expnu*xisq*(np.sqrt(Q) + Q)**2) - Jtilde*omegar*xisqprm_x/(Btilde*exp2mu*expmu*expnu*xisq**2*(2*np.sqrt(Q) + 2*Q)) + Jtilde*omegarprm_x/(Btilde*exp2mu*expmu*expnu*xisq*(2*np.sqrt(Q) + 2*Q)) - Jtilde*expnuprm_x*omegar/(Btilde*exp2mu*expmu*expnu**2*xisq*(2*np.sqrt(Q) + 2*Q)) - Jtilde*expmuprm_x*omegar/(Btilde*exp2mu*expmu**2*expnu*xisq*(2*np.sqrt(Q) + 2*Q)) - Jtilde*exp2muprm_x*omegar/(Btilde*exp2mu**2*expmu*expnu*xisq*(2*np.sqrt(Q) + 2*Q)) + Jtildeprm_x*omegar/(Btilde*exp2mu*expmu*expnu*xisq*(2*np.sqrt(Q) + 2*Q)) - Btildeprm_x*Jtilde*omegar/(Btilde**2*exp2mu*expmu*expnu*xisq*(2*np.sqrt(Q) + 2*Q))
    HssTerm1prm_x = SdotSkerrhat*omegaprm_x + SdotSkerrhatprm_x*omega
    Hssprm_x = HssTerm1prm_x + HssTerm2*HssTerm2coeffprm_x + HssTerm2coeff*HssTerm2prm_x + HssTerm3*HssTerm3coeffprm_x + HssTerm3coeff*HssTerm3prm_x
    HsoTerm2cprm_x = Brtilde*Jtilde*Sdotv*expmu*expnu*pdotxirprm_x*(np.sqrt(Q) + 1) + Brtilde*Jtilde*Sdotv*expmu*expnuprm_x*pdotxir*(np.sqrt(Q) + 1) + Brtilde*Jtilde*Sdotv*expmuprm_x*expnu*pdotxir*(np.sqrt(Q) + 1) + Brtilde*Jtilde*Sdotvprm_x*expmu*expnu*pdotxir*(np.sqrt(Q) + 1) + Brtilde*Jtilde*Qprm_x*Sdotv*expmu*expnu*pdotxir/(2*np.sqrt(Q)) + Brtilde*Jtildeprm_x*Sdotv*expmu*expnu*pdotxir*(np.sqrt(Q) + 1) + Brtildeprm_x*Jtilde*Sdotv*expmu*expnu*pdotxir*(np.sqrt(Q) + 1)
    HsoTerm2bprm_x = Btilde*expmu*expnu*pdotxir*(2*np.sqrt(Q) + 1)*(Jtilde*Sdotv*nurprm_x + Jtilde*Sdotvprm_x*nur + Jtildeprm_x*Sdotv*nur - Sdotn*nucostheta*xisqprm_x - Sdotn*nucosthetaprm_x*xisq - Sdotnprm_x*nucostheta*xisq) + Btilde*expmu*expnu*pdotxirprm_x*(2*np.sqrt(Q) + 1)*(Jtilde*Sdotv*nur - Sdotn*nucostheta*xisq) + Btilde*expmu*expnuprm_x*pdotxir*(2*np.sqrt(Q) + 1)*(Jtilde*Sdotv*nur - Sdotn*nucostheta*xisq) + Btilde*expmuprm_x*expnu*pdotxir*(2*np.sqrt(Q) + 1)*(Jtilde*Sdotv*nur - Sdotn*nucostheta*xisq) + Btilde*Qprm_x*expmu*expnu*pdotxir*(Jtilde*Sdotv*nur - Sdotn*nucostheta*xisq)/np.sqrt(Q) + Btildeprm_x*expmu*expnu*pdotxir*(2*np.sqrt(Q) + 1)*(Jtilde*Sdotv*nur - Sdotn*nucostheta*xisq)
    HsoTerm2aprm_x = Btilde**2*Jtilde*Sdotxi*(-np.sqrt(Q)*(nur*pdotvrprm_x + nurprm_x*pdotvr + pdotn*xisq*(mucosthetaprm_x - nucosthetaprm_x) + pdotn*xisqprm_x*(mucostheta - nucostheta) + pdotnprm_x*xisq*(mucostheta - nucostheta)) - mucostheta*pdotn*xisqprm_x - mucostheta*pdotnprm_x*xisq - mucosthetaprm_x*pdotn*xisq + mur*pdotvrprm_x*(np.sqrt(Q) + 1) + murprm_x*pdotvr*(np.sqrt(Q) + 1) + Qprm_x*mur*pdotvr/(2*np.sqrt(Q)) + Qprm_x*(-nur*pdotvr - pdotn*xisq*(mucostheta - nucostheta))/(2*np.sqrt(Q))) + Btilde**2*Jtilde*Sdotxiprm_x*(np.sqrt(Q)*(-nur*pdotvr - pdotn*xisq*(mucostheta - nucostheta)) - mucostheta*pdotn*xisq + mur*pdotvr*(np.sqrt(Q) + 1)) + Btilde**2*Jtildeprm_x*Sdotxi*(np.sqrt(Q)*(-nur*pdotvr - pdotn*xisq*(mucostheta - nucostheta)) - mucostheta*pdotn*xisq + mur*pdotvr*(np.sqrt(Q) + 1)) + Btilde*Btildeprm_x*Jtilde*Sdotxi*(2*np.sqrt(Q)*(-nur*pdotvr - pdotn*xisq*(mucostheta - nucostheta)) - 2*mucostheta*pdotn*xisq + 2*mur*pdotvr*(np.sqrt(Q) + 1))
    HsoTerm2prm_x = HsoTerm2aprm_x + HsoTerm2bprm_x - HsoTerm2cprm_x
    HsoTerm2coeffprm_x = expnu*(-Qprm_x - Qprm_x/(2*np.sqrt(Q)))/(Btilde**2*exp2mu*xisq*(np.sqrt(Q) + Q)**2) - expnu*xisqprm_x/(Btilde**2*exp2mu*xisq**2*(np.sqrt(Q) + Q)) + expnuprm_x/(Btilde**2*exp2mu*xisq*(np.sqrt(Q) + Q)) - exp2muprm_x*expnu/(Btilde**2*exp2mu**2*xisq*(np.sqrt(Q) + Q)) - 2*Btildeprm_x*expnu/(Btilde**3*exp2mu*xisq*(np.sqrt(Q) + Q))
    HsoTerm1prm_x = SdotSkerrhat*exp2nu*pdotxir*(-Btildeprm_x + expmu*expnuprm_x + expmuprm_x*expnu)/(Btilde**2*np.sqrt(Q)*expmu*xisq) + SdotSkerrhat*exp2nu*pdotxir*xisqprm_x*(Btilde - expmu*expnu)/(Btilde**2*np.sqrt(Q)*expmu*xisq**2) + SdotSkerrhat*exp2nu*pdotxirprm_x*(-Btilde + expmu*expnu)/(Btilde**2*np.sqrt(Q)*expmu*xisq) - SdotSkerrhat*exp2nu*expmuprm_x*pdotxir*(-Btilde + expmu*expnu)/(Btilde**2*np.sqrt(Q)*expmu**2*xisq) + SdotSkerrhat*exp2nuprm_x*pdotxir*(-Btilde + expmu*expnu)/(Btilde**2*np.sqrt(Q)*expmu*xisq) + SdotSkerrhatprm_x*exp2nu*pdotxir*(-Btilde + expmu*expnu)/(Btilde**2*np.sqrt(Q)*expmu*xisq) - Qprm_x*SdotSkerrhat*exp2nu*pdotxir*(-Btilde + expmu*expnu)/(2*Btilde**2*Q**(3/2)*expmu*xisq) - Btildeprm_x*SdotSkerrhat*exp2nu*pdotxir*(-2*Btilde + 2*expmu*expnu)/(Btilde**3*np.sqrt(Q)*expmu*xisq)
    Hsoprm_x = HsoTerm1prm_x + HsoTerm2*HsoTerm2coeffprm_x + HsoTerm2coeff*HsoTerm2prm_x
    Hdprm_x = Hdcoeff*Hdsumprm_x + Hdcoeffprm_x*Hdsum
    Hnsprm_x = np.sqrt(Hnsradicand)*alphaprm_x + betapsumprm_x + Hnsradicandprm_x*alpha/(2*np.sqrt(Hnsradicand))
    Hsprm_x = Hsoprm_x + Hssprm_x
    Heffprm_x = -Hdprm_x + Hnsprm_x + Hsprm_x + 4*dSS*eta*u**3*uprm_x*(S1x**2 + S1y**2 + S1z**2 + S2x**2 + S2y**2 + S2z**2)
    Hrealprm_x = Heffprm_x*eta/np.sqrt(2*eta*(Heff - 1) + 1)
    return np.array([Hrealprm_x])
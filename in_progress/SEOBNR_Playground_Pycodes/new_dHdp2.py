from __future__ import division
import numpy as np
def new_compute_dHdp2(m1, m2, EMgamma, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
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
    prTprm_p2 = csi2*n2
    phat3prm_p2 = n3*prTprm_p2*(1 - 1/csi1)
    phat2prm_p2 = n2*prTprm_p2*(1 - 1/csi1) + 1
    phat1prm_p2 = n1*prTprm_p2*(1 - 1/csi1)
    pdotxirprm_p2 = r*(phat1prm_p2*xi1 + phat2prm_p2*xi2 + phat3prm_p2*xi3)
    pdotnprm_p2 = n1*phat1prm_p2 + n2*phat2prm_p2 + n3*phat3prm_p2
    pdotvrprm_p2 = r*(phat1prm_p2*v1 + phat2prm_p2*v2 + phat3prm_p2*v3)
    pphiprm_p2 = pdotxirprm_p2
    DrSipn2prm_p2 = 2*Deltar*pdotn*pdotnprm_p2/Sigma
    Qprm_p2 = DrSipn2prm_p2 + 2*Qcoeff1*pdotxir*pdotxirprm_p2 + 2*Qcoeff2*pdotvr*pdotvrprm_p2
    Qminus1prm_p2 = Qprm_p2
    sigmacoeffTerm2prm_p2 = eta*(eta*r*(-882*DrSipn2prm_p2 + 204*Qminus1prm_p2 + r*(1620*DrSipn2*DrSipn2prm_p2 - 234*DrSipn2*Qminus1prm_p2 - 234*DrSipn2prm_p2*Qminus1)) + r*(-96*DrSipn2prm_p2 - 436*Qminus1prm_p2 + r*(36*DrSipn2*Qminus1prm_p2 + 36*DrSipn2prm_p2*Qminus1 - 90*Qminus1*Qminus1prm_p2)))/(144*r**2)
    sigmacoeffTerm1prm_p2 = eta*(-36*DrSipn2prm_p2 + 3*Qminus1prm_p2)/12
    sigmacoeffprm_p2 = sigmacoeffTerm1prm_p2 + sigmacoeffTerm2prm_p2
    sigmastarcoeffTerm2prm_p2 = eta*(eta*r*(-324*DrSipn2prm_p2 + 120*Qminus1prm_p2 + r*(720*DrSipn2*DrSipn2prm_p2 + Qminus1*(-126*DrSipn2prm_p2 - 3*Qminus1prm_p2) + Qminus1prm_p2*(-126*DrSipn2 - 3*Qminus1))) + r*(282*DrSipn2prm_p2 + Qminus1*r*(96*DrSipn2prm_p2 - 23*Qminus1prm_p2) + Qminus1prm_p2*r*(96*DrSipn2 - 23*Qminus1) - 206*Qminus1prm_p2))/(72*r**2)
    sigmastarcoeffTerm1prm_p2 = eta*(-30*DrSipn2prm_p2 + 4*Qminus1prm_p2)/12
    sigmastarcoeffprm_p2 = sigmastarcoeffTerm1prm_p2 + sigmastarcoeffTerm2prm_p2
    Deltasigmastar3prm_p2 = sigma3*sigmacoeffprm_p2 + sigmastar3*sigmastarcoeffprm_p2
    Deltasigmastar2prm_p2 = sigma2*sigmacoeffprm_p2 + sigmastar2*sigmastarcoeffprm_p2
    Deltasigmastar1prm_p2 = sigma1*sigmacoeffprm_p2 + sigmastar1*sigmastarcoeffprm_p2
    Sstar3prm_p2 = Deltasigmastar3prm_p2
    Sstar2prm_p2 = Deltasigmastar2prm_p2
    Sstar1prm_p2 = Deltasigmastar1prm_p2
    S3prm_p2 = Sstar3prm_p2
    S2prm_p2 = Sstar2prm_p2
    S1prm_p2 = Sstar1prm_p2
    Sstardotnprm_p2 = Sstar1prm_p2*n1 + Sstar2prm_p2*n2 + Sstar3prm_p2*n3
    SdotSkerrhatprm_p2 = S1prm_p2*Skerrhat1 + S2prm_p2*Skerrhat2 + S3prm_p2*Skerrhat3
    Sdotnprm_p2 = S1prm_p2*n1 + S2prm_p2*n2 + S3prm_p2*n3
    Sdotvprm_p2 = S1prm_p2*v1 + S2prm_p2*v2 + S3prm_p2*v3
    Sdotxiprm_p2 = S1prm_p2*xi1 + S2prm_p2*xi2 + S3prm_p2*xi3
    HdsumTerm2prm_p2 = 6*Sstardotn*Sstardotnprm_p2
    HdsumTerm1prm_p2 = 2*Sstar1*Sstar1prm_p2 + 2*Sstar2*Sstar2prm_p2 + 2*Sstar3*Sstar3prm_p2
    Hdsumprm_p2 = HdsumTerm1prm_p2 - HdsumTerm2prm_p2
    Q4prm_p2 = 8*eta*prT**3*prTprm_p2*u**2*(4 - 3*eta)
    gammappsumprm_p2 = 2*Deltar*pdotn*pdotnprm_p2/Sigma + 2*pdotvr*pdotvrprm_p2/(Sigma*sin2theta) + 2*Sigma*pdotxir*pdotxirprm_p2/(Lambt*sin2theta)
    Hnsradicandprm_p2 = Q4prm_p2 + gammappsumprm_p2
    betapsumprm_p2 = omegatilde*pphiprm_p2/Lambt
    HssTerm3prm_p2 = Btilde**2*(Sdotn*exp2mu*xisq*(-Qprm_p2 - Qprm_p2/(2*np.sqrt(Q))) + Sdotnprm_p2*exp2mu*xisq*(-np.sqrt(Q) - Q) + pdotvr*(-Jtilde*Sdotv*pdotnprm_p2 - Jtilde*Sdotvprm_p2*pdotn + Sdotn*pdotvrprm_p2 + Sdotnprm_p2*pdotvr) + pdotvrprm_p2*(-Jtilde*Sdotv*pdotn + Sdotn*pdotvr)) + expmu*expnu*pdotxir*(Btilde*Jtilde*Sdotxi*pdotnprm_p2 + Btilde*Jtilde*Sdotxiprm_p2*pdotn - Sdotn*expmu*expnu*pdotxirprm_p2 - Sdotnprm_p2*expmu*expnu*pdotxir) + expmu*expnu*pdotxirprm_p2*(Btilde*Jtilde*Sdotxi*pdotn - Sdotn*expmu*expnu*pdotxir)
    HssTerm3coeffprm_p2 = omegacostheta*(-Qprm_p2 - Qprm_p2/(2*np.sqrt(Q)))/(2*Btilde*exp2mu*expmu*expnu*(np.sqrt(Q) + Q)**2)
    HssTerm2prm_p2 = Btilde**2*xisq*(Jtilde*pdotn*(-Jtilde*Sdotv*pdotnprm_p2 - Jtilde*Sdotvprm_p2*pdotn + Sdotn*pdotvrprm_p2 + Sdotnprm_p2*pdotvr) + Jtilde*pdotnprm_p2*(-Jtilde*Sdotv*pdotn + Sdotn*pdotvr) + Sdotv*exp2mu*(Qprm_p2 + Qprm_p2/(2*np.sqrt(Q))) + Sdotvprm_p2*exp2mu*(np.sqrt(Q) + Q)) + expmu*pdotxir*(-Btilde*Sdotxi*expnu*pdotvrprm_p2 - Btilde*Sdotxiprm_p2*expnu*pdotvr + Sdotv*exp2nu*expmu*pdotxirprm_p2 + Sdotvprm_p2*exp2nu*expmu*pdotxir) + expmu*pdotxirprm_p2*(-Btilde*Sdotxi*expnu*pdotvr + Sdotv*exp2nu*expmu*pdotxir)
    HssTerm2coeffprm_p2 = Jtilde*omegar*(-Qprm_p2 - Qprm_p2/(2*np.sqrt(Q)))/(2*Btilde*exp2mu*expmu*expnu*xisq*(np.sqrt(Q) + Q)**2)
    HssTerm1prm_p2 = SdotSkerrhatprm_p2*omega
    Hssprm_p2 = HssTerm1prm_p2 + HssTerm2*HssTerm2coeffprm_p2 + HssTerm2coeff*HssTerm2prm_p2 + HssTerm3*HssTerm3coeffprm_p2 + HssTerm3coeff*HssTerm3prm_p2
    HsoTerm2cprm_p2 = Brtilde*Jtilde*Sdotv*expmu*expnu*pdotxirprm_p2*(np.sqrt(Q) + 1) + Brtilde*Jtilde*Sdotvprm_p2*expmu*expnu*pdotxir*(np.sqrt(Q) + 1) + Brtilde*Jtilde*Qprm_p2*Sdotv*expmu*expnu*pdotxir/(2*np.sqrt(Q))
    HsoTerm2bprm_p2 = Btilde*expmu*expnu*pdotxir*(2*np.sqrt(Q) + 1)*(Jtilde*Sdotvprm_p2*nur - Sdotnprm_p2*nucostheta*xisq) + Btilde*expmu*expnu*pdotxirprm_p2*(2*np.sqrt(Q) + 1)*(Jtilde*Sdotv*nur - Sdotn*nucostheta*xisq) + Btilde*Qprm_p2*expmu*expnu*pdotxir*(Jtilde*Sdotv*nur - Sdotn*nucostheta*xisq)/np.sqrt(Q)
    HsoTerm2aprm_p2 = Btilde**2*Jtilde*Sdotxi*(-np.sqrt(Q)*(nur*pdotvrprm_p2 + pdotnprm_p2*xisq*(mucostheta - nucostheta)) - mucostheta*pdotnprm_p2*xisq + mur*pdotvrprm_p2*(np.sqrt(Q) + 1) + Qprm_p2*mur*pdotvr/(2*np.sqrt(Q)) + Qprm_p2*(-nur*pdotvr - pdotn*xisq*(mucostheta - nucostheta))/(2*np.sqrt(Q))) + Btilde**2*Jtilde*Sdotxiprm_p2*(np.sqrt(Q)*(-nur*pdotvr - pdotn*xisq*(mucostheta - nucostheta)) - mucostheta*pdotn*xisq + mur*pdotvr*(np.sqrt(Q) + 1))
    HsoTerm2prm_p2 = HsoTerm2aprm_p2 + HsoTerm2bprm_p2 - HsoTerm2cprm_p2
    HsoTerm2coeffprm_p2 = expnu*(-Qprm_p2 - Qprm_p2/(2*np.sqrt(Q)))/(Btilde**2*exp2mu*xisq*(np.sqrt(Q) + Q)**2)
    HsoTerm1prm_p2 = SdotSkerrhat*exp2nu*pdotxirprm_p2*(-Btilde + expmu*expnu)/(Btilde**2*np.sqrt(Q)*expmu*xisq) + SdotSkerrhatprm_p2*exp2nu*pdotxir*(-Btilde + expmu*expnu)/(Btilde**2*np.sqrt(Q)*expmu*xisq) - Qprm_p2*SdotSkerrhat*exp2nu*pdotxir*(-Btilde + expmu*expnu)/(2*Btilde**2*Q**(3/2)*expmu*xisq)
    Hsoprm_p2 = HsoTerm1prm_p2 + HsoTerm2*HsoTerm2coeffprm_p2 + HsoTerm2coeff*HsoTerm2prm_p2
    Hdprm_p2 = Hdcoeff*Hdsumprm_p2
    Hnsprm_p2 = betapsumprm_p2 + Hnsradicandprm_p2*alpha/(2*np.sqrt(Hnsradicand))
    Hsprm_p2 = Hsoprm_p2 + Hssprm_p2
    Heffprm_p2 = -Hdprm_p2 + Hnsprm_p2 + Hsprm_p2
    Hrealprm_p2 = Heffprm_p2*eta/np.sqrt(2*eta*(Heff - 1) + 1)
    return np.array([Hrealprm_p2])
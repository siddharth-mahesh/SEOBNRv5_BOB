from __future__ import division
import numpy as np
def new_compute_dHdptheta(m1, m2, eta, r, theta, phi, pr, ptheta, pphi, S1r, S1theta, S1phi, S2r, S2theta, S2phi, KK, k0, k1, dSO, dSS, tortoise, EMgamma):
    M = m1+m2
    mu = m1*m2/M
    eta = mu/M
    u = 1/r
    sigmastar3 = m2/m1*S1phi+m1/m2*S2phi
    sigmastar2 = m2/m1*S1theta+m1/m2*S2theta
    sigmastar1 = m2/m1*S1r+m1/m2*S2r
    sigma3 = S1phi+S2phi##Sid:Snp.pinsshouldbeinpolardirections
    sigma2 = S1theta+S2theta
    sigma1 = S1r+S2r
    Skerr3 = sigma3##Sid:Snp.pinsshouldbeinpolardirections
    Skerr2 = sigma2
    Skerr1 = sigma1
    Skerrmag = np.sqrt(Skerr1*Skerr1+Skerr2*Skerr2+Skerr3*Skerr3)
    Skerrhat3 = Skerr3/Skerrmag
    Skerrhat2 = Skerr2/Skerrmag
    Skerrhat1 = Skerr1/Skerrmag
    a = Skerrmag
    n3 = 0
    n2 = 0
    n1 = 1
    e33 = Skerrhat3##Sid:insphericalpolar,nisjust(1,0,0)
    e32 = Skerrhat2
    e31 = Skerrhat1
    xi3 = e31*n2-e32*n1
    xi2 = e33*n1-e31*n3
    xi1 = e32*n3-e33*n2
    v3 = n1*xi2-n2*xi1##Sid:Testingoutatypohere
    v2 = n3*xi1-n1*xi3
    v1 = n2*xi3-n3*xi2
    costhetaBL = e31*n1+e32*n2+e33*n3
    sin2thetaBL = 1-costhetaBL*costhetaBL
    xisq = sin2thetaBL
    w2 = a*a+r*r
    Sigma = r*r+a*a*costhetaBL*costhetaBL
    sintheta = np.sin(theta)##Sid:costhetaBLsubstitution
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
    Lambt = w2*w2-a*a*Deltat*sin2thetaBL
    csi = np.sqrt(Deltar*Deltat)/w2##Sid:sin2thetaBLsubstitution
    csi1 = 1+(1-abs(1-tortoise))*(csi-1)
    csi2 = 1+(np.true_divide(1,2)-np.true_divide(1,2)*np.sign(np.true_divide(3,2)-tortoise))*(csi-1)
    prT = csi2*(pr*n1+ptheta*n2/r+pphi*n3/r/sintheta)
    phat3 = pphi/r/sintheta+prT*(1-1/csi1)*n3##Sid:Rewritewithsphericalmomenta
    phat2 = ptheta/r+prT*(1-1/csi1)*n2
    phat1 = pr+prT*(1-1/csi1)*n1
    pdotxir = (phat1*xi1+phat2*xi2+phat3*xi3)*r##Sid:Rewritewithsphericalmomenta
    pdotn = phat1*n1+phat2*n2+phat3*n3
    pdotvr = (phat1*v1+phat2*v2+phat3*v3)*r
    Qcoeff2 = 1/(Sigma*sin2thetaBL)
    Qcoeff1 = Sigma/(Lambt*sin2thetaBL)##Sid:sin2thetaBLsubstitution
    DrSipn2 = Deltar*pdotn*pdotn/Sigma##Sid:sin2thetaBLsubstitution
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
    Lambtprime = 4*(a*a+r*r)*r-2*a*a*Deltatprime*sin2thetaBL
    mucostheta = a*a*costhetaBL/Sigma##Sid:sin2thetaBLsubstitution
    nucostheta = a*a*w2*costhetaBL*(w2-Deltat)/(Lambt*Sigma)##Sid:costhetaBLsubstitution
    omegacostheta = -2*a*a*costhetaBL*Deltat*omegatilde/(Lambt*Lambt)##Sid:costhetaBLsubstitution
    mur = r/Sigma-1/np.sqrt(Deltar)##Sid:costhetaBLsubstitution
    nur = r/Sigma+w2*(w2*Deltatprime-4*r*Deltat)/(2*Lambt*Deltat)
    omegar = (Lambt*omegatildeprime-Lambtprime*omegatilde)/(Lambt*Lambt)
    dSO = -74.71-156.*eta+627.5*eta*eta
    sigmacoeffTerm3 = eta*dSO*u*u*u
    sigmacoeffTerm2 = eta/(144*r*r)*(-896+r*(-436*Qminus1-96*DrSipn2+r*(-45*Qminus1*Qminus1+36*Qminus1*DrSipn2))+eta*(-336+r*(204*Qminus1-882*DrSipn2+r*(810*DrSipn2*DrSipn2-234*Qminus1*DrSipn2))))
    sigmacoeffTerm1 = eta/12*(-8*u+3*Qminus1-36*DrSipn2)
    sigmacoeff = sigmacoeffTerm1+sigmacoeffTerm2+sigmacoeffTerm3##Sid:Again,multiplybyuinstead
    sigmastarcoeffTerm2 = eta/(72*r*r)*(706+r*(-206*Qminus1+282*DrSipn2+r*Qminus1*(96*DrSipn2-23*Qminus1))+eta*(-54+r*(120*Qminus1-324*DrSipn2+r*(360*DrSipn2*DrSipn2+Qminus1*(-126*DrSipn2-3*Qminus1)))))
    sigmastarcoeffTerm1 = eta/12*(14*u+4*Qminus1-30*DrSipn2)
    sigmastarcoeff = sigmastarcoeffTerm1+sigmastarcoeffTerm2##Sid:Multiplybyuoncemore
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
    Hdcoeff = np.true_divide(1,2)*(u*u*u)
    Q4 = 2*prT*prT*prT*prT*u*u*(4-3*eta)*eta##Sid:Whydividebyrwhenyoucanmultiplybyu?
    gammappsum = Deltar/Sigma*pdotn*pdotn+1/Sigma*pdotvr*pdotvr/sin2thetaBL+Sigma/Lambt/sin2thetaBL*pdotxir*pdotxir
    Hnsradicand = 1+gammappsum+Q4##Sid:sin2thetaBLsubstitution
    alpha = np.sqrt(Deltat*Sigma/Lambt)
    betapsum = omegatilde*pdotxir/Lambt
    HssTerm3 = expmu*expnu*pdotxir*(Jtilde*pdotn*Sdotxi*Btilde-expmu*expnu*pdotxir*Sdotn)+(pdotvr*(Sdotn*pdotvr-Jtilde*pdotn*Sdotv)-exp2mu*(np.sqrt(Q)+Q)*Sdotn*xisq)*Btilde*Btilde##Sid:Changedtopdotxirtoeliminatepphi
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
    Heff = Hs+Hns-Hd+dSS*eta*u*u*u*u*(S1r*S1r+S1theta*S1theta+S1phi*S1phi+S2r*S2r+S2theta*S2theta+S2phi*S2phi)
    prTprm_ptheta = csi2*n2/r
    phat3prm_ptheta = n3*prTprm_ptheta*(1 - 1/csi1)
    phat2prm_ptheta = n2*prTprm_ptheta*(1 - 1/csi1) + 1/r
    phat1prm_ptheta = n1*prTprm_ptheta*(1 - 1/csi1)
    pdotxirprm_ptheta = r*(phat1prm_ptheta*xi1 + phat2prm_ptheta*xi2 + phat3prm_ptheta*xi3)
    pdotnprm_ptheta = n1*phat1prm_ptheta + n2*phat2prm_ptheta + n3*phat3prm_ptheta
    pdotvrprm_ptheta = r*(phat1prm_ptheta*v1 + phat2prm_ptheta*v2 + phat3prm_ptheta*v3)
    DrSipn2prm_ptheta = 2*Deltar*pdotn*pdotnprm_ptheta/Sigma
    Qprm_ptheta = DrSipn2prm_ptheta + 2*Qcoeff1*pdotxir*pdotxirprm_ptheta + 2*Qcoeff2*pdotvr*pdotvrprm_ptheta
    Qminus1prm_ptheta = Qprm_ptheta
    sigmacoeffTerm2prm_ptheta = eta*(eta*r*(-882*DrSipn2prm_ptheta + 204*Qminus1prm_ptheta + r*(1620*DrSipn2*DrSipn2prm_ptheta - 234*DrSipn2*Qminus1prm_ptheta - 234*DrSipn2prm_ptheta*Qminus1)) + r*(-96*DrSipn2prm_ptheta - 436*Qminus1prm_ptheta + r*(36*DrSipn2*Qminus1prm_ptheta + 36*DrSipn2prm_ptheta*Qminus1 - 90*Qminus1*Qminus1prm_ptheta)))/(144*r**2)
    sigmacoeffTerm1prm_ptheta = eta*(-36*DrSipn2prm_ptheta + 3*Qminus1prm_ptheta)/12
    sigmacoeffprm_ptheta = sigmacoeffTerm1prm_ptheta + sigmacoeffTerm2prm_ptheta
    sigmastarcoeffTerm2prm_ptheta = eta*(eta*r*(-324*DrSipn2prm_ptheta + 120*Qminus1prm_ptheta + r*(720*DrSipn2*DrSipn2prm_ptheta + Qminus1*(-126*DrSipn2prm_ptheta - 3*Qminus1prm_ptheta) + Qminus1prm_ptheta*(-126*DrSipn2 - 3*Qminus1))) + r*(282*DrSipn2prm_ptheta + Qminus1*r*(96*DrSipn2prm_ptheta - 23*Qminus1prm_ptheta) + Qminus1prm_ptheta*r*(96*DrSipn2 - 23*Qminus1) - 206*Qminus1prm_ptheta))/(72*r**2)
    sigmastarcoeffTerm1prm_ptheta = eta*(-30*DrSipn2prm_ptheta + 4*Qminus1prm_ptheta)/12
    sigmastarcoeffprm_ptheta = sigmastarcoeffTerm1prm_ptheta + sigmastarcoeffTerm2prm_ptheta
    Deltasigmastar3prm_ptheta = sigma3*sigmacoeffprm_ptheta + sigmastar3*sigmastarcoeffprm_ptheta
    Deltasigmastar2prm_ptheta = sigma2*sigmacoeffprm_ptheta + sigmastar2*sigmastarcoeffprm_ptheta
    Deltasigmastar1prm_ptheta = sigma1*sigmacoeffprm_ptheta + sigmastar1*sigmastarcoeffprm_ptheta
    Sstar3prm_ptheta = Deltasigmastar3prm_ptheta
    Sstar2prm_ptheta = Deltasigmastar2prm_ptheta
    Sstar1prm_ptheta = Deltasigmastar1prm_ptheta
    S3prm_ptheta = Sstar3prm_ptheta
    S2prm_ptheta = Sstar2prm_ptheta
    S1prm_ptheta = Sstar1prm_ptheta
    Sstardotnprm_ptheta = Sstar1prm_ptheta*n1 + Sstar2prm_ptheta*n2 + Sstar3prm_ptheta*n3
    SdotSkerrhatprm_ptheta = S1prm_ptheta*Skerrhat1 + S2prm_ptheta*Skerrhat2 + S3prm_ptheta*Skerrhat3
    Sdotnprm_ptheta = S1prm_ptheta*n1 + S2prm_ptheta*n2 + S3prm_ptheta*n3
    Sdotvprm_ptheta = S1prm_ptheta*v1 + S2prm_ptheta*v2 + S3prm_ptheta*v3
    Sdotxiprm_ptheta = S1prm_ptheta*xi1 + S2prm_ptheta*xi2 + S3prm_ptheta*xi3
    HdsumTerm2prm_ptheta = 6*Sstardotn*Sstardotnprm_ptheta
    HdsumTerm1prm_ptheta = 2*Sstar1*Sstar1prm_ptheta + 2*Sstar2*Sstar2prm_ptheta + 2*Sstar3*Sstar3prm_ptheta
    Hdsumprm_ptheta = HdsumTerm1prm_ptheta - HdsumTerm2prm_ptheta
    Q4prm_ptheta = 8*eta*prT**3*prTprm_ptheta*u**2*(-3*eta + 4)
    gammappsumprm_ptheta = 2*Deltar*pdotn*pdotnprm_ptheta/Sigma + 2*pdotvr*pdotvrprm_ptheta/(Sigma*sin2thetaBL) + 2*Sigma*pdotxir*pdotxirprm_ptheta/(Lambt*sin2thetaBL)
    Hnsradicandprm_ptheta = Q4prm_ptheta + gammappsumprm_ptheta
    betapsumprm_ptheta = omegatilde*pdotxirprm_ptheta/Lambt
    HssTerm3prm_ptheta = Btilde**2*(Sdotn*exp2mu*xisq*(-Qprm_ptheta - Qprm_ptheta/(2*np.sqrt(Q))) + Sdotnprm_ptheta*exp2mu*xisq*(-np.sqrt(Q) - Q) + pdotvr*(-Jtilde*Sdotv*pdotnprm_ptheta - Jtilde*Sdotvprm_ptheta*pdotn + Sdotn*pdotvrprm_ptheta + Sdotnprm_ptheta*pdotvr) + pdotvrprm_ptheta*(-Jtilde*Sdotv*pdotn + Sdotn*pdotvr)) + expmu*expnu*pdotxir*(Btilde*Jtilde*Sdotxi*pdotnprm_ptheta + Btilde*Jtilde*Sdotxiprm_ptheta*pdotn - Sdotn*expmu*expnu*pdotxirprm_ptheta - Sdotnprm_ptheta*expmu*expnu*pdotxir) + expmu*expnu*pdotxirprm_ptheta*(Btilde*Jtilde*Sdotxi*pdotn - Sdotn*expmu*expnu*pdotxir)
    HssTerm3coeffprm_ptheta = omegacostheta*(-Qprm_ptheta - Qprm_ptheta/(2*np.sqrt(Q)))/(2*Btilde*exp2mu*expmu*expnu*(np.sqrt(Q) + Q)**2)
    HssTerm2prm_ptheta = Btilde**2*xisq*(Jtilde*pdotn*(-Jtilde*Sdotv*pdotnprm_ptheta - Jtilde*Sdotvprm_ptheta*pdotn + Sdotn*pdotvrprm_ptheta + Sdotnprm_ptheta*pdotvr) + Jtilde*pdotnprm_ptheta*(-Jtilde*Sdotv*pdotn + Sdotn*pdotvr) + Sdotv*exp2mu*(Qprm_ptheta + Qprm_ptheta/(2*np.sqrt(Q))) + Sdotvprm_ptheta*exp2mu*(np.sqrt(Q) + Q)) + expmu*pdotxir*(-Btilde*Sdotxi*expnu*pdotvrprm_ptheta - Btilde*Sdotxiprm_ptheta*expnu*pdotvr + Sdotv*exp2nu*expmu*pdotxirprm_ptheta + Sdotvprm_ptheta*exp2nu*expmu*pdotxir) + expmu*pdotxirprm_ptheta*(-Btilde*Sdotxi*expnu*pdotvr + Sdotv*exp2nu*expmu*pdotxir)
    HssTerm2coeffprm_ptheta = Jtilde*omegar*(-Qprm_ptheta - Qprm_ptheta/(2*np.sqrt(Q)))/(2*Btilde*exp2mu*expmu*expnu*xisq*(np.sqrt(Q) + Q)**2)
    HssTerm1prm_ptheta = SdotSkerrhatprm_ptheta*omega
    Hssprm_ptheta = HssTerm1prm_ptheta + HssTerm2*HssTerm2coeffprm_ptheta + HssTerm2coeff*HssTerm2prm_ptheta + HssTerm3*HssTerm3coeffprm_ptheta + HssTerm3coeff*HssTerm3prm_ptheta
    HsoTerm2cprm_ptheta = Brtilde*Jtilde*Sdotv*expmu*expnu*pdotxirprm_ptheta*(np.sqrt(Q) + 1) + Brtilde*Jtilde*Sdotvprm_ptheta*expmu*expnu*pdotxir*(np.sqrt(Q) + 1) + Brtilde*Jtilde*Qprm_ptheta*Sdotv*expmu*expnu*pdotxir/(2*np.sqrt(Q))
    HsoTerm2bprm_ptheta = Btilde*expmu*expnu*pdotxir*(2*np.sqrt(Q) + 1)*(Jtilde*Sdotvprm_ptheta*nur - Sdotnprm_ptheta*nucostheta*xisq) + Btilde*expmu*expnu*pdotxirprm_ptheta*(2*np.sqrt(Q) + 1)*(Jtilde*Sdotv*nur - Sdotn*nucostheta*xisq) + Btilde*Qprm_ptheta*expmu*expnu*pdotxir*(Jtilde*Sdotv*nur - Sdotn*nucostheta*xisq)/np.sqrt(Q)
    HsoTerm2aprm_ptheta = Btilde**2*Jtilde*Sdotxi*(-np.sqrt(Q)*(nur*pdotvrprm_ptheta + pdotnprm_ptheta*xisq*(mucostheta - nucostheta)) - mucostheta*pdotnprm_ptheta*xisq + mur*pdotvrprm_ptheta*(np.sqrt(Q) + 1) + Qprm_ptheta*mur*pdotvr/(2*np.sqrt(Q)) + Qprm_ptheta*(-nur*pdotvr - pdotn*xisq*(mucostheta - nucostheta))/(2*np.sqrt(Q))) + Btilde**2*Jtilde*Sdotxiprm_ptheta*(np.sqrt(Q)*(-nur*pdotvr - pdotn*xisq*(mucostheta - nucostheta)) - mucostheta*pdotn*xisq + mur*pdotvr*(np.sqrt(Q) + 1))
    HsoTerm2prm_ptheta = HsoTerm2aprm_ptheta + HsoTerm2bprm_ptheta - HsoTerm2cprm_ptheta
    HsoTerm2coeffprm_ptheta = expnu*(-Qprm_ptheta - Qprm_ptheta/(2*np.sqrt(Q)))/(Btilde**2*exp2mu*xisq*(np.sqrt(Q) + Q)**2)
    HsoTerm1prm_ptheta = SdotSkerrhat*exp2nu*pdotxirprm_ptheta*(-Btilde + expmu*expnu)/(Btilde**2*np.sqrt(Q)*expmu*xisq) + SdotSkerrhatprm_ptheta*exp2nu*pdotxir*(-Btilde + expmu*expnu)/(Btilde**2*np.sqrt(Q)*expmu*xisq) - Qprm_ptheta*SdotSkerrhat*exp2nu*pdotxir*(-Btilde + expmu*expnu)/(2*Btilde**2*Q**(3/2)*expmu*xisq)
    Hsoprm_ptheta = HsoTerm1prm_ptheta + HsoTerm2*HsoTerm2coeffprm_ptheta + HsoTerm2coeff*HsoTerm2prm_ptheta
    Hdprm_ptheta = Hdcoeff*Hdsumprm_ptheta
    Hnsprm_ptheta = betapsumprm_ptheta + Hnsradicandprm_ptheta*alpha/(2*np.sqrt(Hnsradicand))
    Hsprm_ptheta = Hsoprm_ptheta + Hssprm_ptheta
    Heffprm_ptheta = -Hdprm_ptheta + Hnsprm_ptheta + Hsprm_ptheta
    Hrealprm_ptheta = Heffprm_ptheta*eta/np.sqrt(2*eta*(Heff - 1) + 1)
    return np.array([Hrealprm_ptheta])
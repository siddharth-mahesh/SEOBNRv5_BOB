from __future__ import division
import numpy as np
def ham_first_derivs(m1, m2, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z):
    EMgamma = 0.577215664901532860606512090082402431
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
    L3 = x*p2-y*p1
    L2 = z*p1-x*p3
    L1 = y*p3-z*p2
    Lnorm = np.sqrt(L1*L1+L2*L2+L3*L3)
    Lhat3 = L3/Lnorm
    Lhat2 = L2/Lnorm
    Lhat1 = L1/Lnorm
    S2dotLhat = S2x*Lhat1+S2y*Lhat2+S2z*Lhat3
    S1dotLhat = S1x*Lhat1+S1y*Lhat2+S1z*Lhat3
    S1perp3 = S1z-S1dotLhat*Lhat3
    S1perp2 = S1y-S1dotLhat*Lhat2
    S1perp1 = S1x-S1dotLhat*Lhat1
    S2perp3 = S2z-S2dotLhat*Lhat3
    S2perp2 = S2y-S2dotLhat*Lhat2
    S2perp1 = S2x-S2dotLhat*Lhat1
    Sperp3 = S1perp3+S2perp3
    Sperp2 = S1perp2+S2perp2
    Sperp1 = S1perp1+S2perp1
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
    Dinv = 1+np.log1p(6*eta*u*u+2*(26-3*eta)*eta*u*u*u)
    omegatilde = 2*a*r
    chi = (Skerr1*Lhat1+Skerr2*Lhat2+Skerr3*Lhat3)/(1-2*eta)+np.true_divide(1,2)*(Sperp1*Skerr1+Sperp2*Skerr2+Sperp3*Skerr3)/(Skerrmag*(1.-2.*eta))
    Kchi0 = 267.788*eta*eta*eta-126.687*eta*eta+10.2573*eta+1.7336
    K = -59.1658*chi*chi*chi*eta*eta*eta-0.426958*chi*chi*chi*eta+1.43659*chi*chi*chi+31.1746*chi*chi*eta*eta*eta+6.16466*chi*chi*eta*eta-1.38086*chi*chi-27.5201*chi*eta*eta*eta+17.3736*chi*eta*eta+2.26831*chi*eta-1.62045*chi+Kchi0
    etaKminus1 = eta*K-1
    Delta0 = K*(eta*K-2)
    Delta1 = -2*etaKminus1*(K+Delta0)
    Delta2 = np.true_divide(1,2)*Delta1*(Delta1-4*etaKminus1)-a*a*etaKminus1*etaKminus1*Delta0
    Delta3 = -np.true_divide(1,3)*Delta1*Delta1*Delta1+etaKminus1*Delta1*Delta1+Delta2*Delta1-2*etaKminus1*(Delta2-etaKminus1)-a*a*etaKminus1*etaKminus1*Delta1
    Delta4 = np.true_divide(1,12)*(6*a*a*(Delta1*Delta1-2*Delta2)*etaKminus1*etaKminus1+3*Delta1*Delta1*Delta1*Delta1-8*etaKminus1*Delta1*Delta1*Delta1-12*Delta2*Delta1*Delta1+12*(2*etaKminus1*Delta2+Delta3)*Delta1+12*(np.true_divide(94,3)-np.true_divide(41,32)*np.pi*np.pi)*etaKminus1*etaKminus1+6*(Delta2*Delta2-4*Delta3*etaKminus1))
    Delta5 = etaKminus1*etaKminus1*(np.true_divide(-4237,60)+np.true_divide(128,5)*EMgamma+np.true_divide(2275,512)*np.pi*np.pi-np.true_divide(1,3)*a*a*(Delta1*Delta1*Delta1-3*Delta1*Delta2+3*Delta3)-(Delta1*Delta1*Delta1*Delta1*Delta1-5*Delta1*Delta1*Delta1*Delta2+5*Delta1*Delta2*Delta2+5*Delta1*Delta1*Delta3-5*Delta2*Delta3-5*Delta1*Delta4)/(5*etaKminus1*etaKminus1)+(Delta1*Delta1*Delta1*Delta1-4*Delta1*Delta1*Delta2+2*Delta2*Delta2+4*Delta1*Delta3-4*Delta4)/(2*etaKminus1)+np.true_divide(256,5)*np.log(2)+(np.true_divide(41,32)*np.pi*np.pi-np.true_divide(221,6))*eta)
    Delta5l = etaKminus1*etaKminus1*np.true_divide(64,5)
    logarg = u*(Delta1+u*(Delta2+u*(Delta3+u*(Delta4+u*(Delta5+Delta5l*np.log(u))))))
    Deltaucalib = 1+eta*(Delta0+np.log1p(logarg))
    Deltaucalibprime = -eta*u*u*(Delta1+u*(2*Delta2+u*(3*Delta3+u*(4*Delta4+u*(5*Delta5+Delta5l*(1+5*np.log(u)))))))/(1+logarg)
    Deltaubar = a*a*u*u+(2*u+1/etaKminus1)/etaKminus1
    Deltaubarprime = -2*a*a*u*u*u-2*u*u/(etaKminus1)
    Deltau = np.abs(Deltaubar*Deltaucalib)
    Deltauprime = Deltaubarprime*Deltaucalib+Deltaubar*Deltaucalibprime
    Deltatprime = 2*r*Deltau+r*r*Deltauprime
    Deltat = r*r*Deltau
    Deltar = Deltat*Dinv
    Lambt = np.abs(w2*w2-a*a*Deltat*sin2theta)
    csi = np.sqrt(Deltar*Deltat)/w2
    csi1 = 1+(1-np.abs(1-tortoise))*(csi-1)
    csi2 = 1+(np.true_divide(1,2)-np.true_divide(1,2)*np.sign(np.true_divide(3,2)-tortoise))*(csi-1)
    prT = csi2*(p1*n1+p2*n2+p3*n3)
    phat3 = p3-prT*(1-1/csi1)*n3
    phat2 = p2-prT*(1-1/csi1)*n2
    phat1 = p1-prT*(1-1/csi1)*n1
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
    Lambtprime = 4*(a*a+r*r)*r-a*a*Deltatprime*sin2theta
    mucostheta = a*a*costheta/Sigma
    nucostheta = a*a*w2*costheta*(w2-Deltat)/(Lambt*Sigma)
    omegacostheta = -2*a*a*costheta*Deltat*omegatilde/(Lambt*Lambt)
    mur = r/Sigma-1/np.sqrt(Deltar)
    nur = r/Sigma+w2*(w2*Deltatprime-4*r*Deltat)/(2*Lambt*Deltat)
    omegar = (Lambt*omegatildeprime-Lambtprime*omegatilde)/(Lambt*Lambt)
    dSO = 147.481*chi*chi*chi*eta*eta-568.651*chi*chi*chi*eta+66.1987*chi*chi*chi-343.313*chi*chi*eta+2495.29*chi*eta*eta-44.5324
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
    Hss = -0.5*u3*(sx*sx+sy*sy+sz*sz-3.*sn*sn)
    Hns = betapsum+alpha*np.sqrt(Hnsradicand)
    Hs = omega*SdotSkerrhat+Hso+Hsonl+omegar*Hwr+omegacos*Hwcos
    dSS = 528.511*chi*chi*chi*eta*eta-41.0003*chi*chi*chi*eta+1161.78*chi*chi*eta*eta*eta-326.325*chi*chi*eta*eta+37.1964*chi*eta+706.958*eta*eta*eta-36.0272*eta+6.06807
    Hcalib = dSS*eta*u*u*u*u*(S1x*S1x+S1y*S1y+S1z*S1z+S2x*S2x+S2y*S2y+S2z*S2z)
    Heff = Hs+Hns+Hss+Hcalib
    Hreal = np.sqrt(1+2*eta*(Heff-1))

    Hss_prm = u3*(3.0*sn*sn_prm - 1.0*sx*sx_prm - 1.0*sy*sy_prm - 1.0*sz*sz_prm) + u3_prm*(1.5*sn**2 - 0.5*sx**2 - 0.5*sy**2 - 0.5*sz**2)
    Hs_prm = Hsonl_prm + Hso_prm + Hwcos*omegacos_prm + Hwcos_prm*omegacos + Hwr_prm*omegar
    Heff_prm = Hs_prm
    Hreal_prm = Heff_prm*eta/np.sqrt(eta*(2*Heff - 2) + 1)
    return np.array([Hreal_prmx1, Hreal_prmx2, Hreal_prmx3, Hreal_prmp1, Hreal_prmp2, Hreal_prmp3, Hreal_prmS1x, Hreal_prmS1y, Hreal_prmS1z, Hreal_prmS2x, Hreal_prmS2y, Hreal_prmS2z, csi_prmx, csi_prmy, csi_prmz])

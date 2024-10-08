import numpy as np
def compute_v4P_Hreal(m1=23., m2=23., tortoise=1, x=-6.6254467465750873e+00, y=2.3112723896656919e+00, z=1.4784770976195734e+00, p1=-1.3939842737500083e-01, p2=-4.5283712451822244e-01, p3=1.0873416734209107e-02, S1x=-8.4581203140832748e-04, S1y=9.1636179969106744e-03, S1z=-1.6771127425826671e-03, S2x=-9.5296417537494696e-02, S2y=1.8673844000169849e-02, S2z=3.0840401173698895e-02):
    EMgamma = 0.577215664901532860606512090082402431
    print('EMgamma  = %.16e'% EMgamma )
    M=m1+m2
    print('M = %.16e'% M)
    mu=m1*m2/M
    print('mu = %.16e'% mu)
    eta=mu/M
    print('eta = %.16e'% eta)
    r=np.sqrt(x*x+y*y+z*z)
    print('r = %.16e'% r)
    u=1/r
    print('u = %.16e'% u)
    sigmastar3=(m2/m1*S1z+m1/m2*S2z)/M/M
    print('sigmastar3 = %.16e'% sigmastar3)
    sigmastar2 = (m2/m1*S1y + m1/m2*S2y)/M/M
    print('sigmastar2  = %.16e'% sigmastar2 )
    sigmastar1 = (m2/m1*S1x + m1/m2*S2x)/M/M
    print('sigmastar1  = %.16e'% sigmastar1 )
    sigma3=(S1z+S2z)/M/M
    print('sigma3 = %.16e'% sigma3)
    sigma2 = (S1y + S2y)/M/M
    print('sigma2  = %.16e'% sigma2 )
    sigma1 = (S1x + S2x)/M/M
    print('sigma1  = %.16e'% sigma1 )
    Skerr3=sigma3
    print('Skerr3 = %.16e'% Skerr3)
    Skerr2 = sigma2
    print('Skerr2  = %.16e'% Skerr2 )
    Skerr1 = sigma1
    print('Skerr1  = %.16e'% Skerr1 )
    Skerrmag=np.sqrt(Skerr1*Skerr1+Skerr2*Skerr2+Skerr3*Skerr3)
    print('Skerrmag = %.16e'% Skerrmag)
    Skerrhat3=Skerr3/Skerrmag
    print('Skerrhat3 = %.16e'% Skerrhat3)
    Skerrhat2 = Skerr2/Skerrmag
    print('Skerrhat2  = %.16e'% Skerrhat2 )
    Skerrhat1 = Skerr1/Skerrmag
    print('Skerrhat1  = %.16e'% Skerrhat1 )
    a=Skerrmag
    print('a = %.16e'% a)
    L3=x*p2-y*p1
    print('L3 = %.16e'% L3)
    L2 = z*p1 - x*p3
    print('L2  = %.16e'% L2 )
    L1 = y*p3 - z*p2
    print('L1  = %.16e'% L1 )
    Lnorm=np.sqrt(L1*L1+L2*L2+L3*L3)
    print('Lnorm = %.16e'% Lnorm)
    Lhat3=L3/Lnorm
    print('Lhat3 = %.16e'% Lhat3)
    Lhat2 = L2/Lnorm
    print('Lhat2  = %.16e'% Lhat2 )
    Lhat1 = L1/Lnorm
    print('Lhat1  = %.16e'% Lhat1 )
    S2dotLhat=S2x*Lhat1+S2y*Lhat2+S2z*Lhat3
    print('S2dotLhat = %.16e'% S2dotLhat)
    S1dotLhat = S1x*Lhat1 + S1y*Lhat2 + S1z*Lhat3
    print('S1dotLhat  = %.16e'% S1dotLhat )
    S1perp3=S1z-S1dotLhat*Lhat3
    print('S1perp3 = %.16e'% S1perp3)
    S1perp2 = S1y - S1dotLhat*Lhat2
    print('S1perp2  = %.16e'% S1perp2 )
    S1perp1 = S1x - S1dotLhat*Lhat1
    print('S1perp1  = %.16e'% S1perp1 )
    S2perp3=S2z-S2dotLhat*Lhat3
    print('S2perp3 = %.16e'% S2perp3)
    S2perp2 = S2y - S2dotLhat*Lhat2
    print('S2perp2  = %.16e'% S2perp2 )
    S2perp1 = S2x - S2dotLhat*Lhat1
    print('S2perp1  = %.16e'% S2perp1 )
    Sperp3=S1perp3+S2perp3
    print('Sperp3 = %.16e'% Sperp3)
    Sperp2 = S1perp2 + S2perp2
    print('Sperp2  = %.16e'% Sperp2 )
    Sperp1 = S1perp1 + S2perp1
    print('Sperp1  = %.16e'% Sperp1 )
    n3=z/r
    print('n3 = %.16e'% n3)
    n2 = y/r
    print('n2  = %.16e'% n2 )
    n1 = x/r
    print('n1  = %.16e'% n1 )
    TINYDOUBLE=1e-100
    print('TINYDOUBLE = %.16e'% TINYDOUBLE)
    condition_e3prov_lhs = a
    print('condition_e3prov_lhs  = %.16e'% condition_e3prov_lhs )
    condition_e3prov_rhs = 1e-16
    print('condition_e3prov_rhs  = %.16e'% condition_e3prov_rhs )
    e3prov_gt_bound = np.divide(1,2)*(condition_e3prov_lhs - condition_e3prov_rhs + np.abs(condition_e3prov_lhs - condition_e3prov_rhs))/(condition_e3prov_lhs - condition_e3prov_rhs - TINYDOUBLE)
    print('e3prov_gt_bound  = %.16e'% e3prov_gt_bound )
    e3prov_leq_bound = np.divide(1,2)*(condition_e3prov_lhs - condition_e3prov_rhs - TINYDOUBLE - np.abs(condition_e3prov_lhs - condition_e3prov_rhs - TINYDOUBLE))/(condition_e3prov_lhs - condition_e3prov_rhs - TINYDOUBLE)
    print('e3prov_leq_bound  = %.16e'% e3prov_leq_bound )
    e3prov1 = Skerrhat1*e3prov_gt_bound + e3prov_leq_bound/np.sqrt(3.)
    print('e3prov1  = %.16e'% e3prov1 )
    e3prov2 = Skerrhat2*e3prov_gt_bound + e3prov_leq_bound/np.sqrt(3.)
    print('e3prov2  = %.16e'% e3prov2 )
    e3prov3 = Skerrhat3*e3prov_gt_bound + e3prov_leq_bound/np.sqrt(3.)
    print('e3prov3  = %.16e'% e3prov3 )
    lambdavec3=Lhat1*n2-Lhat2*n3
    print('lambdavec3 = %.16e'% lambdavec3)
    lambdavec2 = Lhat3*n1 - Lhat1*n3
    print('lambdavec2  = %.16e'% lambdavec2 )
    lambdavec1 = Lhat2*n3 - Lhat3*n2
    print('lambdavec1  = %.16e'% lambdavec1 )
    lambdavecnorm=np.sqrt(lambdavec1*lambdavec1+lambdavec2*lambdavec2+lambdavec3*lambdavec3)
    print('lambdavecnorm = %.16e'% lambdavecnorm)
    lambdahat3=lambdavec3/lambdavecnorm
    print('lambdahat3 = %.16e'% lambdahat3)
    lambdahat2 = lambdavec2/lambdavecnorm
    print('lambdahat2  = %.16e'% lambdahat2 )
    lambdahat1 = lambdavec1/lambdavecnorm
    print('lambdahat1  = %.16e'% lambdahat1 )
    lambdahat_dot_e3prov=lambdahat1*e3prov1+lambdahat2*e3prov2+lambdahat3*e3prov3
    print('lambdahat_dot_e3prov = %.16e'% lambdahat_dot_e3prov)
    lambdahat_cross_e3prov3=lambdahat1*e3prov2-lambdahat2*e3prov1
    print('lambdahat_cross_e3prov3 = %.16e'% lambdahat_cross_e3prov3)
    lambdahat_cross_e3prov2 = lambdahat3*e3prov1 - lambdahat1*e3prov3
    print('lambdahat_cross_e3prov2  = %.16e'% lambdahat_cross_e3prov2 )
    lambdahat_cross_e3prov1 = lambdahat2*e3prov3 - lambdahat3*e3prov2
    print('lambdahat_cross_e3prov1  = %.16e'% lambdahat_cross_e3prov1 )
    e3prov_dot_n=e3prov1*n1+e3prov2*n2+e3prov3*n3
    print('e3prov_dot_n = %.16e'% e3prov_dot_n)
    cos_0_1_deg=0.9999983800004374
    print('cos_0_1_deg = %.16e'% cos_0_1_deg)
    sin_0_1_deg = 0.0017999990280001574
    print('sin_0_1_deg  = %.16e'% sin_0_1_deg )
    condition_e3_lhs = 1 - np.abs(e3prov_dot_n)
    print('condition_e3_lhs  = %.16e'% condition_e3_lhs )
    condition_e3_rhs = 1e-8
    print('condition_e3_rhs  = %.16e'% condition_e3_rhs )
    e3_gt_bound = np.divide(1,2)*(condition_e3_lhs - condition_e3_rhs + np.abs(condition_e3_lhs - condition_e3_rhs))/(condition_e3_lhs - condition_e3_rhs - TINYDOUBLE)
    print('e3_gt_bound  = %.16e'% e3_gt_bound )
    e3_leq_bound = np.divide(1,2)*(condition_e3_lhs - condition_e3_rhs - TINYDOUBLE - np.abs(condition_e3_lhs - condition_e3_rhs - TINYDOUBLE))/(condition_e3_lhs - condition_e3_rhs - TINYDOUBLE)
    print('e3_leq_bound  = %.16e'% e3_leq_bound )
    e31 = e3prov1*e3_gt_bound + (e3prov1*cos_0_1_deg + lambdahat_cross_e3prov1*sin_0_1_deg + lambdahat1*lambdahat_dot_e3prov*(1 - cos_0_1_deg))*e3_leq_bound
    print('e31  = %.16e'% e31 )
    e32 = e3prov2*e3_gt_bound + (e3prov2*cos_0_1_deg + lambdahat_cross_e3prov2*sin_0_1_deg + lambdahat2*lambdahat_dot_e3prov*(1 - cos_0_1_deg))*e3_leq_bound
    print('e32  = %.16e'% e32 )
    e33 = e3prov3*e3_gt_bound + (e3prov3*cos_0_1_deg + lambdahat_cross_e3prov3*sin_0_1_deg + lambdahat3*lambdahat_dot_e3prov*(1 - cos_0_1_deg))*e3_leq_bound
    print('e33  = %.16e'% e33 )
    xi3=e31*n2-e32*n1
    print('xi3 = %.16e'% xi3)
    xi2 = -e31*n3 + e33*n1
    print('xi2  = %.16e'% xi2 )
    xi1 = e32*n3 - e33*n2
    print('xi1  = %.16e'% xi1 )
    v3=n1*xi2-n2*xi1
    print('v3 = %.16e'% v3)
    v2 = n3*xi1 - n1*xi3
    print('v2  = %.16e'% v2 )
    v1 = n2*xi3 - n3*xi2
    print('v1  = %.16e'% v1 )
    costheta=e31*n1+e32*n2+e33*n3
    print('costheta = %.16e'% costheta)
    sin2theta=1-costheta*costheta
    print('sin2theta = %.16e'% sin2theta)
    xisq = sin2theta
    print('xisq  = %.16e'% xisq )
    w2=a*a+r*r
    print('w2 = %.16e'% w2)
    Sigma=r*r+a*a*costheta*costheta
    print('Sigma = %.16e'% Sigma)
    Dinv=1+np.log(1+6*eta*u*u+2*(26-3*eta)*eta*u*u*u)
    print('Dinv = %.16e'% Dinv)
    Dinvprime = -u*u*(12*eta*u + 6*(26 - 3*eta)*eta*u*u)/(1 + 6*eta*u*u + 2*(26 - 3*eta)*eta*u*u*u)
    print('Dinvprime  = %.16e'% Dinvprime )
    omegatilde=2*a*r
    print('omegatilde = %.16e'% omegatilde)
    chi=(Skerr1*Lhat1+Skerr2*Lhat2+Skerr3*Lhat3)/(1-2*eta)+np.divide(1,2)*(Sperp1*Skerr1+Sperp2*Skerr2+Sperp3*Skerr3)/(Skerrmag*(1.-2.*eta))
    print('chi = %.16e'% chi)
    Kchi0=267.788*eta*eta*eta-126.687*eta*eta+10.2573*eta+1.7336
    print('Kchi0 = %.16e'% Kchi0)
    K=-59.1658*chi*chi*chi*eta*eta*eta-0.426958*chi*chi*chi*eta+1.43659*chi*chi*chi+31.1746*chi*chi*eta*eta*eta+6.16466*chi*chi*eta*eta-1.38086*chi*chi-27.5201*chi*eta*eta*eta+17.3736*chi*eta*eta+2.26831*chi*eta-1.62045*chi+Kchi0
    print('K = %.16e'% K)
    etaKminus1 = eta*K - 1
    print('etaKminus1  = %.16e'% etaKminus1 )
    Delta0=K*(eta*K-2)
    print('Delta0 = %.16e'% Delta0)
    Delta1=-2*etaKminus1*(K+Delta0)
    print('Delta1 = %.16e'% Delta1)
    Delta2=np.divide(1,2)*Delta1*(Delta1-4*etaKminus1)-a*a*etaKminus1*etaKminus1*Delta0
    print('Delta2 = %.16e'% Delta2)
    Delta3=-np.divide(1,3)*Delta1*Delta1*Delta1+etaKminus1*Delta1*Delta1+Delta2*Delta1-2*etaKminus1*(Delta2-etaKminus1)-a*a*etaKminus1*etaKminus1*Delta1
    print('Delta3 = %.16e'% Delta3)
    Delta4=np.divide(1,12)*(6*a*a*(Delta1*Delta1-2*Delta2)*etaKminus1*etaKminus1+3*Delta1*Delta1*Delta1*Delta1-8*etaKminus1*Delta1*Delta1*Delta1-12*Delta2*Delta1*Delta1+12*(2*etaKminus1*Delta2+Delta3)*Delta1+12*(np.divide(94,3)-np.divide(41,32)*np.pi*np.pi)*etaKminus1*etaKminus1+6*(Delta2*Delta2-4*Delta3*etaKminus1))
    print('Delta4 = %.16e'% Delta4)
    Delta5=etaKminus1*etaKminus1*(np.divide(-4237,60)+np.divide(128,5)*EMgamma+np.divide(2275,512)*np.pi*np.pi-np.divide(1,3)*a*a*(Delta1*Delta1*Delta1-3*Delta1*Delta2+3*Delta3)-(Delta1*Delta1*Delta1*Delta1*Delta1-5*Delta1*Delta1*Delta1*Delta2+5*Delta1*Delta2*Delta2+5*Delta1*Delta1*Delta3-5*Delta2*Delta3-5*Delta1*Delta4)/(5*etaKminus1*etaKminus1)+(Delta1*Delta1*Delta1*Delta1-4*Delta1*Delta1*Delta2+2*Delta2*Delta2+4*Delta1*Delta3-4*Delta4)/(2*etaKminus1)+np.divide(256,5)*np.log(2)+(np.divide(41,32)*np.pi*np.pi-np.divide(221,6))*eta)
    print('Delta5 = %.16e'% Delta5)
    Delta5l=etaKminus1*etaKminus1*np.divide(64,5)
    print('Delta5l = %.16e'% Delta5l)
    logarg=u*(Delta1+u*(Delta2+u*(Delta3+u*(Delta4+u*(Delta5+Delta5l*np.log(u))))))
    print('logarg = %.16e'% logarg)
    Deltaucalib = 1 + eta*(Delta0 + np.log(np.abs(1 + logarg)))
    print('Deltaucalib  = %.16e'% Deltaucalib )
    Deltaucalibprime = -eta*u*u*(Delta1+u*(2*Delta2+u*(3*Delta3+u*(4*Delta4+u*(5*(Delta5+Delta5l*np.log(u)))))))/(1+logarg)
    print('Deltaucalibprime  = %.16e'% Deltaucalibprime )
    Deltaubar=a*a*u*u+(2*u+1/etaKminus1)/etaKminus1
    print('Deltaubar = %.16e'% Deltaubar)
    Deltaubarprime = -2*a*a*u*u*u - 2*u*u/(etaKminus1)
    print('Deltaubarprime  = %.16e'% Deltaubarprime )
    Deltau=np.abs(Deltaubar*Deltaucalib)
    print('Deltau = %.16e'% Deltau)
    Deltauprime = Deltaubarprime*Deltaucalib + Deltaubar*Deltaucalibprime
    print('Deltauprime  = %.16e'% Deltauprime )
    Deltatprime=2*r*Deltau+r*r*Deltauprime
    print('Deltatprime = %.16e'% Deltatprime)
    Deltat=r*r*Deltau
    print('Deltat = %.16e'% Deltat)
    Deltar=Deltat*Dinv
    print('Deltar = %.16e'% Deltar)
    Deltarprime = Deltatprime*Dinv + Deltat*Dinvprime
    print('Deltarprime  = %.16e'% Deltarprime )
    Lambdat=np.abs(w2*w2-a*a*Deltat*sin2theta)
    print('Lambdat = %.16e'% Lambdat)
    csi=np.sqrt(np.abs(Deltar*Deltat))/w2
    print('csi = %.16e'% csi)
    csiprime = (Deltatprime*Deltar + Deltarprime*Deltat)/(2*np.sqrt(Deltar*Deltat)*w2) - 2.*r*np.sqrt(Deltat*Deltar)/(w2*w2)
    print('csiprime  = %.16e'% csiprime )
    csi1=1+(1-np.abs(1-tortoise))*(csi-1)
    print('csi1 = %.16e'% csi1)
    csi2=1+(np.divide(1,2)-np.divide(1,2)*np.sign(np.divide(3,2)-tortoise))*(csi-1)
    print('csi2 = %.16e'% csi2)
    prT=csi2*(p1*n1+p2*n2+p3*n3)
    print('prT = %.16e'% prT)
    phat3=p3-prT*(1-1/csi1)*n3
    print('phat3 = %.16e'% phat3)
    phat2 = p2 - prT*(1 - 1/csi1)*n2
    print('phat2  = %.16e'% phat2 )
    phat1 = p1 - prT*(1 - 1/csi1)*n1
    print('phat1  = %.16e'% phat1 )
    pdotxir=(phat1*xi1+phat2*xi2+phat3*xi3)*r
    print('pdotxir = %.16e'% pdotxir)
    pdotn=phat1*n1+phat2*n2+phat3*n3
    print('pdotn = %.16e'% pdotn)
    pdotvr=(phat1*v1+phat2*v2+phat3*v3)*r
    print('pdotvr = %.16e'% pdotvr)
    pphi=pdotxir
    print('pphi = %.16e'% pphi)
    Qcoeff2=1/(Sigma*sin2theta)
    print('Qcoeff2 = %.16e'% Qcoeff2)
    Qcoeff1=Sigma/(Lambdat*sin2theta)
    print('Qcoeff1 = %.16e'% Qcoeff1)
    DrSipn2=Deltar*pdotn*pdotn/Sigma
    print('DrSipn2 = %.16e'% DrSipn2)
    Q=1+DrSipn2+Qcoeff1*pdotxir*pdotxir+Qcoeff2*pdotvr*pdotvr
    print('Q = %.16e'% Q)
    Qminus1 = Q - 1
    print('Qminus1  = %.16e'% Qminus1 )
    Jtilde=np.sqrt(Deltar)
    print('Jtilde = %.16e'% Jtilde)
    exp2mu=Sigma
    print('exp2mu = %.16e'% exp2mu)
    expmu = np.sqrt(exp2mu)
    print('expmu  = %.16e'% expmu )
    Brtilde=(np.sqrt(Deltar)*Deltatprime-2*Deltat)/(2*np.sqrt(Deltar*Deltat))
    print('Brtilde = %.16e'% Brtilde)
    Btilde=np.sqrt(Deltat)
    print('Btilde = %.16e'% Btilde)
    exp2nu=Deltat*Sigma/Lambdat
    print('exp2nu = %.16e'% exp2nu)
    expnu = np.sqrt(exp2nu)
    print('expnu  = %.16e'% expnu )
    omega=omegatilde/Lambdat
    print('omega = %.16e'% omega)
    omegatildeprime=2*a
    print('omegatildeprime = %.16e'% omegatildeprime)
    Lambdatprime=4*(a*a+r*r)*r-a*a*Deltatprime*sin2theta
    print('Lambdatprime = %.16e'% Lambdatprime)
    mucostheta=a*a*costheta/Sigma
    print('mucostheta = %.16e'% mucostheta)
    nucostheta=a*a*w2*costheta*(w2-Deltat)/(Lambdat*Sigma)
    print('nucostheta = %.16e'% nucostheta)
    omegacostheta=-2*a*a*costheta*Deltat*omegatilde/(Lambdat*Lambdat)
    print('omegacostheta = %.16e'% omegacostheta)
    mur=r/Sigma-1/np.sqrt(Deltar)
    print('mur = %.16e'% mur)
    nur=r/Sigma+w2*(w2*Deltatprime-4*r*Deltat)/(2*Lambdat*Deltat)
    print('nur = %.16e'% nur)
    omegar=(Lambdat*omegatildeprime-Lambdatprime*omegatilde)/(Lambdat*Lambdat)
    print('omegar = %.16e'% omegar)
    dSO=147.481*chi*chi*chi*eta*eta-568.651*chi*chi*chi*eta+66.1987*chi*chi*chi-343.313*chi*chi*eta+2495.29*chi*eta*eta-44.5324
    print('dSO = %.16e'% dSO)
    sigmacoeffTerm3 = eta*dSO*u*u*u
    print('sigmacoeffTerm3  = %.16e'% sigmacoeffTerm3 )
    sigmacoeffTerm2=(-56.0/9.0*u*u+(-2.0/3.0*DrSipn2*u*u-109.0/36.0*Qminus1*u*u+(DrSipn2*Qminus1*u*u/4.0-5.0/16.0*Qminus1*Qminus1*u*u)*r)*r+(-7.0/3.0*u*u+(-49.0/8.0*DrSipn2*u*u+17.0/12.0*Qminus1*u*u+(45.0/8.0*DrSipn2*DrSipn2*u*u-13.0/8.0*DrSipn2*Qminus1*u*u)*r)*r)*eta)*eta
    print('sigmacoeffTerm2 = %.16e'% sigmacoeffTerm2)
    sigmacoeffTerm1=eta/12*(-8/r+3*Qminus1-36*DrSipn2)
    print('sigmacoeffTerm1 = %.16e'% sigmacoeffTerm1)
    sigmacoeff=sigmacoeffTerm1+sigmacoeffTerm2+sigmacoeffTerm3
    print('sigmacoeff = %.16e'% sigmacoeff)
    sigmastarcoeffTerm2=eta/(72*r*r)*(706+r*(-206*Qminus1+282*DrSipn2+r*Qminus1*(96*DrSipn2-23*Qminus1))+eta*(-54+r*(120*Qminus1-324*DrSipn2+r*(360*DrSipn2*DrSipn2+Qminus1*(-126*DrSipn2-3*Qminus1)))))
    print('sigmastarcoeffTerm2 = %.16e'% sigmastarcoeffTerm2)
    sigmastarcoeffTerm1=eta/12*(14/r+4*Qminus1-30*DrSipn2)
    print('sigmastarcoeffTerm1 = %.16e'% sigmastarcoeffTerm1)
    sigmastarcoeff=sigmastarcoeffTerm1+sigmastarcoeffTerm2
    print('sigmastarcoeff = %.16e'% sigmastarcoeff)
    Deltasigmastar3=sigmastar3*sigmastarcoeff+sigma3*sigmacoeff
    print('Deltasigmastar3 = %.16e'% Deltasigmastar3)
    Deltasigmastar2 = sigmastar2*sigmastarcoeff + sigma2*sigmacoeff
    print('Deltasigmastar2  = %.16e'% Deltasigmastar2 )
    Deltasigmastar1 = sigmastar1*sigmastarcoeff + sigma1*sigmacoeff
    print('Deltasigmastar1  = %.16e'% Deltasigmastar1 )
    Sstar3=sigmastar3+Deltasigmastar3
    print('Sstar3 = %.16e'% Sstar3)
    Sstar2 = sigmastar2 + Deltasigmastar2
    print('Sstar2  = %.16e'% Sstar2 )
    Sstar1 = sigmastar1 + Deltasigmastar1
    print('Sstar1  = %.16e'% Sstar1 )
    S3 = Sstar3
    print('S3  = %.16e'% S3 )
    S2 = Sstar2
    print('S2  = %.16e'% S2 )
    S1 = Sstar1
    print('S1  = %.16e'% S1 )
    Sstardotn=Sstar1*n1+Sstar2*n2+Sstar3*n3
    print('Sstardotn = %.16e'% Sstardotn)
    Sdote3=S1*e31+S2*e32+S3*e33
    print('Sdote3 = %.16e'% Sdote3)
    Sdotn=S1*n1+S2*n2+S3*n3
    print('Sdotn = %.16e'% Sdotn)
    Sdotv=S1*v1+S2*v2+S3*v3
    print('Sdotv = %.16e'% Sdotv)
    Sdotxi=S1*xi1+S2*xi2+S3*xi3
    print('Sdotxi = %.16e'% Sdotxi)
    HdsumTerm2=3*Sstardotn*Sstardotn
    print('HdsumTerm2 = %.16e'% HdsumTerm2)
    HdsumTerm1=Sstar1*Sstar1+Sstar2*Sstar2+Sstar3*Sstar3
    print('HdsumTerm1 = %.16e'% HdsumTerm1)
    Hdsum=HdsumTerm1-HdsumTerm2
    print('Hdsum = %.16e'% Hdsum)
    Hdcoeff=np.divide(1,2)/(r*r*r)
    print('Hdcoeff = %.16e'% Hdcoeff)
    Q4=2*prT*prT*prT*prT*u*u*(4-3*eta)*eta
    print('Q4 = %.16e'% Q4)
    gammappsum=Deltar/Sigma*pdotn*pdotn+1/Sigma*pdotvr*pdotvr/sin2theta+Sigma/Lambdat/sin2theta*pdotxir*pdotxir
    print('gammappsum = %.16e'% gammappsum)
    Hnsradicand=1+gammappsum+Q4
    print('Hnsradicand = %.16e'% Hnsradicand)
    alpha=np.sqrt(Deltat*Sigma/Lambdat)
    print('alpha = %.16e'% alpha)
    betapsum=omegatilde*pphi/Lambdat
    print('betapsum = %.16e'% betapsum)
    HssTerm3=expmu*expnu*pdotxir*(Jtilde*pdotn*Sdotxi*Btilde-expmu*expnu*pdotxir*Sdotn)+(pdotvr*(Sdotn*pdotvr-Jtilde*pdotn*Sdotv)-exp2mu*(np.sqrt(Q)+Q)*Sdotn*xisq)*Btilde*Btilde
    print('HssTerm3 = %.16e'% HssTerm3)
    HssTerm3coeff=omegacostheta/(2*exp2mu*expmu*expnu*Btilde*(Q+np.sqrt(Q)))
    print('HssTerm3coeff = %.16e'% HssTerm3coeff)
    HssTerm2=expmu*pdotxir*(expmu*exp2nu*pdotxir*Sdotv-expnu*pdotvr*Sdotxi*Btilde)+xisq*Btilde*Btilde*(exp2mu*(np.sqrt(Q)+Q)*Sdotv+Jtilde*pdotn*(pdotvr*Sdotn-Jtilde*pdotn*Sdotv))
    print('HssTerm2 = %.16e'% HssTerm2)
    HssTerm2coeff=Jtilde*omegar/(2*exp2mu*expmu*expnu*Btilde*(Q+np.sqrt(Q))*xisq)
    print('HssTerm2coeff = %.16e'% HssTerm2coeff)
    HssTerm1=omega*Sdote3
    print('HssTerm1 = %.16e'% HssTerm1)
    Hss=HssTerm1+HssTerm2coeff*HssTerm2+HssTerm3coeff*HssTerm3
    print('Hss = %.16e'% Hss)
    HsoTerm2c=Jtilde*Brtilde*expmu*expnu*pdotxir*(np.sqrt(Q)+1)*Sdotv
    print('HsoTerm2c = %.16e'% HsoTerm2c)
    HsoTerm2b=expmu*expnu*pdotxir*(2*np.sqrt(Q)+1)*(Jtilde*nur*Sdotv-nucostheta*Sdotn*xisq)*Btilde
    print('HsoTerm2b = %.16e'% HsoTerm2b)
    HsoTerm2a=Sdotxi*Jtilde*(mur*pdotvr*(np.sqrt(Q)+1)-mucostheta*pdotn*xisq-np.sqrt(Q)*(nur*pdotvr+(mucostheta-nucostheta)*pdotn*xisq))*Btilde*Btilde
    print('HsoTerm2a = %.16e'% HsoTerm2a)
    HsoTerm2=HsoTerm2a+HsoTerm2b-HsoTerm2c
    print('HsoTerm2 = %.16e'% HsoTerm2)
    HsoTerm2coeff=expnu/(exp2mu*Deltat*(Q+np.sqrt(Q))*xisq)
    print('HsoTerm2coeff = %.16e'% HsoTerm2coeff)
    HsoTerm1=exp2nu*(expmu*expnu-Btilde)*pdotxir*Sdote3/(expmu*Deltat*np.sqrt(Q)*xisq)
    print('HsoTerm1 = %.16e'% HsoTerm1)
    Hso=HsoTerm1+HsoTerm2coeff*HsoTerm2
    print('Hso = %.16e'% Hso)
    Hd=Hdcoeff*Hdsum
    print('Hd = %.16e'% Hd)
    Hns=betapsum+alpha*np.sqrt(Hnsradicand)
    print('Hns = %.16e'% Hns)
    Hs=Hso+Hss
    print('Hs = %.16e'% Hs)
    dSS=528.511*chi*chi*chi*eta*eta-41.0003*chi*chi*chi*eta+1161.78*chi*chi*eta*eta*eta-326.325*chi*chi*eta*eta+37.1964*chi*eta+706.958*eta*eta*eta-36.0272*eta+6.06807
    print('dSS = %.16e'% dSS)
    Heff = Hs + Hns - Hd + dSS*eta*u*u*u*u*(S1x*S1x + S1y*S1y + S1z*S1z + S2x*S2x + S2y*S2y + S2z*S2z)
    print('Heff  = %.16e'% Heff )
    Hreal=np.sqrt(1+2*eta*(np.abs(Heff)-1))
    print('Hreal = %.16e'% Hreal)
    return Hreal
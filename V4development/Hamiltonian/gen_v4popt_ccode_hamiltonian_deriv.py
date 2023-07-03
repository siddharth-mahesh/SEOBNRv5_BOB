import sympy as sp
from outputC import outputC, nrpyAbs
def gen_v4popt_ham_deriv(varind):
    m1,m2,tortoise = sp.symbols('m1 m2 tortoise',real=True)
    x,y,z,p1,p2,p3,S1x,S1y,S1z,S2x,S2y,S2z = sp.symbols("x y z p1 p2 p3 S1x S1y S1z S2x S2y S2z", real=True)
    vartodiff = [x,y,z,p1,p2,p3,S1x,S1y,S1z,S2x,S2y,S2z][varind]
    EMgamma = 0.577215664901532860606512090082402431
    M=m1+m2
    mu=m1*m2/M
    eta=mu/M
    r=sp.sqrt(x*x+y*y+z*z)
    u=1/r
    sigmastar3=m2/m1*S1z+m1/m2*S2z
    sigmastar2 = m2/m1*S1y + m1/m2*S2y
    sigmastar1 = m2/m1*S1x + m1/m2*S2x
    sigma3=S1z+S2z
    sigma2 = S1y + S2y
    sigma1 = S1x + S2x
    Skerr3=sigma3
    Skerr2 = sigma2
    Skerr1 = sigma1
    Skerrmag=sp.sqrt(Skerr1*Skerr1+Skerr2*Skerr2+Skerr3*Skerr3)
    Skerrhat3=Skerr3/Skerrmag
    Skerrhat2 = Skerr2/Skerrmag
    Skerrhat1 = Skerr1/Skerrmag
    a=Skerrmag
    L3=x*p2-y*p1
    L2 = z*p1 - x*p3
    L1 = y*p3 - z*p2
    Lnorm=sp.sqrt(L1*L1+L2*L2+L3*L3)
    Lhat3=L3/Lnorm
    Lhat2 = L2/Lnorm
    Lhat1 = L1/Lnorm
    S2dotLhat=S2x*Lhat1+S2y*Lhat2+S2z*Lhat3
    S1dotLhat = S1x*Lhat1 + S1y*Lhat2 + S1z*Lhat3
    S1perp3=S1z-S1dotLhat*Lhat3
    S1perp2 = S1y - S1dotLhat*Lhat2
    S1perp1 = S1x - S1dotLhat*Lhat1
    S2perp3=S2z-S2dotLhat*Lhat3
    S2perp2 = S2y - S2dotLhat*Lhat2
    S2perp1 = S2x - S2dotLhat*Lhat1
    Sperp3=S1perp3+S2perp3
    Sperp2 = S1perp2 + S2perp2
    Sperp1 = S1perp1 + S2perp1
    n3=z/r
    n2 = y/r
    n1 = x/r
    lambdavec3=Lhat1*n2-Lhat2*n3
    lambdavec2 = Lhat3*n1 - Lhat1*n3
    lambdavec1 = Lhat2*n3 - Lhat3*n2
    lambdavecnorm=sp.sqrt(lambdavec1*lambdavec1+lambdavec2*lambdavec2+lambdavec3*lambdavec3)
    lambdahat3=lambdavec3/lambdavecnorm
    lambdahat2 = lambdavec2/lambdavecnorm
    lambdahat1 = lambdavec1/lambdavecnorm
    lambdahat_dot_Skerrhat=lambdahat1*Skerrhat1+lambdahat2*Skerrhat2+lambdahat3*Skerrhat3
    lambdahat_cross_Skerrhat3=lambdahat1*Skerrhat2-lambdahat2*Skerrhat1
    lambdahat_cross_Skerrhat2 = lambdahat3*Skerrhat1 - lambdahat1*Skerrhat3
    lambdahat_cross_Skerrhat1 = lambdahat2*Skerrhat3 - lambdahat3*Skerrhat2
    Skerrhat_dot_n=Skerrhat1*n1+Skerrhat2*n2+Skerrhat3*n3
    cos_0_1_deg=0.9999983800004374
    sin_0_1_deg = 0.0017999990280001574
    TINYDOUBLE = 1e-100
    condition_lhs = 1 - nrpyAbs(Skerrhat_dot_n)
    condition_rhs = 1e-8
    greaterthan_bound = sp.Rational(1,2)*(condition_lhs - condition_rhs + nrpyAbs(condition_lhs - condition_rhs))/(condition_lhs - condition_rhs - TINYDOUBLE)
    lesserthan_equal_bound = sp.Rational(1,2)*(condition_lhs - condition_rhs - TINYDOUBLE - nrpyAbs(condition_lhs - condition_rhs - TINYDOUBLE))/(condition_lhs - condition_rhs - TINYDOUBLE)
    e31 = Skerrhat1*greaterthan_bound + (Skerrhat1*cos_0_1_deg + lambdahat_cross_Skerrhat1*sin_0_1_deg + lambdahat1*lambdahat_dot_Skerrhat*(1 - cos_0_1_deg))*lesserthan_equal_bound
    e32 = Skerrhat2*greaterthan_bound + (Skerrhat2*cos_0_1_deg + lambdahat_cross_Skerrhat2*sin_0_1_deg + lambdahat2*lambdahat_dot_Skerrhat*(1 - cos_0_1_deg))*lesserthan_equal_bound
    e33 = Skerrhat3*greaterthan_bound + (Skerrhat3*cos_0_1_deg + lambdahat_cross_Skerrhat3*sin_0_1_deg + lambdahat3*lambdahat_dot_Skerrhat*(1 - cos_0_1_deg))*lesserthan_equal_bound
    xi3=e31*n2-e32*n1
    xi2 = -e31*n3 + e33*n1
    xi1 = e32*n3 - e33*n2
    v3=n1*xi2-n2*xi1
    v2 = n3*xi1 - n1*xi3
    v1 = n2*xi3 - n3*xi2
    costheta=e31*n1+e32*n2+e33*n3
    sin2theta=1-costheta*costheta
    xisq = sin2theta
    w2=a*a+r*r
    Sigma=r*r+a*a*costheta*costheta
    Dinv=1+sp.log(1+6*eta*u*u+2*(26-3*eta)*eta*u*u*u)
    Dinvprime = -u*u*(12*eta*u + 6*(26 - 3*eta)*eta*u*u)/(1 + 6*eta*u*u + 2*(26 - 3*eta)*eta*u*u*u)
    omegatilde=2*a*r
    chi=(Skerr1*Lhat1+Skerr2*Lhat2+Skerr3*Lhat3)/(1-2*eta)+sp.Rational(1,2)*(Sperp1*Skerr1+Sperp2*Skerr2+Sperp3*Skerr3)/(Skerrmag*(1.-2.*eta))
    Kchi0=267.788*eta*eta*eta-126.687*eta*eta+10.2573*eta+1.7336
    K=-59.1658*chi*chi*chi*eta*eta*eta-0.426958*chi*chi*chi*eta+1.43659*chi*chi*chi+31.1746*chi*chi*eta*eta*eta+6.16466*chi*chi*eta*eta-1.38086*chi*chi-27.5201*chi*eta*eta*eta+17.3736*chi*eta*eta+2.26831*chi*eta-1.62045*chi+Kchi0
    etaKminus1 = eta*K - 1
    Delta0=K*(eta*K-2)
    Delta1=-2*etaKminus1*(K+Delta0)
    Delta2=sp.Rational(1,2)*Delta1*(Delta1-4*etaKminus1)-a*a*etaKminus1*etaKminus1*Delta0
    Delta3=-sp.Rational(1,3)*Delta1*Delta1*Delta1+etaKminus1*Delta1*Delta1+Delta2*Delta1-2*etaKminus1*(Delta2-etaKminus1)-a*a*etaKminus1*etaKminus1*Delta1
    Delta4=sp.Rational(1,12)*(6*a*a*(Delta1*Delta1-2*Delta2)*etaKminus1*etaKminus1+3*Delta1*Delta1*Delta1*Delta1-8*etaKminus1*Delta1*Delta1*Delta1-12*Delta2*Delta1*Delta1+12*(2*etaKminus1*Delta2+Delta3)*Delta1+12*(sp.Rational(94,3)-sp.Rational(41,32)*sp.pi*sp.pi)*etaKminus1*etaKminus1+6*(Delta2*Delta2-4*Delta3*etaKminus1))
    Delta5=etaKminus1*etaKminus1*(sp.Rational(-4237,60)+sp.Rational(128,5)*EMgamma+sp.Rational(2275,512)*sp.pi*sp.pi-sp.Rational(1,3)*a*a*(Delta1*Delta1*Delta1-3*Delta1*Delta2+3*Delta3)-(Delta1*Delta1*Delta1*Delta1*Delta1-5*Delta1*Delta1*Delta1*Delta2+5*Delta1*Delta2*Delta2+5*Delta1*Delta1*Delta3-5*Delta2*Delta3-5*Delta1*Delta4)/(5*etaKminus1*etaKminus1)+(Delta1*Delta1*Delta1*Delta1-4*Delta1*Delta1*Delta2+2*Delta2*Delta2+4*Delta1*Delta3-4*Delta4)/(2*etaKminus1)+sp.Rational(256,5)*sp.log(2)+(sp.Rational(41,32)*sp.pi*sp.pi-sp.Rational(221,6))*eta)
    Delta5l=etaKminus1*etaKminus1*sp.Rational(64,5)
    logarg=u*(Delta1+u*(Delta2+u*(Delta3+u*(Delta4+u*(Delta5+Delta5l*sp.log(u))))))
    Deltaucalib = 1 + eta*(Delta0 + sp.log(nrpyAbs(1 + logarg)))
    Deltaucalibprime = -eta*u*u*(Delta1+u*(2*Delta2+u*(3*Delta3+u*(4*Delta4+u*(5*(Delta5+Delta5l*sp.log(u)))))))/(1+logarg)
    Deltaubar=a*a*u*u+(2*u+1/etaKminus1)/etaKminus1
    Deltaubarprime = -2*a*a*u*u*u - 2*u*u/(etaKminus1)
    Deltau=nrpyAbs(Deltaubar*Deltaucalib)
    Deltauprime = Deltaubarprime*Deltaucalib + Deltaubar*Deltaucalibprime
    Deltatprime=2*r*Deltau+r*r*Deltauprime
    Deltat=r*r*Deltau
    Deltar=Deltat*Dinv
    Hreal = Deltar
    Deltarprime = Deltatprime*Dinv + Deltat*Dinvprime
    Lambdat=nrpyAbs(w2*w2-a*a*Deltat*sin2theta)
    csi=sp.sqrt(nrpyAbs(Deltar*Deltat))/w2
    csiprime = (Deltatprime*Deltar + Deltarprime*Deltat)/(2*sp.sqrt(Deltar*Deltat)*w2) - 2.*r*sp.sqrt(Deltat*Deltar)/(w2*w2)
    csi1=1+(1-nrpyAbs(1-tortoise))*(csi-1)
    csi2=1+(sp.Rational(1,2)-sp.Rational(1,2)*sp.sign(sp.Rational(3,2)-tortoise))*(csi-1)
    prT=csi2*(p1*n1+p2*n2+p3*n3)
    phat3=p3-prT*(1-1/csi1)*n3
    phat2 = p2 - prT*(1 - 1/csi1)*n2
    phat1 = p1 - prT*(1 - 1/csi1)*n1
    pdotxir=(phat1*xi1+phat2*xi2+phat3*xi3)*r
    pdotn=phat1*n1+phat2*n2+phat3*n3
    pdotvr=(phat1*v1+phat2*v2+phat3*v3)*r
    pphi=pdotxir
    Qcoeff2=1/(Sigma*sin2theta)
    Qcoeff1=Sigma/(Lambdat*sin2theta)
    DrSipn2=Deltar*pdotn*pdotn/Sigma
    Q=1+DrSipn2+Qcoeff1*pdotxir*pdotxir+Qcoeff2*pdotvr*pdotvr
    Qminus1 = Q - 1
    Jtilde=sp.sqrt(Deltar)
    exp2mu=Sigma
    expmu = sp.sqrt(exp2mu)
    Brtilde=(sp.sqrt(Deltar)*Deltatprime-2*Deltat)/(2*sp.sqrt(Deltar*Deltat))
    Btilde=sp.sqrt(Deltat)
    exp2nu=Deltat*Sigma/Lambdat
    expnu = sp.sqrt(exp2nu)
    omega=omegatilde/Lambdat
    omegatildeprime=2*a
    Lambdatprime=4*(a*a+r*r)*r-a*a*Deltatprime*sin2theta
    mucostheta=a*a*costheta/Sigma
    nucostheta=a*a*w2*costheta*(w2-Deltat)/(Lambdat*Sigma)
    omegacostheta=-2*a*a*costheta*Deltat*omegatilde/(Lambdat*Lambdat)
    mur=r/Sigma-1/sp.sqrt(Deltar)
    nur=r/Sigma+w2*(w2*Deltatprime-4*r*Deltat)/(2*Lambdat*Deltat)
    omegar=(Lambdat*omegatildeprime-Lambdatprime*omegatilde)/(Lambdat*Lambdat)
    dSO=147.481*chi*chi*chi*eta*eta-568.651*chi*chi*chi*eta+66.1987*chi*chi*chi-343.313*chi*chi*eta+2495.29*chi*eta*eta-44.5324
    sigmacoeffTerm3 = eta*dSO*u*u*u
    sigmacoeffTerm2=eta/(144*r*r)*(-896+r*(-436*Qminus1-96*DrSipn2+r*(-45*Qminus1*Qminus1+36*Qminus1*DrSipn2))+eta*(-336+r*(204*Qminus1-882*DrSipn2+r*(810*DrSipn2*DrSipn2-234*Qminus1*DrSipn2))))
    sigmacoeffTerm1=eta/12*(-8/r+3*Qminus1-36*DrSipn2)
    sigmacoeff=sigmacoeffTerm1+sigmacoeffTerm2+sigmacoeffTerm3
    sigmastarcoeffTerm2=eta/(72*r*r)*(706+r*(-206*Qminus1+282*DrSipn2+r*Qminus1*(96*DrSipn2-23*Qminus1))+eta*(-54+r*(120*Qminus1-324*DrSipn2+r*(360*DrSipn2*DrSipn2+Qminus1*(-126*DrSipn2-3*Qminus1)))))
    sigmastarcoeffTerm1=eta/12*(14/r+4*Qminus1-30*DrSipn2)
    sigmastarcoeff=sigmastarcoeffTerm1+sigmastarcoeffTerm2
    Deltasigmastar3=sigmastar3*sigmastarcoeff+sigma3*sigmacoeff
    Deltasigmastar2 = sigmastar2*sigmastarcoeff + sigma2*sigmacoeff
    Deltasigmastar1 = sigmastar1*sigmastarcoeff + sigma1*sigmacoeff
    Sstar3=sigmastar3+Deltasigmastar3
    Sstar2 = sigmastar2 + Deltasigmastar2
    Sstar1 = sigmastar1 + Deltasigmastar1
    S3 = Sstar3
    S2 = Sstar2
    S1 = Sstar1
    Sstardotn=Sstar1*n1+Sstar2*n2+Sstar3*n3
    SdotSkerrhat=S1*Skerrhat1+S2*Skerrhat2+S3*Skerrhat3
    Sdotn=S1*n1+S2*n2+S3*n3
    Sdotv=S1*v1+S2*v2+S3*v3
    Sdotxi=S1*xi1+S2*xi2+S3*xi3
    HdsumTerm2=3*Sstardotn*Sstardotn
    HdsumTerm1=Sstar1*Sstar1+Sstar2*Sstar2+Sstar3*Sstar3
    Hdsum=HdsumTerm1-HdsumTerm2
    Hdcoeff=sp.Rational(1,2)/(r*r*r)
    Q4=2*prT*prT*prT*prT*u*u*(4-3*eta)*eta
    gammappsum=Deltar/Sigma*pdotn*pdotn+1/Sigma*pdotvr*pdotvr/sin2theta+Sigma/Lambdat/sin2theta*pdotxir*pdotxir
    Hnsradicand=1+gammappsum+Q4
    alpha=sp.sqrt(Deltat*Sigma/Lambdat)
    betapsum=omegatilde*pphi/Lambdat
    HssTerm3=expmu*expnu*pdotxir*(Jtilde*pdotn*Sdotxi*Btilde-expmu*expnu*pdotxir*Sdotn)+(pdotvr*(Sdotn*pdotvr-Jtilde*pdotn*Sdotv)-exp2mu*(sp.sqrt(Q)+Q)*Sdotn*xisq)*Btilde*Btilde
    HssTerm3coeff=omegacostheta/(2*exp2mu*expmu*expnu*Btilde*(Q+sp.sqrt(Q)))
    HssTerm2=expmu*pdotxir*(expmu*exp2nu*pdotxir*Sdotv-expnu*pdotvr*Sdotxi*Btilde)+xisq*Btilde*Btilde*(exp2mu*(sp.sqrt(Q)+Q)*Sdotv+Jtilde*pdotn*(pdotvr*Sdotn-Jtilde*pdotn*Sdotv))
    HssTerm2coeff=Jtilde*omegar/(2*exp2mu*expmu*expnu*Btilde*(Q+sp.sqrt(Q))*xisq)
    HssTerm1=omega*SdotSkerrhat
    Hss=HssTerm1+HssTerm2coeff*HssTerm2+HssTerm3coeff*HssTerm3
    HsoTerm2c=Jtilde*Brtilde*expmu*expnu*pdotxir*(sp.sqrt(Q)+1)*Sdotv
    HsoTerm2b=expmu*expnu*pdotxir*(2*sp.sqrt(Q)+1)*(Jtilde*nur*Sdotv-nucostheta*Sdotn*xisq)*Btilde
    HsoTerm2a=Sdotxi*Jtilde*(mur*pdotvr*(sp.sqrt(Q)+1)-mucostheta*pdotn*xisq-sp.sqrt(Q)*(nur*pdotvr+(mucostheta-nucostheta)*pdotn*xisq))*Btilde*Btilde
    HsoTerm2=HsoTerm2a+HsoTerm2b-HsoTerm2c
    HsoTerm2coeff=expnu/(exp2mu*Btilde*Btilde*(Q+sp.sqrt(Q))*xisq)
    HsoTerm1=exp2nu*(expmu*expnu-Btilde)*pdotxir*SdotSkerrhat/(expmu*Btilde*Btilde*sp.sqrt(Q)*xisq)
    Hso=HsoTerm1+HsoTerm2coeff*HsoTerm2
    Hd=Hdcoeff*Hdsum
    Hns=betapsum+alpha*sp.sqrt(Hnsradicand)
    Hs=Hso+Hss
    dSS=528.511*chi*chi*chi*eta*eta-41.0003*chi*chi*chi*eta+1161.78*chi*chi*eta*eta*eta-326.325*chi*chi*eta*eta+37.1964*chi*eta+706.958*eta*eta*eta-36.0272*eta+6.06807
    Heff = Hs + Hns - Hd + dSS*eta*u*u*u*u*(S1x*S1x + S1y*S1y + S1z*S1z + S2x*S2x + S2y*S2y + S2z*S2z)
    #Hreal=sp.sqrt(1+2*eta*(nrpyAbs(Heff)-1))
    dHreal = sp.diff(Hreal,vartodiff)
    return outputC(dHreal,"Hreal",filename="returnstring",params="CSE_sorting=none,includebraces=False,preindent=2,outCverbose=False,declareoutputvars=True,GoldenKernelsEnable=False")

gen_v4popt_ham_deriv(0)
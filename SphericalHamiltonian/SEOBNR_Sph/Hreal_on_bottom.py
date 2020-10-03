import numpy as np
def compute_Hreal(m1=23., m2=10., EMgamma=0.577215664901532860606512090082402431, tortoise=1, r = 21.29681018601393,theta = 1.5707963267948966,phi = 0.0,pr = 0.0,ptheta = 9.019549949487391e-21,pphi = 4.973638129861999,S1r = 0.004857667584940312, S1theta = 0.01457311842632286, S1phi = 0.009715161660389764, S2r = 0.003673094582185491,S2theta = -0.005509696538546906,S2phi =  -0.004591302628615413):
    M = m1 + m2
    mu=m1*m2/M
    eta=mu/M
    u=1/r
    sigmastar3=m2/m1*S1phi+m1/m2*S2phi
    sigmastar2 = m2/m1*S1theta + m1/m2*S2theta
    sigmastar1 = m2/m1*S1r + m1/m2*S2r
    ##sigmastar3=m2/m1*S1z+m1/m2*S2z
    ##sigmastar2 = m2/m1*S1y + m1/m2*S2y
    ##sigmastar1 = m2/m1*S1x + m1/m2*S2x
    sigma3=S1phi+S2phi##Sid:Spinsshouldbeinpolardirections
    sigma2 = S1theta + S2theta
    sigma1 = S1r + S2r
    ##sigma3=S1z+S2z
    ##sigma2 = S1y + S2y
    ##sigma1 = S1x + S2x
    Skerr3=sigma3##Sid:Spinsshouldbeinpolardirections
    Skerr2 = sigma2
    Skerr1 = sigma1
    Skerrmag=np.sqrt(Skerr1*Skerr1+Skerr2*Skerr2+Skerr3*Skerr3)
    Skerrhat3=Skerr3/Skerrmag
    Skerrhat2 = Skerr2/Skerrmag
    Skerrhat1 = Skerr1/Skerrmag
    a=Skerrmag
    n3=0
    n2 = 0
    n1 = 1
    ##n3=z/r
    ##n2 = y/r
    ##n1 = x/r
    e33=Skerrhat3##Sid:insphericalpolar,nisjust(1,0,0)
    e32 = Skerrhat2
    e31 = Skerrhat1
    xi3=e31*n2-e32*n1
    xi2 = e33*n1 - e31*n3
    xi1 = e32*n3 - e33*n2
    v3=n1*xi2-n2*xi1##Sid:Testingoutatypohere
    v2 = n3*xi1 - n1*xi3
    v1 = n2*xi3 - n3*xi2
    costhetaBL=e31*n1+e32*n2+e33*n3
    sin2thetaBL=1-costhetaBL*costhetaBL
    xisq = sin2thetaBL
    w2=a*a+r*r
    Sigma=r*r+a*a*costhetaBL*costhetaBL
    ##Sigma=r*r+a*a*costheta*costheta
    sintheta=np.sin(theta)##Sid:costhetaBLsubstitution
    Dinv=1+np.log(1+6*eta*u*u+2*(26-3*eta)*eta*u*u*u)
    omegatilde=2*a*r
    K=1.712-1.803949138004582*eta-39.77229225266885*eta*eta+103.16588921239249*eta*eta*eta
    etaKminus1 = eta*K - 1
    Delta0=K*(eta*K-2)
    Delta1 = -2*etaKminus1*(K + Delta0)
    Delta2 = np.divide(1,2)*Delta1*(Delta1 - 4*etaKminus1) - a*a*etaKminus1*etaKminus1*Delta0
    Delta3=-np.divide(1,3)*Delta1*Delta1*Delta1+etaKminus1*Delta1*Delta1+Delta2*Delta1-2*etaKminus1*(Delta2-etaKminus1)-a*a*etaKminus1*etaKminus1*Delta1
    Delta4=np.divide(1,12)*(6*a*a*(Delta1*Delta1-2*Delta2)*etaKminus1*etaKminus1+3*Delta1*Delta1*Delta1*Delta1-8*etaKminus1*Delta1*Delta1*Delta1-12*Delta2*Delta1*Delta1+12*(2*etaKminus1*Delta2+Delta3)*Delta1+12*(np.divide(94,3)-np.divide(41,32)*np.pi*np.pi)*etaKminus1*etaKminus1+6*(Delta2*Delta2-4*Delta3*etaKminus1))
    Delta5=etaKminus1*etaKminus1*((np.divide(-4237,60)+np.divide(128,5)*EMgamma+np.divide(2275,512)*np.pi*np.pi-np.divide(1,3)*a*a*(Delta1*Delta1*Delta1-3*Delta1*Delta2+3*Delta3)-(Delta1*Delta1*Delta1*Delta1*Delta1-5*Delta1*Delta1*Delta1*Delta2+5*Delta1*Delta2*Delta2+5*Delta1*Delta1*Delta3-5*Delta2*Delta3-5*Delta1*Delta4)/(5*etaKminus1*etaKminus1)+(Delta1*Delta1*Delta1*Delta1-4*Delta1*Delta1*Delta2+2*Delta2*Delta2+4*Delta1*Delta3-4*Delta4)/(2*etaKminus1)+np.divide(256,5)*np.log(2)))
    Delta5l = etaKminus1*etaKminus1*np.divide(64,5)
    logarg=u*(Delta1+u*(Delta2+u*(Delta3+u*(Delta4+u*(Delta5+Delta5l*np.log(u))))))
    Deltaucalib = 1 + eta*(Delta0 + np.log(1 + logarg))
    Deltaucalibprm=-eta*u*u*(Delta1+u*(2*Delta2+u*(3*Delta3+u*(4*Delta4+u*(5*(Delta5+Delta5l*np.log(u)))))))/(1+logarg)
    Deltaubar=a*a*u*u+(2*u+1/etaKminus1)/etaKminus1
    Deltaubarprm = -2*a*a*u*u*u - 2*u*u/(etaKminus1)
    Deltau=Deltaubar*Deltaucalib
    Deltauprm = Deltaubarprm*Deltaucalib + Deltaubar*Deltaucalibprm
    Deltatprm=2*r*Deltau+r*r*Deltauprm
    Deltat=r*r*Deltau
    Deltar=Deltat*Dinv
    Lambdat=w2*w2-a*a*Deltat*sin2thetaBL
    ##Lambdat=w2*w2-a*a*Deltat*sin2theta
    csi=np.sqrt(Deltar*Deltat)/w2##Sid:sin2thetaBLsubstitution
    csi1=1+(1-np.abs(1-tortoise))*(csi-1)
    csi2=1+(np.divide(1,2)-np.divide(1,2)*np.sign(np.divide(3,2)-tortoise))*(csi-1)
    prT=csi2*(pr*n1+ptheta*n2/r+pphi*n3/r/sintheta)
    ##prT=csi2*(p1*n1+p2*n2+p3*n3)
    phat3=pphi/r/sintheta+prT*(1-1/csi1)*n3##Sid:Rewritewithsphericalmomenta
    phat2 = ptheta/r + prT*(1 - 1/csi1)*n2
    phat1 = pr + prT*(1 - 1/csi1)*n1
    sintheta = np.sin(theta)
    ##phat3=p3+prT*(1-1/csi1)*n3
    ##phat2 = p2 + prT*(1 - 1/csi1)*n2
    ##phat1 = p1 + prT*(1 - 1/csi1)*n1
    pdotxir=pphi##Sid:Rewritewithsphericalmomenta
    #pdotxir = (phat1*xi1 + phat2*xi2 + phat3*xi3)*r
    pdotn=phat1*n1+phat2*n2+phat3*n3
    pdotvr=(phat1*v1+phat2*v2+phat3*v3)*r
    Qcoeff2=1/(Sigma*sin2thetaBL)
    ##Qcoeff2=1/(Sigma*sin2theta)
    Qcoeff1=Sigma/(Lambdat*sin2thetaBL)##Sid:sin2thetaBLsubstitution
    ##Qcoeff1=Sigma/(Lambdat*sin2theta)
    DrSipn2=Deltar*pdotn*pdotn/Sigma##Sid:sin2thetaBLsubstitution
    Q=1+DrSipn2+Qcoeff1*pdotxir*pdotxir+Qcoeff2*pdotvr*pdotvr
    Qminus1 = Q - 1
    Jtilde=np.sqrt(Deltar)
    exp2mu=Sigma
    expmu = np.sqrt(exp2mu)
    Brtilde=(np.sqrt(Deltar)*Deltatprm-2*Deltat)/(2*np.sqrt(Deltar*Deltat))
    Btilde=np.sqrt(Deltat)
    exp2nu=Deltat*Sigma/Lambdat
    expnu = np.sqrt(exp2nu)
    omega=omegatilde/Lambdat
    omegatildeprm=2*a
    Lambdatprm=4*(a*a+r*r)*r-2*a*a*Deltatprm*sin2thetaBL
    ##Lambdatprm=4*(a*a+r*r)*r-2*a*a*Deltatprm*sin2theta
    mucostheta=a*a*costhetaBL/Sigma##Sid:sin2thetaBLsubstitution
    ##mucostheta=a*a*costheta/Sigma
    nucostheta=a*a*w2*costhetaBL*(w2-Deltat)/(Lambdat*Sigma)##Sid:costhetaBLsubstitution
    ##nucostheta=a*a*w2*costheta*(w2-Deltat)/(Lambdat*Sigma)
    omegacostheta=-2*a*a*costhetaBL*Deltat*omegatilde/(Lambdat*Lambdat)##Sid:costhetaBLsubstitution
    ##omegacostheta=-2*a*a*costheta*Deltat*omegatilde/(Lambdat*Lambdat)
    mur=r/Sigma-1/np.sqrt(Deltar)##Sid:costhetaBLsubstitution
    nur=r/Sigma+w2*(w2*Deltatprm-4*r*Deltat)/(2*Lambdat*Deltat)
    omegar=(Lambdat*omegatildeprm-Lambdatprm*omegatilde)/(Lambdat*Lambdat)
    dSO=-74.71-156.*eta+627.5*eta*eta
    sigmacoeffTerm3 = eta*dSO*u*u*u
    sigmacoeffTerm2=eta/(144*r*r)*(-896+r*(-436*Qminus1-96*DrSipn2+r*(-45*Qminus1*Qminus1+36*Qminus1*DrSipn2))+eta*(-336+r*(204*Qminus1-882*DrSipn2+r*(810*DrSipn2*DrSipn2-234*Qminus1*DrSipn2))))
    sigmacoeffTerm1=eta/12*(-8*u+3*Qminus1-36*DrSipn2)
    ##sigmacoeffTerm1=eta/12*(-8/r+3*Qminus1-36*DrSipn2)
    sigmacoeff=sigmacoeffTerm1+sigmacoeffTerm2+sigmacoeffTerm3##Sid:Again,multiplybyuinstead
    sigmastarcoeffTerm2=eta/(72*r*r)*(706+r*(-206*Qminus1+282*DrSipn2+r*Qminus1*(96*DrSipn2-23*Qminus1))+eta*(-54+r*(120*Qminus1-324*DrSipn2+r*(360*DrSipn2*DrSipn2+Qminus1*(-126*DrSipn2-3*Qminus1)))))
    sigmastarcoeffTerm1=eta/12*(14*u+4*Qminus1-30*DrSipn2)
    ##sigmastarcoeffTerm1=eta/12*(14/r+4*Qminus1-30*DrSipn2)
    sigmastarcoeff=sigmastarcoeffTerm1+sigmastarcoeffTerm2##Sid:Multiplybyuoncemore
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
    Hdcoeff=np.divide(1,2)*(u*u*u)
    ##Hdcoeff=np.divide(1,2)/(r*r*r)
    ## Also, where did the eta go?-> I checked with LAL and we compute Hcap = H/mu so the mu in the numerator disappears and the whole thing is in M = 1units so there is no eta term.
    Q4=2*prT*prT*prT*prT*u*u*(4-3*eta)*eta##Sid:Whydividebyrwhenyoucanmultiplybyu?
    gammappsum=Deltar/Sigma*pdotn*pdotn+1/Sigma*pdotvr*pdotvr/sin2thetaBL+Sigma/Lambdat/sin2thetaBL*pdotxir*pdotxir
    ##gammappsum=Deltar/Sigma*pdotn*pdotn+1/Sigma*pdotvr*pdotvr/sin2theta+Sigma/Lambdat/sin2theta*pdotxir*pdotxir
    Hnsradicand=1+gammappsum+Q4##Sid:sin2thetaBLsubstitution
    alpha=np.sqrt(Deltat*Sigma/Lambdat)
    betapsum=omegatilde*pdotxir/Lambdat
    ##betapsum=omegatilde*pphi/Lambdat
    HssTerm3=expmu*expnu*pdotxir*(Jtilde*pdotn*Sdotxi*Btilde-expmu*expnu*pdotxir*Sdotn)+(pdotvr*(Sdotn*pdotvr-Jtilde*pdotn*Sdotv)-exp2mu*(np.sqrt(Q)+Q)*Sdotn*xisq)*Btilde*Btilde##Sid:Changedtopdotxirtoeliminatepphi
    HssTerm3coeff=omegacostheta/(2*exp2mu*expmu*expnu*Btilde*(Q+np.sqrt(Q)))
    HssTerm2=expmu*pdotxir*(expmu*exp2nu*pdotxir*Sdotv-expnu*pdotvr*Sdotxi*Btilde)+xisq*Btilde*Btilde*(exp2mu*(np.sqrt(Q)+Q)*Sdotv+Jtilde*pdotn*(pdotvr*Sdotn-Jtilde*pdotn*Sdotv))
    HssTerm2coeff=Jtilde*omegar/(2*exp2mu*expmu*expnu*Btilde*(Q+np.sqrt(Q))*xisq)
    HssTerm1=omega*SdotSkerrhat
    Hss=HssTerm1+HssTerm2coeff*HssTerm2+HssTerm3coeff*HssTerm3
    HsoTerm2c=Jtilde*Brtilde*expmu*expnu*pdotxir*(np.sqrt(Q)+1)*Sdotv
    HsoTerm2b=expmu*expnu*pdotxir*(2*np.sqrt(Q)+1)*(Jtilde*nur*Sdotv-nucostheta*Sdotn*xisq)*Btilde
    HsoTerm2a=Sdotxi*Jtilde*(mur*pdotvr*(np.sqrt(Q)+1)-mucostheta*pdotn*xisq-np.sqrt(Q)*(nur*pdotvr+(mucostheta-nucostheta)*pdotn*xisq))*Btilde*Btilde
    HsoTerm2=HsoTerm2a+HsoTerm2b-HsoTerm2c
    HsoTerm2coeff=expnu/(exp2mu*Btilde*Btilde*(Q+np.sqrt(Q))*xisq)
    HsoTerm1=exp2nu*(expmu*expnu-Btilde)*pdotxir*SdotSkerrhat/(expmu*Btilde*Btilde*np.sqrt(Q)*xisq)
    Hso=HsoTerm1+HsoTerm2coeff*HsoTerm2
    Hd=Hdcoeff*Hdsum
    Hns=betapsum+alpha*np.sqrt(Hnsradicand)
    Hs=Hso+Hss
    dSS=8.127-154.2*eta+830.8*eta*eta
    Heff = Hs + Hns - Hd + dSS*eta*u*u*u*u*(S1r*S1r + S1theta*S1theta + S1phi*S1phi + S2r*S2r + S2theta*S2theta + S2phi*S2phi)
    #dSS=8.127-154.2*eta+830.8*eta*eta
    #Heff = Hs + Hns - Hd + dSS*eta*u*u*u*u*(S1x*S1x + S1y*S1y + S1z*S1z + S2x*S2x + S2y*S2y + S2z*S2z)
    Hreal=np.sqrt(1+2*eta*(Heff-1))##Sid:Spinvectorsshouldbelistedinpolarcomponents
    return Hreal
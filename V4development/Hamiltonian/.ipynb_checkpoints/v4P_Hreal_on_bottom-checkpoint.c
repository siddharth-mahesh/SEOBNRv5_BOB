#include math.h 
DOUBLE compute_v4P_Hreal(DOUBLE m1,DOUBLE m2,INT tortoise,DOUBLE x,DOUBLE y,DOUBLE z,DOUBLE p1,DOUBLE p2,DOUBLE p3,DOUBLE S1x,DOUBLE S1y,DOUBLE S1z,DOUBLE S2x,DOUBLE S2y,DOUBLE S2z){
    DOUBLE EMgamma = 0.577215664901532860606512090082402431 ;
    DOUBLE M=m1+m2 ;
    DOUBLE mu=m1*m2/M ;
    DOUBLE eta=mu/M ;
    DOUBLE r=sqrt(x*x+y*y+z*z) ;
    DOUBLE u=1/r ;
    DOUBLE sigmastar3=m2/m1*S1z+m1/m2*S2z ;
    DOUBLE sigmastar2 = m2/m1*S1y + m1/m2*S2y ;
    DOUBLE sigmastar1 = m2/m1*S1x + m1/m2*S2x ;
    DOUBLE sigma3=S1z+S2z ;
    DOUBLE sigma2 = S1y + S2y ;
    DOUBLE sigma1 = S1x + S2x ;
    DOUBLE Skerr3=sigma3 ;
    DOUBLE Skerr2 = sigma2 ;
    DOUBLE Skerr1 = sigma1 ;
    DOUBLE Skerrmag=sqrt(Skerr1*Skerr1+Skerr2*Skerr2+Skerr3*Skerr3) ;
    DOUBLE Skerrhat3=Skerr3/Skerrmag ;
    DOUBLE Skerrhat2 = Skerr2/Skerrmag ;
    DOUBLE Skerrhat1 = Skerr1/Skerrmag ;
    DOUBLE a=Skerrmag ;
    DOUBLE L3=x*p2-y*p1 ;
    DOUBLE L2 = z*p1 - x*p3 ;
    DOUBLE L1 = y*p3 - z*p2 ;
    DOUBLE Lnorm=sqrt(L1*L1+L2*L2+L3*L3) ;
    DOUBLE Lhat3=L3/Lnorm ;
    DOUBLE Lhat2 = L2/Lnorm ;
    DOUBLE Lhat1 = L1/Lnorm ;
    DOUBLE S2dotLhat=S2x*Lhat1+S2y*Lhat2+S2z*Lhat3 ;
    DOUBLE S1dotLhat = S1x*Lhat1 + S1y*Lhat2 + S1z*Lhat3 ;
    DOUBLE S1perp3=S1z-S1dotLhat*Lhat3 ;
    DOUBLE S1perp2 = S1y - S1dotLhat*Lhat2 ;
    DOUBLE S1perp1 = S1x - S1dotLhat*Lhat1 ;
    DOUBLE S2perp3=S2z-S2dotLhat*Lhat3 ;
    DOUBLE S2perp2 = S2y - S2dotLhat*Lhat2 ;
    DOUBLE S2perp1 = S2x - S2dotLhat*Lhat1 ;
    DOUBLE Sperp3=S1perp3+S2perp3 ;
    DOUBLE Sperp2 = S1perp2 + S2perp2 ;
    DOUBLE Sperp1 = S1perp1 + S2perp1 ;
    DOUBLE n3=z/r ;
    DOUBLE n2 = y/r ;
    DOUBLE n1 = x/r ;
    DOUBLE e33=Skerrhat3 ;
    DOUBLE e32 = Skerrhat2 ;
    DOUBLE e31 = Skerrhat1 ;
    DOUBLE xi3=e31*n2-e32*n1 ;
    DOUBLE xi2 = -e31*n3 + e33*n1 ;
    DOUBLE xi1 = e32*n3 - e33*n2 ;
    DOUBLE v3=n1*xi2-n2*xi1 ;
    DOUBLE v2 = n3*xi1 - n1*xi3 ;
    DOUBLE v1 = n2*xi3 - n3*xi2 ;
    DOUBLE costheta=e31*n1+e32*n2+e33*n3 ;
    DOUBLE sin2theta=1-costheta*costheta ;
    DOUBLE xisq = sin2theta ;
    DOUBLE w2=a*a+r*r ;
    DOUBLE Sigma=r*r+a*a*costheta*costheta ;
    DOUBLE Dinv=1+log(1+6*eta*u*u+2*(26-3*eta)*eta*u*u*u) ;
    DOUBLE omegatilde=2*a*r ;
    DOUBLE chi=(Skerr1*Lhat1+Skerr2*Lhat2+Skerr3*Lhat3)/(1-2*eta)+np.divide(1,2)*(Sperp1*Skerr1+Sperp2*Skerr2+Sperp3*Skerr3)/(Skerrmag*(1.-2.*eta)) ;
    DOUBLE Kchi0=267.788*eta*eta*eta-126.687*eta*eta+10.2573*eta+1.7336 ;
    DOUBLE K=-59.1658*chi*chi*chi*eta*eta*eta-0.426958*chi*chi*chi*eta+1.43659*chi*chi*chi+31.1746*chi*chi*eta*eta*eta+6.16466*chi*chi*eta*eta-1.38086*chi*chi-27.5201*chi*eta*eta*eta+17.3736*chi*eta*eta+2.26831*chi*eta-1.62045*chi+Kchi0 ;
    DOUBLE etaKminus1 = eta*K - 1 ;
    DOUBLE Delta0=K*(eta*K-2) ;
    DOUBLE Delta1=-2*etaKminus1*(K+Delta0) ;
    DOUBLE Delta2=np.divide(1,2)*Delta1*(Delta1-4*etaKminus1)-a*a*etaKminus1*etaKminus1*Delta0 ;
    DOUBLE Delta3=-np.divide(1,3)*Delta1*Delta1*Delta1+etaKminus1*Delta1*Delta1+Delta2*Delta1-2*etaKminus1*(Delta2-etaKminus1)-a*a*etaKminus1*etaKminus1*Delta1 ;
    DOUBLE Delta4=np.divide(1,12)*(6*a*a*(Delta1*Delta1-2*Delta2)*etaKminus1*etaKminus1+3*Delta1*Delta1*Delta1*Delta1-8*etaKminus1*Delta1*Delta1*Delta1-12*Delta2*Delta1*Delta1+12*(2*etaKminus1*Delta2+Delta3)*Delta1+12*(np.divide(94,3)-np.divide(41,32)*M_PI*M_PI)*etaKminus1*etaKminus1+6*(Delta2*Delta2-4*Delta3*etaKminus1)) ;
    DOUBLE Delta5=etaKminus1*etaKminus1*(np.divide(-4237,60)+np.divide(128,5)*EMgamma+np.divide(2275,512)*M_PI*M_PI-np.divide(1,3)*a*a*(Delta1*Delta1*Delta1-3*Delta1*Delta2+3*Delta3)-(Delta1*Delta1*Delta1*Delta1*Delta1-5*Delta1*Delta1*Delta1*Delta2+5*Delta1*Delta2*Delta2+5*Delta1*Delta1*Delta3-5*Delta2*Delta3-5*Delta1*Delta4)/(5*etaKminus1*etaKminus1)+(Delta1*Delta1*Delta1*Delta1-4*Delta1*Delta1*Delta2+2*Delta2*Delta2+4*Delta1*Delta3-4*Delta4)/(2*etaKminus1)+np.divide(256,5)*log(2)+(np.divide(41,32)*M_PI*M_PI-np.divide(221,6))*eta) ;
    DOUBLE Delta5l=etaKminus1*etaKminus1*np.divide(64,5) ;
    DOUBLE logarg=u*(Delta1+u*(Delta2+u*(Delta3+u*(Delta4+u*(Delta5+Delta5l*log(u)))))) ;
    DOUBLE Deltaucalib = 1 + eta*(Delta0 + log1p( logarg)) ;
    DOUBLE Deltaucalibprime=-eta*u*u*(Delta1+u*(2*Delta2+u*(3*Delta3+u*(4*Delta4+u*(5*Delta5+Delta5l*(1+5*log(u)))))))/(1+logarg) ;
    DOUBLE Deltaubar=a*a*u*u+(2*u+1/etaKminus1)/etaKminus1 ;
    DOUBLE Deltaubarprime = -2*a*a*u*u*u - 2*u*u/(etaKminus1) ;
    DOUBLE Deltau=fabs(Deltaubar*Deltaucalib) ;
    DOUBLE Deltauprime = Deltaubarprime*Deltaucalib + Deltaubar*Deltaucalibprime ;
    DOUBLE Deltatprime=2*r*Deltau+r*r*Deltauprime ;
    DOUBLE Deltat=r*r*Deltau ;
    DOUBLE Deltar=Deltat*Dinv ;
    DOUBLE Lambdat=fabs(w2*w2-a*a*Deltat*sin2theta) ;
    DOUBLE csi=sqrt(Deltar*Deltat)/w2 ;
    DOUBLE csi1=1+(1-fabs(1-tortoise))*(csi-1) ;
    DOUBLE csi2=1+(np.divide(1,2)-np.divide(1,2)*copysign(np.divide(3,2)-tortoise))*(csi-1) ;
    DOUBLE prT=csi2*(p1*n1+p2*n2+p3*n3) ;
    DOUBLE phat3=p3-prT*(1-1/csi1)*n3 ;
    DOUBLE phat2 = p2 - prT*(1 - 1/csi1)*n2 ;
    DOUBLE phat1 = p1 - prT*(1 - 1/csi1)*n1 ;
    DOUBLE pdotxir=(phat1*xi1+phat2*xi2+phat3*xi3)*r ;
    DOUBLE pdotn=phat1*n1+phat2*n2+phat3*n3 ;
    DOUBLE pdotvr=(phat1*v1+phat2*v2+phat3*v3)*r ;
    DOUBLE pphi=pdotxir ;
    DOUBLE Qcoeff2=1/(Sigma*sin2theta) ;
    DOUBLE Qcoeff1=Sigma/(Lambdat*sin2theta) ;
    DOUBLE DrSipn2=Deltar*pdotn*pdotn/Sigma ;
    DOUBLE Q=1+DrSipn2+Qcoeff1*pdotxir*pdotxir+Qcoeff2*pdotvr*pdotvr ;
    DOUBLE Qminus1 = Q - 1 ;
    DOUBLE Jtilde=sqrt(Deltar) ;
    DOUBLE exp2mu=Sigma ;
    DOUBLE expmu = sqrt(exp2mu) ;
    DOUBLE Brtilde=(sqrt(Deltar)*Deltatprime-2*Deltat)/(2*sqrt(Deltar*Deltat)) ;
    DOUBLE Btilde=sqrt(Deltat) ;
    DOUBLE exp2nu=Deltat*Sigma/Lambdat ;
    DOUBLE expnu = sqrt(exp2nu) ;
    DOUBLE omega=omegatilde/Lambdat ;
    DOUBLE omegatildeprime=2*a ;
    DOUBLE Lambdatprime=4*(a*a+r*r)*r-a*a*Deltatprime*sin2theta ;
    DOUBLE mucostheta=a*a*costheta/Sigma ;
    DOUBLE nucostheta=a*a*w2*costheta*(w2-Deltat)/(Lambdat*Sigma) ;
    DOUBLE omegacostheta=-2*a*a*costheta*Deltat*omegatilde/(Lambdat*Lambdat) ;
    DOUBLE mur=r/Sigma-1/sqrt(Deltar) ;
    DOUBLE nur=r/Sigma+w2*(w2*Deltatprime-4*r*Deltat)/(2*Lambdat*Deltat) ;
    DOUBLE omegar=(Lambdat*omegatildeprime-Lambdatprime*omegatilde)/(Lambdat*Lambdat) ;
    DOUBLE dSO=147.481*chi*chi*chi*eta*eta-568.651*chi*chi*chi*eta+66.1987*chi*chi*chi-343.313*chi*chi*eta+2495.29*chi*eta*eta-44.5324 ;
    DOUBLE sigmacoeffTerm3 = eta*dSO*u*u*u ;
    DOUBLE sigmacoeffTerm2=eta/(144*r*r)*(-896+r*(-436*Qminus1-96*DrSipn2+r*(-45*Qminus1*Qminus1+36*Qminus1*DrSipn2))+eta*(-336+r*(204*Qminus1-882*DrSipn2+r*(810*DrSipn2*DrSipn2-234*Qminus1*DrSipn2)))) ;
    DOUBLE sigmacoeffTerm1=eta/12*(-8/r+3*Qminus1-36*DrSipn2) ;
    DOUBLE sigmacoeff=sigmacoeffTerm1+sigmacoeffTerm2+sigmacoeffTerm3 ;
    DOUBLE sigmastarcoeffTerm2=eta/(72*r*r)*(706+r*(-206*Qminus1+282*DrSipn2+r*Qminus1*(96*DrSipn2-23*Qminus1))+eta*(-54+r*(120*Qminus1-324*DrSipn2+r*(360*DrSipn2*DrSipn2+Qminus1*(-126*DrSipn2-3*Qminus1))))) ;
    DOUBLE sigmastarcoeffTerm1=eta/12*(14/r+4*Qminus1-30*DrSipn2) ;
    DOUBLE sigmastarcoeff=sigmastarcoeffTerm1+sigmastarcoeffTerm2 ;
    DOUBLE Deltasigmastar3=sigmastar3*sigmastarcoeff+sigma3*sigmacoeff ;
    DOUBLE Deltasigmastar2 = sigmastar2*sigmastarcoeff + sigma2*sigmacoeff ;
    DOUBLE Deltasigmastar1 = sigmastar1*sigmastarcoeff + sigma1*sigmacoeff ;
    DOUBLE Sstar3=sigmastar3+Deltasigmastar3 ;
    DOUBLE Sstar2 = sigmastar2 + Deltasigmastar2 ;
    DOUBLE Sstar1 = sigmastar1 + Deltasigmastar1 ;
    DOUBLE S3 = Sstar3 ;
    DOUBLE S2 = Sstar2 ;
    DOUBLE S1 = Sstar1 ;
    DOUBLE Sstardotn=Sstar1*n1+Sstar2*n2+Sstar3*n3 ;
    DOUBLE SdotSkerrhat=S1*Skerrhat1+S2*Skerrhat2+S3*Skerrhat3 ;
    DOUBLE Sdotn=S1*n1+S2*n2+S3*n3 ;
    DOUBLE Sdotv=S1*v1+S2*v2+S3*v3 ;
    DOUBLE Sdotxi=S1*xi1+S2*xi2+S3*xi3 ;
    DOUBLE HdsumTerm2=3*Sstardotn*Sstardotn ;
    DOUBLE HdsumTerm1=Sstar1*Sstar1+Sstar2*Sstar2+Sstar3*Sstar3 ;
    DOUBLE Hdsum=HdsumTerm1-HdsumTerm2 ;
    DOUBLE Hdcoeff=np.divide(1,2)/(r*r*r) ;
    DOUBLE Q4=2*prT*prT*prT*prT*u*u*(4-3*eta)*eta ;
    DOUBLE gammappsum=Deltar/Sigma*pdotn*pdotn+1/Sigma*pdotvr*pdotvr/sin2theta+Sigma/Lambdat/sin2theta*pdotxir*pdotxir ;
    DOUBLE Hnsradicand=1+gammappsum+Q4 ;
    DOUBLE alpha=sqrt(Deltat*Sigma/Lambdat) ;
    DOUBLE betapsum=omegatilde*pphi/Lambdat ;
    DOUBLE HssTerm3=expmu*expnu*pdotxir*(Jtilde*pdotn*Sdotxi*Btilde-expmu*expnu*pdotxir*Sdotn)+(pdotvr*(Sdotn*pdotvr-Jtilde*pdotn*Sdotv)-exp2mu*(sqrt(Q)+Q)*Sdotn*xisq)*Btilde*Btilde ;
    DOUBLE HssTerm3coeff=omegacostheta/(2*exp2mu*expmu*expnu*Btilde*(Q+sqrt(Q))) ;
    DOUBLE HssTerm2=expmu*pdotxir*(expmu*exp2nu*pdotxir*Sdotv-expnu*pdotvr*Sdotxi*Btilde)+xisq*Btilde*Btilde*(exp2mu*(sqrt(Q)+Q)*Sdotv+Jtilde*pdotn*(pdotvr*Sdotn-Jtilde*pdotn*Sdotv)) ;
    DOUBLE HssTerm2coeff=Jtilde*omegar/(2*exp2mu*expmu*expnu*Btilde*(Q+sqrt(Q))*xisq) ;
    DOUBLE HssTerm1=omega*SdotSkerrhat ;
    DOUBLE Hss=HssTerm1+HssTerm2coeff*HssTerm2+HssTerm3coeff*HssTerm3 ;
    DOUBLE HsoTerm2c=Jtilde*Brtilde*expmu*expnu*pdotxir*(sqrt(Q)+1)*Sdotv ;
    DOUBLE HsoTerm2b=expmu*expnu*pdotxir*(2*sqrt(Q)+1)*(Jtilde*nur*Sdotv-nucostheta*Sdotn*xisq)*Btilde ;
    DOUBLE HsoTerm2a=Sdotxi*Jtilde*(mur*pdotvr*(sqrt(Q)+1)-mucostheta*pdotn*xisq-sqrt(Q)*(nur*pdotvr+(mucostheta-nucostheta)*pdotn*xisq))*Btilde*Btilde ;
    DOUBLE HsoTerm2=HsoTerm2a+HsoTerm2b-HsoTerm2c ;
    DOUBLE HsoTerm2coeff=expnu/(exp2mu*Btilde*Btilde*(Q+sqrt(Q))*xisq) ;
    DOUBLE HsoTerm1=exp2nu*(expmu*expnu-Btilde)*pdotxir*SdotSkerrhat/(expmu*Btilde*Btilde*sqrt(Q)*xisq) ;
    DOUBLE Hso=HsoTerm1+HsoTerm2coeff*HsoTerm2 ;
    DOUBLE Hd=Hdcoeff*Hdsum ;
    DOUBLE Hns=betapsum+alpha*sqrt(Hnsradicand) ;
    DOUBLE Hs=Hso+Hss ;
    DOUBLE dSS=528.511*chi*chi*chi*eta*eta-41.0003*chi*chi*chi*eta+1161.78*chi*chi*eta*eta*eta-326.325*chi*chi*eta*eta+37.1964*chi*eta+706.958*eta*eta*eta-36.0272*eta+6.06807 ;
    DOUBLE Heff = Hs + Hns - Hd + dSS*eta*u*u*u*u*(S1x*S1x + S1y*S1y + S1z*S1z + S2x*S2x + S2y*S2y + S2z*S2z) ;
    DOUBLE Hreal=sqrt(1+2*eta*(Heff-1)) ;
    return Hreal;
}
#include math.h 
double compute_v4P_Hreal(double m1,double m2,INT tortoise,double x,double y,double z,double p1,double p2,double p3,double S1x,double S1y,double S1z,double S2x,double S2y,double S2z){
    double EMgamma = 0.577215664901532860606512090082402431 ;
    double M=m1+m2 ;
    double mu=m1*m2/M ;
    double eta=mu/M ;
    double r=sqrt(x*x+y*y+z*z) ;
    double u=1/r ;
    double sigmastar3=m2/m1*S1z+m1/m2*S2z ;
    double sigmastar2 = m2/m1*S1y + m1/m2*S2y ;
    double sigmastar1 = m2/m1*S1x + m1/m2*S2x ;
    double sigma3=S1z+S2z ;
    double sigma2 = S1y + S2y ;
    double sigma1 = S1x + S2x ;
    double Skerr3=sigma3 ;
    double Skerr2 = sigma2 ;
    double Skerr1 = sigma1 ;
    double Skerrmag=sqrt(Skerr1*Skerr1+Skerr2*Skerr2+Skerr3*Skerr3) ;
    double Skerrhat3=Skerr3/Skerrmag ;
    double Skerrhat2 = Skerr2/Skerrmag ;
    double Skerrhat1 = Skerr1/Skerrmag ;
    double a=Skerrmag ;
    double L3=x*p2-y*p1 ;
    double L2 = z*p1 - x*p3 ;
    double L1 = y*p3 - z*p2 ;
    double Lnorm=sqrt(L1*L1+L2*L2+L3*L3) ;
    double Lhat3=L3/Lnorm ;
    double Lhat2 = L2/Lnorm ;
    double Lhat1 = L1/Lnorm ;
    double S2dotLhat=S2x*Lhat1+S2y*Lhat2+S2z*Lhat3 ;
    double S1dotLhat = S1x*Lhat1 + S1y*Lhat2 + S1z*Lhat3 ;
    double S1perp3=S1z-S1dotLhat*Lhat3 ;
    double S1perp2 = S1y - S1dotLhat*Lhat2 ;
    double S1perp1 = S1x - S1dotLhat*Lhat1 ;
    double S2perp3=S2z-S2dotLhat*Lhat3 ;
    double S2perp2 = S2y - S2dotLhat*Lhat2 ;
    double S2perp1 = S2x - S2dotLhat*Lhat1 ;
    double Sperp3=S1perp3+S2perp3 ;
    double Sperp2 = S1perp2 + S2perp2 ;
    double Sperp1 = S1perp1 + S2perp1 ;
    double n3=z/r ;
    double n2 = y/r ;
    double n1 = x/r ;
    double e33=Skerrhat3 ;
    double e32 = Skerrhat2 ;
    double e31 = Skerrhat1 ;
    double xi3=e31*n2-e32*n1 ;
    double xi2 = -e31*n3 + e33*n1 ;
    double xi1 = e32*n3 - e33*n2 ;
    double v3=n1*xi2-n2*xi1 ;
    double v2 = n3*xi1 - n1*xi3 ;
    double v1 = n2*xi3 - n3*xi2 ;
    double costheta=e31*n1+e32*n2+e33*n3 ;
    double sin2theta=1-costheta*costheta ;
    double xisq = sin2theta ;
    double w2=a*a+r*r ;
    double Sigma=r*r+a*a*costheta*costheta ;
    double Dinv=1+log(1+6*eta*u*u+2*(26-3*eta)*eta*u*u*u) ;
    double omegatilde=2*a*r ;
    double chi=(Skerr1*Lhat1+Skerr2*Lhat2+Skerr3*Lhat3)/(1-2*eta)+np.divide(1,2)*(Sperp1*Skerr1+Sperp2*Skerr2+Sperp3*Skerr3)/(Skerrmag*(1.-2.*eta)) ;
    double Kchi0=267.788*eta*eta*eta-126.687*eta*eta+10.2573*eta+1.7336 ;
    double K=-59.1658*chi*chi*chi*eta*eta*eta-0.426958*chi*chi*chi*eta+1.43659*chi*chi*chi+31.1746*chi*chi*eta*eta*eta+6.16466*chi*chi*eta*eta-1.38086*chi*chi-27.5201*chi*eta*eta*eta+17.3736*chi*eta*eta+2.26831*chi*eta-1.62045*chi+Kchi0 ;
    double etaKminus1 = eta*K - 1 ;
    double Delta0=K*(eta*K-2) ;
    double Delta1=-2*etaKminus1*(K+Delta0) ;
    double Delta2=np.divide(1,2)*Delta1*(Delta1-4*etaKminus1)-a*a*etaKminus1*etaKminus1*Delta0 ;
    double Delta3=-np.divide(1,3)*Delta1*Delta1*Delta1+etaKminus1*Delta1*Delta1+Delta2*Delta1-2*etaKminus1*(Delta2-etaKminus1)-a*a*etaKminus1*etaKminus1*Delta1 ;
    double Delta4=np.divide(1,12)*(6*a*a*(Delta1*Delta1-2*Delta2)*etaKminus1*etaKminus1+3*Delta1*Delta1*Delta1*Delta1-8*etaKminus1*Delta1*Delta1*Delta1-12*Delta2*Delta1*Delta1+12*(2*etaKminus1*Delta2+Delta3)*Delta1+12*(np.divide(94,3)-np.divide(41,32)*M_PI*M_PI)*etaKminus1*etaKminus1+6*(Delta2*Delta2-4*Delta3*etaKminus1)) ;
    double Delta5=etaKminus1*etaKminus1*(np.divide(-4237,60)+np.divide(128,5)*EMgamma+np.divide(2275,512)*M_PI*M_PI-np.divide(1,3)*a*a*(Delta1*Delta1*Delta1-3*Delta1*Delta2+3*Delta3)-(Delta1*Delta1*Delta1*Delta1*Delta1-5*Delta1*Delta1*Delta1*Delta2+5*Delta1*Delta2*Delta2+5*Delta1*Delta1*Delta3-5*Delta2*Delta3-5*Delta1*Delta4)/(5*etaKminus1*etaKminus1)+(Delta1*Delta1*Delta1*Delta1-4*Delta1*Delta1*Delta2+2*Delta2*Delta2+4*Delta1*Delta3-4*Delta4)/(2*etaKminus1)+np.divide(256,5)*log(2)+(np.divide(41,32)*M_PI*M_PI-np.divide(221,6))*eta) ;
    double Delta5l=etaKminus1*etaKminus1*np.divide(64,5) ;
    double logarg=u*(Delta1+u*(Delta2+u*(Delta3+u*(Delta4+u*(Delta5+Delta5l*log(u)))))) ;
    double Deltaucalib = 1 + eta*(Delta0 + log1p( logarg)) ;
    double Deltaucalibprime=-eta*u*u*(Delta1+u*(2*Delta2+u*(3*Delta3+u*(4*Delta4+u*(5*Delta5+Delta5l*(1+5*log(u)))))))/(1+logarg) ;
    double Deltaubar=a*a*u*u+(2*u+1/etaKminus1)/etaKminus1 ;
    double Deltaubarprime = -2*a*a*u*u*u - 2*u*u/(etaKminus1) ;
    double Deltau=fabs(Deltaubar*Deltaucalib) ;
    double Deltauprime = Deltaubarprime*Deltaucalib + Deltaubar*Deltaucalibprime ;
    double Deltatprime=2*r*Deltau+r*r*Deltauprime ;
    double Deltat=r*r*Deltau ;
    double Deltar=Deltat*Dinv ;
    double Lambdat=fabs(w2*w2-a*a*Deltat*sin2theta) ;
    double csi=sqrt(Deltar*Deltat)/w2 ;
    double csi1=1+(1-fabs(1-(double)tortoise))*(csi-1) ;
    double csi2=1+(np.divide(1,2)-np.divide(1,2)*copysign(np.divide(3,2)-(double)tortoise))*(csi-1) ;
    double prT=csi2*(p1*n1+p2*n2+p3*n3) ;
    double phat3=p3-prT*(1-1/csi1)*n3 ;
    double phat2 = p2 - prT*(1 - 1/csi1)*n2 ;
    double phat1 = p1 - prT*(1 - 1/csi1)*n1 ;
    double pdotxir=(phat1*xi1+phat2*xi2+phat3*xi3)*r ;
    double pdotn=phat1*n1+phat2*n2+phat3*n3 ;
    double pdotvr=(phat1*v1+phat2*v2+phat3*v3)*r ;
    double pphi=pdotxir ;
    double Qcoeff2=1/(Sigma*sin2theta) ;
    double Qcoeff1=Sigma/(Lambdat*sin2theta) ;
    double DrSipn2=Deltar*pdotn*pdotn/Sigma ;
    double Q=1+DrSipn2+Qcoeff1*pdotxir*pdotxir+Qcoeff2*pdotvr*pdotvr ;
    double Qminus1 = Q - 1 ;
    double Jtilde=sqrt(Deltar) ;
    double exp2mu=Sigma ;
    double expmu = sqrt(exp2mu) ;
    double Brtilde=(sqrt(Deltar)*Deltatprime-2*Deltat)/(2*sqrt(Deltar*Deltat)) ;
    double Btilde=sqrt(Deltat) ;
    double exp2nu=Deltat*Sigma/Lambdat ;
    double expnu = sqrt(exp2nu) ;
    double omega=omegatilde/Lambdat ;
    double omegatildeprime=2*a ;
    double Lambdatprime=4*(a*a+r*r)*r-a*a*Deltatprime*sin2theta ;
    double mucostheta=a*a*costheta/Sigma ;
    double nucostheta=a*a*w2*costheta*(w2-Deltat)/(Lambdat*Sigma) ;
    double omegacostheta=-2*a*a*costheta*Deltat*omegatilde/(Lambdat*Lambdat) ;
    double mur=r/Sigma-1/sqrt(Deltar) ;
    double nur=r/Sigma+w2*(w2*Deltatprime-4*r*Deltat)/(2*Lambdat*Deltat) ;
    double omegar=(Lambdat*omegatildeprime-Lambdatprime*omegatilde)/(Lambdat*Lambdat) ;
    double dSO=147.481*chi*chi*chi*eta*eta-568.651*chi*chi*chi*eta+66.1987*chi*chi*chi-343.313*chi*chi*eta+2495.29*chi*eta*eta-44.5324 ;
    double sigmacoeffTerm3 = eta*dSO*u*u*u ;
    double sigmacoeffTerm2=eta/(144*r*r)*(-896+r*(-436*Qminus1-96*DrSipn2+r*(-45*Qminus1*Qminus1+36*Qminus1*DrSipn2))+eta*(-336+r*(204*Qminus1-882*DrSipn2+r*(810*DrSipn2*DrSipn2-234*Qminus1*DrSipn2)))) ;
    double sigmacoeffTerm1=eta/12*(-8/r+3*Qminus1-36*DrSipn2) ;
    double sigmacoeff=sigmacoeffTerm1+sigmacoeffTerm2+sigmacoeffTerm3 ;
    double sigmastarcoeffTerm2=eta/(72*r*r)*(706+r*(-206*Qminus1+282*DrSipn2+r*Qminus1*(96*DrSipn2-23*Qminus1))+eta*(-54+r*(120*Qminus1-324*DrSipn2+r*(360*DrSipn2*DrSipn2+Qminus1*(-126*DrSipn2-3*Qminus1))))) ;
    double sigmastarcoeffTerm1=eta/12*(14/r+4*Qminus1-30*DrSipn2) ;
    double sigmastarcoeff=sigmastarcoeffTerm1+sigmastarcoeffTerm2 ;
    double Deltasigmastar3=sigmastar3*sigmastarcoeff+sigma3*sigmacoeff ;
    double Deltasigmastar2 = sigmastar2*sigmastarcoeff + sigma2*sigmacoeff ;
    double Deltasigmastar1 = sigmastar1*sigmastarcoeff + sigma1*sigmacoeff ;
    double Sstar3=sigmastar3+Deltasigmastar3 ;
    double Sstar2 = sigmastar2 + Deltasigmastar2 ;
    double Sstar1 = sigmastar1 + Deltasigmastar1 ;
    double S3 = Sstar3 ;
    double S2 = Sstar2 ;
    double S1 = Sstar1 ;
    double Sstardotn=Sstar1*n1+Sstar2*n2+Sstar3*n3 ;
    double SdotSkerrhat=S1*Skerrhat1+S2*Skerrhat2+S3*Skerrhat3 ;
    double Sdotn=S1*n1+S2*n2+S3*n3 ;
    double Sdotv=S1*v1+S2*v2+S3*v3 ;
    double Sdotxi=S1*xi1+S2*xi2+S3*xi3 ;
    double HdsumTerm2=3*Sstardotn*Sstardotn ;
    double HdsumTerm1=Sstar1*Sstar1+Sstar2*Sstar2+Sstar3*Sstar3 ;
    double Hdsum=HdsumTerm1-HdsumTerm2 ;
    double Hdcoeff=np.divide(1,2)/(r*r*r) ;
    double Q4=2*prT*prT*prT*prT*u*u*(4-3*eta)*eta ;
    double gammappsum=Deltar/Sigma*pdotn*pdotn+1/Sigma*pdotvr*pdotvr/sin2theta+Sigma/Lambdat/sin2theta*pdotxir*pdotxir ;
    double Hnsradicand=1+gammappsum+Q4 ;
    double alpha=sqrt(Deltat*Sigma/Lambdat) ;
    double betapsum=omegatilde*pphi/Lambdat ;
    double HssTerm3=expmu*expnu*pdotxir*(Jtilde*pdotn*Sdotxi*Btilde-expmu*expnu*pdotxir*Sdotn)+(pdotvr*(Sdotn*pdotvr-Jtilde*pdotn*Sdotv)-exp2mu*(sqrt(Q)+Q)*Sdotn*xisq)*Btilde*Btilde ;
    double HssTerm3coeff=omegacostheta/(2*exp2mu*expmu*expnu*Btilde*(Q+sqrt(Q))) ;
    double HssTerm2=expmu*pdotxir*(expmu*exp2nu*pdotxir*Sdotv-expnu*pdotvr*Sdotxi*Btilde)+xisq*Btilde*Btilde*(exp2mu*(sqrt(Q)+Q)*Sdotv+Jtilde*pdotn*(pdotvr*Sdotn-Jtilde*pdotn*Sdotv)) ;
    double HssTerm2coeff=Jtilde*omegar/(2*exp2mu*expmu*expnu*Btilde*(Q+sqrt(Q))*xisq) ;
    double HssTerm1=omega*SdotSkerrhat ;
    double Hss=HssTerm1+HssTerm2coeff*HssTerm2+HssTerm3coeff*HssTerm3 ;
    double HsoTerm2c=Jtilde*Brtilde*expmu*expnu*pdotxir*(sqrt(Q)+1)*Sdotv ;
    double HsoTerm2b=expmu*expnu*pdotxir*(2*sqrt(Q)+1)*(Jtilde*nur*Sdotv-nucostheta*Sdotn*xisq)*Btilde ;
    double HsoTerm2a=Sdotxi*Jtilde*(mur*pdotvr*(sqrt(Q)+1)-mucostheta*pdotn*xisq-sqrt(Q)*(nur*pdotvr+(mucostheta-nucostheta)*pdotn*xisq))*Btilde*Btilde ;
    double HsoTerm2=HsoTerm2a+HsoTerm2b-HsoTerm2c ;
    double HsoTerm2coeff=expnu/(exp2mu*Btilde*Btilde*(Q+sqrt(Q))*xisq) ;
    double HsoTerm1=exp2nu*(expmu*expnu-Btilde)*pdotxir*SdotSkerrhat/(expmu*Btilde*Btilde*sqrt(Q)*xisq) ;
    double Hso=HsoTerm1+HsoTerm2coeff*HsoTerm2 ;
    double Hd=Hdcoeff*Hdsum ;
    double Hns=betapsum+alpha*sqrt(Hnsradicand) ;
    double Hs=Hso+Hss ;
    double dSS=528.511*chi*chi*chi*eta*eta-41.0003*chi*chi*chi*eta+1161.78*chi*chi*eta*eta*eta-326.325*chi*chi*eta*eta+37.1964*chi*eta+706.958*eta*eta*eta-36.0272*eta+6.06807 ;
    double Heff = Hs + Hns - Hd + dSS*eta*u*u*u*u*(S1x*S1x + S1y*S1y + S1z*S1z + S2x*S2x + S2y*S2y + S2z*S2z) ;
    double Hreal=sqrt(1+2*eta*(Heff-1)) ;
    return Hreal;
}
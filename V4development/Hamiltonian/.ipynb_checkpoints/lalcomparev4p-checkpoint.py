import numpy as np
def compute_Hreal_and_csi(m1=23., m2=10., EMgamma=0.577215664901532860606512090082402431, tortoise=1, x=2.129681018601393e+01, y=0.000000000000000e+00, z=0.000000000000000e+00, p1=0.000000000000000e+00, p2=2.335391115580442e-01, p3=-4.235164736271502e-22, S1x=4.857667584940312e-03, S1y=9.715161660389764e-03, S1z=-1.457311842632286e-02, S2x=3.673094582185491e-03, S2y=-4.591302628615413e-03, S2z=5.509696538546906e-03):
    EMgamma = 0.577215664901532860606512090082402431
    M=m1+m2
    mu=m1*m2/M
    eta=mu/M
    
    r=np.sqrt(x*x+y*y+z*z)
    # lal version:246-247
    #r2 = x->data[0]*x->data[0] + x->data[1]*x->data[1] + x->data[2]*x->data[2];
    #r  = sqrt(r2);
    
    u=1/r
    #lal version: 248
    #u  = 1./r;
    
    sigmastar3=m2/m1*S1z+m1/m2*S2z
    sigmastar2 = m2/m1*S1y + m1/m2*S2y
    sigmastar1 = m2/m1*S1x + m1/m2*S2x

    sigma2 = S1y + S2y
    sigma1 = S1x + S2x
    Skerr3=sigma3
    Skerr2 = sigma2
    Skerr1 = sigma1
    
    Skerrmag=np.sqrt(Skerr1*Skerr1+Skerr2*Skerr2+Skerr3*Skerr3)
    Skerrhat3=Skerr3/Skerrmag
    Skerrhat2 = Skerr2/Skerrmag
    Skerrhat1 = Skerr1/Skerrmag
    a=Skerrmag
    #lal version:266-267
    #a2 = sKerr_x*sKerr_x + sKerr_y*sKerr_y + sKerr_z*sKerr_z;
    #a  = sqrt( a2 );
    
    
    L3=x*p2-y*p1
    L2 = z*p1 - x*p3
    L1 = y*p3 - z*p2
    #lal version: 160
    #cross_product(x->data, p->data, L); // Note that L = r x p is invariant under tortoise transform
    
    Lnorm=np.sqrt(L1*L1+L2*L2+L3*L3)
    #lal version: 161
    #REAL8 L_mag = sqrt(inner_product(L, L));
    
    Lhat3=L3/Lnorm
    Lhat2 = L2/Lnorm
    Lhat1 = L1/Lnorm
    #lal version: 162-165
    #for (UINT4 jj = 0; jj < 3; jj++)
    #{
    #  Lhat[jj] = L[jj] / L_mag;
    #}
    
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
    #lal version: 167-177
    #REAL8 tempS1_p = inner_product(s1Vec->data, Lhat);
    #REAL8 tempS2_p = inner_product(s2Vec->data, Lhat);
    #REAL8 S1_perp[3] = {0, 0, 0};
    #REAL8 S2_perp[3] = {0, 0, 0};
    #REAL8 S_perp[3] = {0,0,0};
    #for (UINT4 jj = 0; jj < 3; jj++)
    #{
    # S1_perp[jj] = s1Vec->data[jj] - tempS1_p * Lhat[jj];
    # S2_perp[jj] = s2Vec->data[jj] - tempS2_p * Lhat[jj];
    # S_perp[jj] = S1_perp[jj]+S2_perp[jj];
    #}
    
    
    n3=z/r
    n2 = y/r
    n1 = x/r
    #lal version: 254-256
    #nx = x->data[0] *u;
    #ny = x->data[1] *u;
    #nz = x->data[2] *u;
    
    e33=Skerrhat3
    e32 = Skerrhat2
    e31 = Skerrhat1
    #lal version: 269-275
    #if(a !=0.)
    #{
    # const REAL8 inva = 1./a;
    # e3_x = sKerr_x * inva;
    # e3_y = sKerr_y * inva;
    # e3_z = sKerr_z * inva;
    #}
    
    xi3=e31*n2-e32*n1
    xi2 = -e31*n3 + e33*n1
    xi1 = e32*n3 - e33*n2
    #lal version: 315-317
    #xi_x = -e3_z*ny + e3_y*nz;
    #xi_y =  e3_z*nx - e3_x*nz;
    #xi_z = -e3_y*nx + e3_x*ny;
    
    v3=n1*xi2-n2*xi1
    v2 = n3*xi1 - n1*xi3
    v1 = n2*xi3 - n3*xi2
    #lal version: 319-321
    #vx = -nz*xi_y + ny*xi_z;
    #vy =  nz*xi_x - nx*xi_z;
    #vz = -ny*xi_x + nx*xi_y;
    
    costheta=e31*n1+e32*n2+e33*n3
    sin2theta=1-costheta*costheta
    xisq = sin2theta
    #lal version: 311-313
    #costheta = e3_x*nx + e3_y*ny + e3_z*nz;
    #
    #xi2=1. - costheta*costheta;
    
    w2=a*a+r*r
    Sigma=r*r+a*a*costheta*costheta
    #lal version: 323-324
    #w2 = r2 + a2;
    #rho2 = r2 + a2*costheta*costheta;
    # note: rho2 (lal) == Sigma (nrpy)
    
    Dinv=1+np.log(1+6*eta*u*u+2*(26-3*eta)*eta*u*u*u)
    #lal version: 359
    #D = 1. + log1p(6.*eta*u2 + 2.*(26. - 3.*eta)*eta*u3);
    
    omegatilde=2*a*r
    #lal version: 365 (coeffs map to zero)
    #ww=2.*a*r + coeffs->b3*eta*a2*a*u + coeffs->bb3*eta*a*u;
    
    chi=(Skerr1*Lhat1+Skerr2*Lhat2+Skerr3*Lhat3)/(1-2*eta)+np.divide(1,2)*(Sperp1*Skerr1+Sperp2*Skerr2+Sperp3*Skerr3)/(Skerrmag*(1.-2.*eta))
    #lal version: 179-187
    #REAL8 S_con = 0.0;
    #if (sKerr_norm>1e-6){
    # S_con = sigmaKerr->data[0] * Lhat[0] + sigmaKerr->data[1] * Lhat[1] + sigmaKerr->data[2] * Lhat[2];
    # S_con /= (1 - 2 * eta);
    # // Last division by 2 is to ensure the spin oebys the Kerr bound.
    # S_con += inner_product(S_perp, sigmaKerr->data) / sKerr_norm / (1 - 2 * eta) / 2.;
    #}
    #
    #REAL8 chi = S_con;
    
    Kchi0=267.788*eta*eta*eta-126.687*eta*eta+10.2573*eta+1.7336
    K=-59.1658*chi*chi*chi*eta*eta*eta-0.426958*chi*chi*chi*eta+1.43659*chi*chi*chi+31.1746*chi*chi*eta*eta*eta+6.16466*chi*chi*eta*eta-1.38086*chi*chi-27.5201*chi*eta*eta*eta+17.3736*chi*eta*eta+2.26831*chi*eta-1.62045*chi+Kchi0
    etaKminus1 = eta*K - 1
    #lal version: [SpinPrecEOBHCoeffs.c] 226-233, [SpinEOBHamiltonian.h]14-29. coeffnmK -> coeff of eta^nchi^m
    #static const REAL8 coeff00K = 1.7336;
    #static const REAL8 coeff01K = -1.62045;
    #static const REAL8 coeff02K = -1.38086;
    #static const REAL8 coeff03K = 1.43659;
    #static const REAL8 coeff10K = 10.2573;
    #static const REAL8 coeff11K = 2.26831;
    #static const REAL8 coeff12K = 0;
    #static const REAL8 coeff13K = -0.426958;
    #static const REAL8 coeff20K = -126.687;
    #static const REAL8 coeff21K = 17.3736;
    #static const REAL8 coeff22K = 6.16466;
    #static const REAL8 coeff23K = 0;
    #static const REAL8 coeff30K = 267.788;
    #static const REAL8 coeff31K = -27.5201;
    #static const REAL8 coeff32K = 31.1746;
    #static const REAL8 coeff33K = -59.1658;
    
    #lal version:[SpinPrecEOBHCoeffs.c] 236-250, kn in lal == Deltan in nrpy, KK in lal == K in nrpy, m1PlusEtaKK in lal 
    Delta0=K*(eta*K-2)
    #coeffs->k0 = k0 = KK * (m1PlusEtaKK - 1.);
    
    Delta1=-2*etaKminus1*(K+Delta0)
    #coeffs->k1 = k1 = -2. * (k0 + KK) * m1PlusEtaKK;
    #k1p2 = k1 * k1;
    #k1p3 = k1 * k1p2;
    
    Delta2=np.divide(1,2)*Delta1*(Delta1-4*etaKminus1)-a*a*etaKminus1*etaKminus1*Delta0
    #coeffs->k2 = k2 = (k1 * (k1 - 4. * m1PlusEtaKK)) * 0.5 - a * a * k0 * m1PlusEtaKK * m1PlusEtaKK;
    
    Delta3=-np.divide(1,3)*Delta1*Delta1*Delta1+etaKminus1*Delta1*Delta1+Delta2*Delta1-2*etaKminus1*(Delta2-etaKminus1)-a*a*etaKminus1*etaKminus1*Delta1
    #coeffs->k3 = k3 = -(k1 * k1) * k1 * third + k1 * k2 + (k1 * k1) * m1PlusEtaKK - 2. * (k2 - m1PlusEtaKK) * m1PlusEtaKK - a * a * k1 * (m1PlusEtaKK * m1PlusEtaKK);
    
    Delta4=np.divide(1,12)*(6*a*a*(Delta1*Delta1-2*Delta2)*etaKminus1*etaKminus1+3*Delta1*Delta1*Delta1*Delta1-8*etaKminus1*Delta1*Delta1*Delta1-12*Delta2*Delta1*Delta1+12*(2*etaKminus1*Delta2+Delta3)*Delta1+12*(np.divide(94,3)-np.divide(41,32)*np.pi*np.pi)*etaKminus1*etaKminus1+6*(Delta2*Delta2-4*Delta3*etaKminus1))
    #coeffs->k4 = k4 = ((24. / 96.) * (k1 * k1) * (k1 * k1) - (96. / 96.) * (k1 * k1) * k2 + (48. / 96.) * k2 * k2 - (64. / 96.) * (k1 * k1) * k1 * m1PlusEtaKK + (48. / 96.) * (a * a) * (k1 * k1 - 2. * k2) * (m1PlusEtaKK * m1PlusEtaKK) +
    #                  (96. / 96.) * k1 * (k3 + 2. * k2 * m1PlusEtaKK) - m1PlusEtaKK * ((192. / 96.) * k3 + m1PlusEtaKK * (-(3008. / 96.) + (123. / 96.) * LAL_PI * LAL_PI)));
    
    Delta5=etaKminus1*etaKminus1*(np.divide(-4237,60)+np.divide(128,5)*EMgamma+np.divide(2275,512)*np.pi*np.pi-np.divide(1,3)*a*a*(Delta1*Delta1*Delta1-3*Delta1*Delta2+3*Delta3)-(Delta1*Delta1*Delta1*Delta1*Delta1-5*Delta1*Delta1*Delta1*Delta2+5*Delta1*Delta2*Delta2+5*Delta1*Delta1*Delta3-5*Delta2*Delta3-5*Delta1*Delta4)/(5*etaKminus1*etaKminus1)+(Delta1*Delta1*Delta1*Delta1-4*Delta1*Delta1*Delta2+2*Delta2*Delta2+4*Delta1*Delta3-4*Delta4)/(2*etaKminus1)+np.divide(256,5)*np.log(2)+(np.divide(41,32)*np.pi*np.pi-np.divide(221,6))*eta)
    #coeffs->k5 = k5 = m1PlusEtaKK * m1PlusEtaKK * (-4237. / 60. + 128. / 5. * LAL_GAMMA + 2275. * LAL_PI * LAL_PI / 512. - third * (a * a) * (k1p3 - 3. * (k1 * k2) + 3. * k3) - ((k1p3 * k1p2) - 5. * (k1p3 * k2) + 5. * k1 * k2 * k2 + 5. * k1p2 * k3 - 5. * k2 * k3 - 5. * k1 * k4) * fifth * invm1PlusEtaKK * invm1PlusEtaKK + ((k1p2 * k1p2) - 4. * (k1p2 * k2) + 2. * k2 * k2 + 4. * k1 * k3 - 4. * k4) * 0.5 * invm1PlusEtaKK + (256. / 5.) * ln2 + (41. * LAL_PI * LAL_PI / 32. - 221. / 6.) * eta);
    
    Delta5l=etaKminus1*etaKminus1*np.divide(64,5)
    #coeffs->k5l = k5l = (m1PlusEtaKK * m1PlusEtaKK) * (64. / 5.);        
    
    logarg=u*(Delta1+u*(Delta2+u*(Delta3+u*(Delta4+u*(Delta5+Delta5l*np.log(u))))))
    Deltaucalib = 1 + eta*(Delta0 + np.log(1 + logarg))
    Deltaucalibprime=-eta*u*u*(Delta1+u*(2*Delta2+u*(3*Delta3+u*(4*Delta4+u*(5*(Delta5+Delta5l*np.log(u)))))))/(1+logarg)        
    Deltaubar=a*a*u*u+(2*u+1/etaKminus1)/etaKminus1    
    Deltaubarprime = -2*a*a*u*u*u - 2*u*u/(etaKminus1)
    Deltau=np.abs(Deltaubar*Deltaucalib)
    Deltauprime = Deltaubarprime*Deltaucalib + Deltaubar*Deltaucalibprime
    #lal version: 329-343,346-348
    #bulk = invm1PlusetaKK*(invm1PlusetaKK + (2.*u)) + a2*u2;
    #/* Eq. 5.73 of BB1 */
    #// use ln(u) = log_2(u)/log_2(e) and the fact that log2 is faster than ln
    #// this relies on the compiler evaluating the expression at compile time.
    #// which apparently not all do so in stead of 1./log2(exp(1.)) I use the
    #// result returned by Maple.
    #const REAL8 invlog_2e = 0.69314718055994530941723212145817656807550013436026;
    #logu = log2(u)*invlog_2e;
    #const REAL8 logarg = coeffs->k1*u + coeffs->k2*u2 + coeffs->k3*u3 + coeffs->k4*u4
    #                                          + coeffs->k5*u5 + coeffs->k5l*u5*logu;
    #logTerms = 1. + eta*coeffs->k0 + eta*log1p(fabs(1. + logarg) - 1.);
    #if(debugPK)XLAL_PRINT_INFO( "bulk = %.16e, logTerms = %.16e\n", bulk, logTerms );
    #/* Eq. 5.73 of BB1 */
    #deltaU = fabs(bulk*logTerms);
    #deltaU_u = 2.*(invm1PlusetaKK + a2*u)*logTerms +
    #       bulk * (eta*(coeffs->k1 + u*(2.*coeffs->k2 + u*(3.*coeffs->k3 + u*(4.*coeffs->k4 + 5.*(coeffs->k5+coeffs->k5l*logu)*u)))))
    #       / (1. + logarg);
    #BUG!!!!! logu not differentiated in k5l*u5*logu
    
    
    Deltatprime=2*r*Deltau+r*r*Deltauprime
    Deltat=r*r*Deltau
    #lalversion: 344,350
    #deltaT = r2*deltaU;
    #deltaT_r = 2.*r*deltaU - deltaU_u;
    
    Deltar=Deltat*Dinv
    #lalversion: 361
    #deltaR = deltaT*D;
    
    Lambdat=np.abs(w2*w2-a*a*Deltat*sin2theta)
    #lalversion: 352
    #Lambda = fabs(w2*w2 - a2*deltaT*xi2);
    
    csi=np.sqrt(Deltar*Deltat)/w2
    csi1=1+(1-np.abs(1-tortoise))*(csi-1)
    csi2=1+(np.divide(1,2)-np.divide(1,2)*np.sign(np.divide(3,2)-tortoise))*(csi-1)
    #lalversion: 370-374
    #csi = sqrt( fabs(deltaT * deltaR) )/ w2;
    #// non-unity only for tortoise==1
    #const REAL8 csi1 = 1.0 + (1.-fabs(1.-tortoise)) * (csi - 1.0);
    #// non-unity only for tortoise==2
    #const REAL8 csi2 = 1.0 + (0.5-copysign(0.5, 1.5-tortoise)) * (csi - 1.0);
    
    
    prT=csi2*(p1*n1+p2*n2+p3*n3)
    phat3=p3-prT*(1-1/csi1)*n3
    phat2 = p2 - prT*(1 - 1/csi1)*n2
    phat1 = p1 - prT*(1 - 1/csi1)*n1
    #lalversion: 380-384
    #prT = (p->data[0]*nx + p->data[1]*ny + p->data[2]*nz)*csi2;
    #/* p->data is BL momentum vector; tmpP is tortoise momentum vector */
    #tmpP[0] = p->data[0] - nx * prT * (1. - 1./csi1);
    #tmpP[1] = p->data[1] - ny * prT * (1. - 1./csi1);
    #tmpP[2] = p->data[2] - nz * prT * (1. - 1./csi1);
    
    pdotxir=(phat1*xi1+phat2*xi2+phat3*xi3)*r
    pdotn=phat1*n1+phat2*n2+phat3*n3
    pdotvr=(phat1*v1+phat2*v2+phat3*v3)*r
    pphi=pdotxir
    #lal version: 386-392
    #pxir = (tmpP[0]*xi_x + tmpP[1]*xi_y + tmpP[2]*xi_z) * r;
    #pvr  = (tmpP[0]*vx + tmpP[1]*vy + tmpP[2]*vz) * r;
    #pn   = tmpP[0]*nx + tmpP[1]*ny + tmpP[2]*nz;
    #
    #pr = pn;
    #pf = pxir;
    #ptheta2 = pvr * pvr *invxi2;
    
    Qcoeff2=1/(Sigma*sin2theta)
    Qcoeff1=Sigma/(Lambdat*sin2theta)
    DrSipn2=Deltar*pdotn*pdotn/Sigma
    Q=1+DrSipn2+Qcoeff1*pdotxir*pdotxir+Qcoeff2*pdotvr*pdotvr
    Qminus1 = Q - 1
    #lal version:446-451 
    #Q = 1. + pvr*pvr*invrho2*invxi2 + pxir*pxir*rho2*invLambda*invxi2 + pn*pn*deltaR*invrho2;
    #pn2 = pr * pr * deltaR * invrho2;
    #pp  = Q - 1.;
    
    Jtilde=np.sqrt(Deltar)
    #lal version: 418
    #const REAL8 sqrtdeltaR = sqrt(deltaR);

    exp2mu=Sigma
    expmu = np.sqrt(exp2mu)
    #lal version: 427
    #const REAL8 expMU = sqrt(rho2);
    
    Brtilde=(np.sqrt(Deltar)*Deltatprime-2*Deltat)/(2*np.sqrt(Deltar*Deltat))
    #lal version: 437
    #BR = (-deltaT*invsqrtdeltaR + deltaT_r*0.5)*invsqrtdeltaT;
    
    Btilde=np.sqrt(Deltat)
    #lal version: 415
    #B = sqrt(deltaT);
    
    exp2nu=Deltat*Sigma/Lambdat
    expnu = np.sqrt(exp2nu)
    #lal version: 426
    #const REAL8 expnu = sqrt(deltaT*rho2*invLambda);
    
    omega=omegatilde/Lambdat
    #w = ww*invLambda;
    
    omegatildeprime=2*a
    #lal version: 435
    #ww_r=2.*a - (a2*a*coeffs->b3*eta)*u2 - coeffs->bb3*eta*a*u2;
    
    Lambdatprime=4*(a*a+r*r)*r-2*a*a*Deltatprime*sin2theta
    #lal version:433
    #Lambda_r = 4.*r*w2 - a2*deltaT_r*xi2;
    #BUG spotted - we're off by a factor of 2
    
    mucostheta=a*a*costheta/Sigma
    #lal version: 444
    #mucos = (a2*costheta)*invrho2;
    
    nucostheta=a*a*w2*costheta*(w2-Deltat)/(Lambdat*Sigma)
    #lal version: 443 
    #nucos = (a2*costheta)*w2*(w2-deltaT)*(invrho2*invLambda);
    
    omegacostheta=-2*a*a*costheta*Deltat*omegatilde/(Lambdat*Lambdat)
    #lal version: 442
    #wcos  = -2.*(a2*costheta)*deltaT*ww*(invLambda*invLambda);
    
    mur=r/Sigma-1/np.sqrt(Deltar)
    #lal version: 440
    #mur = (r*invrho2 - invsqrtdeltaR);
    
    nur=r/Sigma+w2*(w2*Deltatprime-4*r*Deltat)/(2*Lambdat*Deltat)
    #lal version: 439
    #nur = (r*invrho2 + (w2 * (-4.*r*deltaT + w2*deltaT_r) ) * 0.5*invdeltaT*invLambda );
    
    omegar=(Lambdat*omegatildeprime-Lambdatprime*omegatilde)/(Lambdat*Lambdat)
    #lal version: 438
    #wr = (-Lambda_r*ww + Lambda*ww_r)*(invLambda*invLambda);
    
    dSO=147.481*chi*chi*chi*eta*eta-568.651*chi*chi*chi*eta+66.1987*chi*chi*chi-343.313*chi*chi*eta+2495.29*chi*eta*eta-44.5324
    sigmacoeffTerm3 = eta*dSO*u*u*u
    #lal version: 494-499
    #deltaSigmaStar_x += coeffs->d1 * eta * sigmaStar->data[0] * u3;
    #deltaSigmaStar_y += coeffs->d1 * eta * sigmaStar->data[1] * u3;
    #deltaSigmaStar_z += coeffs->d1 * eta * sigmaStar->data[2] * u3;
    #deltaSigmaStar_x += coeffs->d1v2 * eta * sigmaKerr->data[0] * u3;
    #deltaSigmaStar_y += coeffs->d1v2 * eta * sigmaKerr->data[1] * u3;
    #deltaSigmaStar_z += coeffs->d1v2 * eta * sigmaKerr->data[2] * u3;
    
    sigmacoeffTerm2=eta/(144*r*r)*(-896+r*(-436*Qminus1-96*DrSipn2+r*(-45*Qminus1*Qminus1+36*Qminus1*DrSipn2))+eta*(-336+r*(204*Qminus1-882*DrSipn2+r*(810*DrSipn2*DrSipn2-234*Qminus1*DrSipn2))))
    #lal version: 482-487
    #sMultiplier2 = (-56.0/9.0*u2+(-2.0/3.0*pn2*u2-109.0/36.0*pp*u2
    #                             +(pn2*pp*u2/4.0-5.0/16.0*pp*pp*u2)*r)*r
    #                           +(-7.0/3.0*u2+(-49.0/8.0*pn2*u2+17.0/12.0*pp*u2
    #                                          +(45.0/8.0* pn2*pn2*u2
    #                                            -13.0/8.0*pn2*pp*u2)*r)*r)*eta)
    #              *eta;
    
    sigmacoeffTerm1=eta/12*(-8/r+3*Qminus1-36*DrSipn2)
    #lal version: 458
    #deltaSigmaStar_x=eta*((-8. - 3.*r*(12.*pn2 - pp))*sKerr_x + (14. + (- 30.*pn2 + 4.*pp)*r)*sStar_x)*(1./12.)*u;

    sigmacoeff=sigmacoeffTerm1+sigmacoeffTerm2+sigmacoeffTerm3
    
    sigmastarcoeffTerm2=eta/(72*r*r)*(706+r*(-206*Qminus1+282*DrSipn2+r*Qminus1*(96*DrSipn2-23*Qminus1))+eta*(-54+r*(120*Qminus1-324*DrSipn2+r*(360*DrSipn2*DrSipn2+Qminus1*(-126*DrSipn2-3*Qminus1)))))
    #lal version: 472-475
    #sMultiplier1 = (-706.0+(206.0*pp-282.0*pn2+(-96.0*pn2*pp+23.0*pp*pp)*r)*r
    #               +(54.0+( -120.0*pp+324.0*pn2+(-360.0*pn2*pn2+126.0*pn2*pp
    #                                              +3.0*pp*pp)*r)*r)*eta)*eta*u2
    #              *(-1./72.0);
    
    sigmastarcoeffTerm1=eta/12*(14/r+4*Qminus1-30*DrSipn2)
    #lal version: 458
    #deltaSigmaStar_x = eta*((-8. - 3.*r*(12.*pn2 - pp))*sKerr_x + (14. + (- 30.*pn2 + 4.*pp)*r)*sStar_x)*(1./12.)*u;    
    
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
    #lal version 505-507
    #sx = sStar_x + deltaSigmaStar_x;
    #sy = sStar_y + deltaSigmaStar_y;
    #sz = sStar_z + deltaSigmaStar_z;
    
    Sstardotn=Sstar1*n1+Sstar2*n2+Sstar3*n3
    SdotSkerrhat=S1*Skerrhat1+S2*Skerrhat2+S3*Skerrhat3
    
    Sdotn=S1*n1+S2*n2+S3*n3
    Sdotv=S1*v1+S2*v2+S3*v3
    Sdotxi=S1*xi1+S2*xi2+S3*xi3
    #lal version: 510-512
    #sxi = sx*xi_x + sy*xi_y + sz*xi_z;
    #sv  = sx*vx + sy*vy + sz*vz;
    #sn  = sx*nx + sy*ny + sz*nz;    
    
    HdsumTerm2=3*Sstardotn*Sstardotn
    HdsumTerm1=Sstar1*Sstar1+Sstar2*Sstar2+Sstar3*Sstar3
    Hdsum=HdsumTerm1-HdsumTerm2
    Hdcoeff=np.divide(1,2)/(r*r*r)
    #Hd = Hdcoeff*Hdsum
    #lal version: 532
    #Hss = -0.5*u3 * (sx*sx + sy*sy + sz*sz - 3.*sn*sn);
    
    Q4=2*prT*prT*prT*prT*u*u*(4-3*eta)*eta
    gammappsum=Deltar/Sigma*pdotn*pdotn+1/Sigma*pdotvr*pdotvr/sin2theta+Sigma/Lambdat/sin2theta*pdotxir*pdotxir
    Hnsradicand=1+gammappsum+Q4
    alpha=np.sqrt(Deltat*Sigma/Lambdat)
    betapsum=omegatilde*pphi/Lambdat
    #Hns=betapsum+alpha*np.sqrt(Hnsradicand)
    #lal version: 403-404
    #Hns = sqrt((1. + ((prT*prT)*(prT*prT))*qq*u2 + ptheta2*invrho2 + pf*pf*rho2*invLambda*invxi2 + pr*pr*deltaR*invrho2)
    #          * (rho2*deltaT) * invLambda) + pf*ww*invLambda;
    
    HssTerm3=expmu*expnu*pdotxir*(Jtilde*pdotn*Sdotxi*Btilde-expmu*expnu*pdotxir*Sdotn)+(pdotvr*(Sdotn*pdotvr-Jtilde*pdotn*Sdotv)-exp2mu*(np.sqrt(Q)+Q)*Sdotn*xisq)*Btilde*Btilde
    HssTerm3coeff=omegacostheta/(2*exp2mu*expmu*expnu*Btilde*(Q+np.sqrt(Q)))
    #lal version: 530(term->)
    
    
    HssTerm2=expmu*pdotxir*(expmu*exp2nu*pdotxir*Sdotv-expnu*pdotvr*Sdotxi*Btilde)+xisq*Btilde*Btilde*(exp2mu*(np.sqrt(Q)+Q)*Sdotv+Jtilde*pdotn*(pdotvr*Sdotn-Jtilde*pdotn*Sdotv))
    HssTerm2coeff=Jtilde*omegar/(2*exp2mu*expmu*expnu*Btilde*(Q+np.sqrt(Q))*xisq)
    
    
    HssTerm1=omega*SdotSkerrhat
    #lal version: 530(term1 -> w*s3)
    #Hs = w*s3 + Hwr*wr + Hwcos*wcos + HSOL + HSONL;
    
    Hss=HssTerm1+HssTerm2coeff*HssTerm2+HssTerm3coeff*HssTerm3
    
    
    HsoTerm2c=Jtilde*Brtilde*expmu*expnu*pdotxir*(np.sqrt(Q)+1)*Sdotv
    HsoTerm2b=expmu*expnu*pdotxir*(2*np.sqrt(Q)+1)*(Jtilde*nur*Sdotv-nucostheta*Sdotn*xisq)*Btilde
    HsoTerm2a=Sdotxi*Jtilde*(mur*pdotvr*(np.sqrt(Q)+1)-mucostheta*pdotn*xisq-np.sqrt(Q)*(nur*pdotvr+(mucostheta-nucostheta)*pdotn*xisq))*Btilde*Btilde
    HsoTerm2=HsoTerm2a+HsoTerm2b-HsoTerm2c
    HsoTerm2coeff=expnu/(exp2mu*Btilde*Btilde*(Q+np.sqrt(Q))*xisq)
    #HSONL(lal) = HsoTerm2*HsoTerm2coeff
    #lal version:526 - 528
    #HSONL = ((expnu*(invexpMU*invexpMU))*(-(B*expMU*expnu*nucos*pxir*(1. + 2.*sqrtQ)*sn*xi2) +
    #     (-(BR*(expMU*expnu)*pxir*(1. + sqrtQ)*sv) + B*((expMU*expnu)*nur*pxir*(1. + 2.*sqrtQ)*sv + B*mur*pvr*sxi +
    #     B*sxi*(-(mucos*pn*xi2) + sqrtQ*(mur*pvr - nur*pvr + (-mucos + nucos)*pn*xi2))))*sqrtdeltaR))*invxi2/(deltaT*(sqrtQ + Q));
    
    HsoTerm1=exp2nu*(expmu*expnu-Btilde)*pdotxir*SdotSkerrhat/(expmu*Btilde*Btilde*np.sqrt(Q)*xisq)
    #lal version: 524 
    #HSOL = ((expnu*expnu*invexpMU)*(-B + (expMU*expnu))*pxir*s3)/(deltaT*sqrtQ)*invxi2;
    
    Hso=HsoTerm1+HsoTerm2coeff*HsoTerm2
    
    Hd=Hdcoeff*Hdsum
    
    Hns=betapsum+alpha*np.sqrt(Hnsradicand)
    
    Hs=Hso+Hss
    
    dSS=528.511*chi*chi*chi*eta*eta-41.0003*chi*chi*chi*eta+1161.78*chi*chi*eta*eta*eta-326.325*chi*chi*eta*eta+37.1964*chi*eta+706.958*eta*eta*eta-36.0272*eta+6.06807
    
    Heff = Hs + Hns - Hd + dSS*eta*u*u*u*u*(S1x*S1x + S1y*S1y + S1z*S1z + S2x*S2x + S2y*S2y + S2z*S2z)
    
    Hreal=np.sqrt(1+2*eta*(Heff-1))
    return [Hreal, csi] 
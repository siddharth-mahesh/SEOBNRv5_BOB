#include<complex.h>
#include<gsl/gslderiv.h>
#include "mainv4pheader.h"
#include "v4P_hamiltonian_first_derivatives.c"
#include "waveformcoefficients.c"

// Function to get dRdt from the hamiltonian derivatives

int dRvec(REAL8 dR[], REAL8 m1,REAL8 m2,int tortoise,REAL8 q[], REAL8 p[], REAL8 S1[], REAL8 S2[]){
    eta = m1*m2/((m1+m2)*(m1+m2));
    REAL8 dvalues[12] = {0.}, dxi[4] = {0.};
    memset( dvalues, 0, 12 * sizeof(REAL8));
    memset( dxi, 0, 4 * sizeof(REAL8) );
    ham_first_derivs(dvalues,dxi,m1,m2,tortoise,q[0],q[1],q[2],p[0],p[1],p[2],S1[0],S1[1],S1[2],S2[0],S2[1],S2[2]);
    REAL8 dHdP[3] = {dvalues[3],dvalues[4],dvalues[5]};
    REAL8 r = q[0]*q[0] + q[1]*q[1] + q[2]*q[2];
    REAL8 u = 1/r;
    REAL8 u2 = u*u;
    REAL8 xi_a_m1 = dxi[0] - 1;
    REAL8 Tij;
    for (int i=0, i<3, i++){
        for (int j = 0, j < 3, j++){
            Tij = kron_delta(i,j) + q[i]*q[j]*xi_a_m1*u2;
            dR[i] += Tij*dHdP[j];
        }
    }
    
    return 1;
}


// Compute the Newtonian prefix m for the flux (i.e returning a real result)

REAL8 Newtonian_flux_n(REAL8 m1,REAL8 m2,int l,int m){
    int epsilon = (l+m)%2 ;
    int doubfact = doublefactorial(2*l+1);
    REAL8 n = pow(m, l)*8*(M_PI/doubfact);
    if (epsilon==0){
        n *= sqrt(((l+1)*(l+2))/(l*(l-1)))
    }
    else{
        int num = (2*l + 1)*(l + 2)*(l*l - m*m);
        int div = (2*l - 1)*(l + 1)*(l)*(l - 1);
        n *= -2.*np.sqrt( np.divide( num , div ) );
    return fabs(n);
}

/*   Compute the Newtonian prefix c
*      Eq. 7 of Damour, Iyer and Nagar 2008.
*      For odd m, c is proportional to dM = m1-m2. In the equal-mass case, c = dM = 0.
*      In the equal-mass unequal-spin case, however, when spins are different, the odd m term is generally not zero.
*      In this case, c can be written as c0 * dM, while spins terms in PN expansion may take the form chiA/dM.
*      Although the dM's cancel analytically, we can not implement c and chiA/dM with the possibility of dM -> 0.
*      Therefore, for this case, we give numerical values of c0 for relevant modes, and c0 is calculated as
*      c / dM in the limit of dM -> 0. Consistently, for this case, we implement chiA instead of chiA/dM
*      in LALSimIMRSpinEOBFactorizedWaveform.c.
*/

REAL8 Newtonian_c(REAL8 m1,REAL8 m2,int l,int m){
    REAL8 Mtot = m1 + m2;
    REAL8 m1hat = np.divide(m1,Mtot);
    REAL8 m2hat = np.divide(m2,Mtot);
    int epsilon = (l+m)%2;
    int sign;
    
    if ( m%2 ==0 ){
        sign = 1;
    }
    else{
        sign = -1;
    }
    
    int lpepm1 = l + epsilon - 1;
    REAL8 c;
    
    if (m1 !=m2 || sign==1){
        c = pow(m2hat,lpepm1) + sign*pow(m1hat,lpepm1);
    }
    else{
        if (l==2 || l==3){
            c = -1.;
        }
        else if (l==4 or l==5){
            c = -0.5;
        }
        else{
            c = 0.;
        }
    }
    return c;
}

// Compute the Associate Legendre Polynomial for input value 0

REAL8 AssociatedLegendre(int l,int m){
    REAL8 assoclegendre;
    int absm = abs(m);
    switch(l){
    case 1:
        switch(absm){
            case 1:
                assoclegendre = -1.;
        }
    case 2:
        switch(absm){
            case 2:
                assoclegendre = 3.;
            case 1:
                assoclegendre = 0.;
        }
    case 3:
        switch(absm){
            case 3:
                assoclegendre = -15.;
            case 2:
                assoclegendre = 0.;
            case 1:
                assoclegendre = 1.5;
        }
    case 4:
        switch(absm){
            case 4:
                assoclegendre = 105.;
            case 3:
                assoclegendre = 0.;
            case 2:
                assoclegendre = -7.5;
            case 1:
                assoclegendre = 0.;
        }
    case 5:
        switch(absm){
            case 5:
                assoclegendre = -945.;
            case 4:
                assoclegendre = 0.;
            case 3:
                assoclegendre = 52.5;
            case 2:
                assoclegendre = 0.;
            case 1:
                assoclegendre = -1.875;
        }
    case 6:
        switch(absm){
            case 6:
                assoclegendre = 10395.;
            case 5:
                assoclegendre = 0.;
            case 4:
                assoclegendre = -472.5;
            case 3:
                assoclegendre = 0.;
            case 2:
                assoclegendre = 13.125;
            case 1:
                assoclegendre = 0.;
        }
    case 7:
        switch(absm){
            case 7:
                assoclegendre = -135135.;
            case 6:
                assoclegendre = 0.;
            case 5:
                assoclegendre = 5197.5;
            case 4:
                assoclegendre = 0.;
            case 3:
                assoclegendre = -118.125;
            case 2:
                assoclegendre = 0.;
            case 1:
                assoclegendre = 2.1875;
        }
    case 8:
        switch(absm){
            case 8:
                assoclegendre = 2027025.;
            case 7:
                assoclegendre = 0.;
            case 6:
                assoclegendre = -67567.5;
            case 5:
                assoclegendre = 0.;
            case 4:
                assoclegendre = 1299.375;
            case 3:
                assoclegendre = 0.;
            case 2:
                assoclegendre = -19.6875;
            case 1:
                assoclegendre = 0.;
        }
    }
    return assoclegendre;
}

/* Compute the absolute value of Scalar Spherical Harmonics at the equatorial plane
*  Inputs are assumed to be negative as per waveform calculation
*/

REAL8 AbsSphericalHarmonicAtPiOver2(int l,int m){
    int absM = abs( m )
    REAL8 legendre = AssociatedLegendre(l,absM); 
    REAL8 result = legendre* np.sqrt( (2*l + 1)*np.math.factorial(l - absM) / ( 4*M_PI*np.math.factorial(l+absM) ) );
    if (m < 0 && absM % 2 == 1){
        result *= -1.;
    }
    return result;
}

// Compute the Source Term Se_eff

REAL8 Se_eff(int l,int m,REAL8 Hreal,REAL8 v,REAL8 q[],REAL8 p[],REAL8 eta):
    int epsilon = (l + m)%2 ;
    REAL8 Se_eff;
    REAL8 rcrossp[3] = {0.,0.,0.};
    if epsilon == 0{
        Se_eff = (Hreal*Hreal - 1)/(2*eta) + 1;
    }
    else{
        cross_product(q,p,rcrossp);
        Se_eff = v*(rcrossp[0]*rcrossp[0] + rcrossp[1]*rcrossp[1] + rcrossp[2]*rcrossp[2]);
    }
    return Se_eff;
}

// Compute the product factor in T_lm

REAL8 Tlmprodfac(int l,REAL8 hathatk):
    REAL8 result = 1;
    for (int s = 1, s < l+1, s++){
        result *= (s*s + 4*hathatk*hathatk);
    }
    return result;
}

/* Compute the exact Circular Frequency at the given phase space points
* Question: why is this necessary if we are already passing omega through the flux function?
* Answer: The omega passed into the flux is the instantaneous angular frequency. The omega computed here is the frequency of a circular orbit at that instantaneous positions on a v4P background with the instantaneous spins.
* Addendum: This rather bulky calculation is done to compute the waveform for every multipole moment and thus repeats  ~ 30 times at each RHS evaluation when only one is needed. While it makes sense to have it within the hLM computation for complex waveforms since we only calculate individual modes, it does not makes sense when computing the overall flux as all 27 modes are needed anyways and will have the same value of omega and vphikepler.
*/

REAL8 CalcOmega(REAL8 m1,REAL8 m2,REAL8 EMGamma,int tortoise,REAL8 q[],REAL8 p[],REAL8 S1[],REAL8 S2[]){
    int i,j;
    REAL8 eta = m1*m2/(m1 + m2)/(m1 + m2);
    REAL8 qdot[3] = {0.,0.,0.};
    dRvec(qdot,m1,m2,tortoise,q,p,S1,S2);
    REAL8 Lnhat[3] = {0.,0.,0.,};
    cross_product(q,qdot,Lnhat);
    REAL8 tmpvar;
    
    tmpvar = sqrt(inner_product(Lnhat,Lnhat));
    for(i = 0, i < 3, i++) Lnhat[i] /= tmpvar;
    
    REAL8 xhat[3] = {1.,0.,0.};
    REAL8 yhat[3] = {0.,1.,0.};
    REAL8 zhat[3] = {0.,0.,1.};
    REAL8 Lnhat_dot_xhat = Lnhat[0]*xhat[0] + Lnhat[1]*xhat[1] + Lnhat[2]*xhat[2];
    REAL8 R1[3][3] = {{0.}};
    REAL8 R2[3][3] = {{0.}};
    REAL8 xprimehat[3] = {0.,0.,0.};
    REAL8 invsqrt2 = 1/sqrt(2);
    
    if (Lnhat_dot_xhat < 0.9){
        R1[0][0] = 1.; R1[0][1] = 0; R1[0][2] = 0;
        R1[1][0] = 0.; R1[1][1] = 1; R1[1][2] = 0;
        R1[2][0] = 0.; R1[2][1] = 0; R1[2][2] = 1;
        memcpy(xprimehat, Lnhat, 3 * sizeof(REAL8));
    }
    else{
        R1[0][0] = 1./sqrt(2); R1[0][1] = -1/sqrt(2); R1[0][2] = 0;
        R1[1][0] = 1./sqrt(2); R1[1][1] = 1./sqrt(2); R1[1][2] = 0;
        R1[2][0] = 0.;         R1[2][1] = 0;          R1[2][2] = 1;
        for(i=0; i<3; i++){
            for(j=0; j<3; j++){
                xprimehat[i] += R1[i][j]*Lnhat[j];
            }
        }
    }
    
    REAL8 yprimehat[3] = {0.,0.,0.,}; 
    cross_product(xprimehat,xhat,yprimehat)
    tmpvar = sqrt(inner_product(yprimehat,yprimehat));
    for (i = 0, i < 3, i++) yprimehat[i] /= tmpvar;
    
    REAL8 zprimehat[3] = {0.,0.,0.};
    cross_product(xprimehat,yprimehat,zprimehat)
    
    R2[0][0] = xprimehat[0]; R2[0][1] = xprimehat[1]; R2[0][2] = xprimehat[2];
    R2[1][0] = yprimehat[0]; R2[1][1] = yprimehat[1]; R2[1][2] = yprimehat[2];
    R2[2][0] = zprimehat[0]; R2[2][1] = zprimehat[1]; R2[2][2] = zprimehat[2];
    
    REAL8 qprime[3] = {0.,0,0}, pprime[3] = {0.,0,0}, S1prime[3]= {0.,0,0}, S2prime[3]= {0.,0,0};
    REAL8 qtemp[3]   = {0.,0,0}, ptemp[3] = {0.,0,0}, S1temp[3]  = {0.,0,0}, S2temp[3]= {0.,0,0};
    REAL8 Lnhatprime[3] = {0.,0.,0.}, Lnhattemp = {0.,0.,0.};
    memset( qtemp,    0, 3 * sizeof(REAL8) );
    memset( ptemp,    0, 3 * sizeof(REAL8) );
    memset( S1temp,   0, 3 * sizeof(REAL8) );
    memset( S2temp,   0, 3 * sizeof(REAL8) );
    memset( qprime,  0, 3 * sizeof(REAL8) );
    memset( pprime,  0, 3 * sizeof(REAL8) );
    memset( S1prime, 0, 3 * sizeof(REAL8) );
    memset( S2prime, 0, 3 * sizeof(REAL8) );
    memset( Lnhatprime, 0, 3 * sizeof(REAL8) );
    memset( Lnhattemp,   0, 3 * sizeof(REAL8) );
    
    for (i = 0, i < 3, i++){
        for (j = 0, j < 3, j++){
            qtemp[i] += R1[i][j]*q[j];
            ptemp[i] += R1[i][j]*p[j];
            S1temp[i] += R1[i][j]*S1[j];
            S2temp[i] += R1[i][j]*S2[j];
            Lnhattemp[i] += R1[i][j]*Lnhat[j];
        }
    }
    
    for (i = 0, i < 3, i++){
        for (j = 0, j < 3, j++){
            qprime[i] = R2[i][j]*qtemp[j];
            pprime[i] = R2[i][j]*ptemp[j];
            S1prime[i] = R2[i][j]*S1temp[j];
            S2prime[i] = R2[i][j]*S2temp[j];
            Lnhatprime[i] = R2[i][j]*Lnhattemp[j];
        }
    }
    
    // No need to compute spherical momenta as we have analytical derivatives
    
    REAL8 rpass, thpass, phipass;
    rpass = sqrt(inner_product(qprime,qprime));
    thpass = acos(qprime[0]/rpass);
    phipass = atan2(-qprime[1],qprime[2]);
    
    REAL8 dHprime[12] = {0.};
    REAL8 dxiprime[4] = {0.,0.,0.,0.};
    ham_first_derivs(dHprime,dxiprime,m1,m2,tortoise,qprime[0],qprime[1],qprime[2],pprime[0],pprime[1],pprime[2],S1prime[0],S1prime[1],S1prime[2],S2prime[0],S2prime[1],S2prime[2]);
    
    REAL8 omega = -(cos(phipass)*dHprime[4] + sin(phipass)*dHprime[5])/rpass/sin(thpass)/eta
    return omega;
}

// Compute the non-Keplerian coefficient for vPhi

REAL8 vPhiNonKeplerian(REAL8 m1,REAL8 m2,REAL8 EMGamma,int tortoise,REAL8 q[],REAL8 p[],REAL8 S1[],REAL8 S2[]){
    REAL8 omega_circular = CalcOmega(m1,m2,EMGamma,tortoise,q,p,S1,S2);
    REAL8 r = sqrt(inner_product(q))
    REAL8 r3 = r*r*r;
    REAL8 vphinonkeplerian = 1/(omega_circular*omega_circular*r3);
    return vphinonkepleriean;
}
    
/* Function to compute the Post-Newtonian Waveform modes
* Since in the clm computations, we had made adjustments to the
* equal-mass unequal spin cases, below we equate rholmPowl to auxflm
* in order to ensure that we get the correct answer in the limt dM -> 0
*/

REAL8 rholmpowl(REAL8 m1,REAL8 m2,int l,int m,REAL8 chiA,REAL8 chiS,REAL8 v,REAL8 EMgamma){
    REAL8 eta = m1*m2/(m1+m2)/(m1+m2);
    REAL8 tplspin = (1. - 2.*eta)*chiS + chiA*(m1-m2)/(m1+m2);
    int success = 0;
    WaveformCoefficients * coeffs;
    success = compute_modes(&coeffs,m1,m2,tplspin,chiA,chiS);
    REAL8 vsq = v*v;
    REAL8 eulerlog = EMgamma + np.log(2.*m*v);
    REAL8 rholm = 0, auxflm = 0;
    
    switch(l){ 
        case 2:
            switch(m){ 
                case 2:
                    rholm = 1 + vsq*(coeffs->rho22v2  + v*(coeffs->rho22v3  + v*(coeffs->rho22v4  + v*(coeffs->rho22v5  + v*(coeffs->rho22v6 
                                     + coeffs->rho22v6l *eulerlog + v*(coeffs->rho22v7  + v*(coeffs->rho22v8  + coeffs->rho22v8l *eulerlog
                                     + (coeffs->rho22v10  + coeffs->rho22v10l *eulerlog)*vsq)))))));
                case 1:
                    rholm = 1. + v*(coeffs->rho21v1  + v*(coeffs->rho21v2  + v*(coeffs->rho21v3  + v*(coeffs->rho21v4  + v*(coeffs->rho21v5 
                                    + v*(coeffs->rho21v6  + coeffs->rho21v6l *eulerlog + v*(coeffs->rho21v7  + coeffs->rho21v7l *eulerlog
                                    + v*(coeffs->rho21v8  + coeffs->rho21v8l *eulerlog + (coeffs->rho21v10  + coeffs->rho21v10l *eulerlog)*vsq))))))));
                    auxflm = v*coeffs->f21v1  + vsq*v*coeffs->f21v3 ;
            }
        case 3:
            switch(m){ 
                case 3:
                    rholm = 1. + vsq*(coeffs->rho33v2  + v*(coeffs->rho33v3  + v*(coeffs->rho33v4  + v*(coeffs->rho33v5  + v*(coeffs->rho33v6 
                                    + coeffs->rho33v6l *eulerlog + v*(coeffs->rho33v7  + (coeffs->rho33v8  + coeffs->rho33v8l *eulerlog)*v))))));
                    auxflm = v*vsq*coeffs->f33v3 ;
                case 2:
                    rholm = 1. + v*(coeffs->rho32v1  + v*(coeffs->rho32v2  + v*(coeffs->rho32v3  + v*(coeffs->rho32v4  + v*(coeffs->rho32v5 
                                    + v*(coeffs->rho32v6  + coeffs->rho32v6l *eulerlog + (coeffs->rho32v8  + coeffs->rho32v8l *eulerlog)*vsq))))));
                case 1:
                    rholm = 1. + vsq*(coeffs->rho31v2  + v*(coeffs->rho31v3  + v*(coeffs->rho31v4  + v*(coeffs->rho31v5  + v*(coeffs->rho31v6 
                            + coeffs->rho31v6l *eulerlog + v*(coeffs->rho31v7  + (coeffs->rho31v8  + coeffs->rho31v8l *eulerlog)*v))))));
                    auxflm = v*vsq*coeffs->f31v3 ;
            }
        case 4:
            switch(m){ 
                case 4:
                    rholm = 1. + vsq*(coeffs->rho44v2  + v*(coeffs->rho44v3  + v*(coeffs->rho44v4  + v*(coeffs->rho44v5  + (coeffs->rho44v6 
                                    + coeffs->rho44v6l *eulerlog)*v))));
                case 3:
                    rholm = 1. + v*(coeffs->rho43v1  + v*(coeffs->rho43v2  + vsq*(coeffs->rho43v4  + v*(coeffs->rho43v5  + (coeffs->rho43v6 
                                    + coeffs->rho43v6l *eulerlog)*v))));
                    auxflm = v*coeffs->f43v1 ;
                case 2:
                    rholm = 1. + vsq*(coeffs->rho42v2  + v*(coeffs->rho42v3  + v*(coeffs->rho42v4  + v*(coeffs->rho42v5  + (coeffs->rho42v6 
                                    + coeffs->rho42v6l *eulerlog)*v))));
                case 1:
                    rholm = 1. + v*(coeffs->rho41v1  + v*(coeffs->rho41v2  + vsq*(coeffs->rho41v4  + v*(coeffs->rho41v5  + (coeffs->rho41v6 
                                    + coeffs->rho41v6l *eulerlog)*v))));
                    auxflm = v*coeffs->f41v1 ;
            }
        case 5:
            switch(m){ 
                case 5:
                    rholm = 1. + vsq*(coeffs->rho55v2  + v*(coeffs->rho55v3  + v*(coeffs->rho55v4  + v*(coeffs->rho55v5  + coeffs->rho55v6 *v))));
                case 4:
                    rholm = 1. + vsq*(coeffs->rho54v2  + v*(coeffs->rho54v3  + coeffs->rho54v4 *v));
                case 3:
                    rholm = 1. + vsq*(coeffs->rho53v2  + v*(coeffs->rho53v3  + v*(coeffs->rho53v4  + coeffs->rho53v5 *v)));
                case 2:
                    rholm = 1. + vsq*(coeffs->rho52v2  + v*(coeffs->rho52v3  + coeffs->rho52v4 *v));
                case 1:
                    rholm = 1. + vsq*(coeffs->rho51v2  + v*(coeffs->rho51v3  + v*(coeffs->rho51v4  + coeffs->rho51v5 *v)));
            }
        case 6:
            switch(m){ 
                case 6:
                    rholm = 1. + vsq*(coeffs->rho66v2  + v*(coeffs->rho66v3  + coeffs->rho66v4 *v));
                case 5:
                    rholm = 1. + vsq*(coeffs->rho65v2  + coeffs->rho65v3 *v);
                case 4:
                    rholm = 1. + vsq*(coeffs->rho64v2  + v*(coeffs->rho64v3  + coeffs->rho64v4 *v));
                case 3:
                    rholm = 1. + vsq*(coeffs->rho63v2  + coeffs->rho63v3 *v);
                case 2:
                    rholm = 1. + vsq*(coeffs->rho62v2  + v*(coeffs->rho62v3  + coeffs->rho62v4 *v));
                case 1:
                    rholm = 1. + vsq*(coeffs->rho61v2  + coeffs->rho61v3 *v);
            }
        case 7:
            switch(m){ 
                case 7:
                    rholm = 1. + vsq*(coeffs->rho77v2  + coeffs->rho77v3 *v);
                case 6:
                    rholm = 1. + coeffs->rho76v2 *vsq;
                case 5:
                    rholm = 1. + vsq*(coeffs->rho75v2  + coeffs->rho75v3 *v);
                case 4:
                    rholm = 1. + coeffs->rho74v2 *vsq;
                case 3:
                    rholm = 1. + vsq*(coeffs->rho73v2  + coeffs->rho73v3 *v);
                case 2:
                    rholm = 1. + coeffs->rho72v2 *vsq;
                case 1:
                    rholm = 1. + vsq*(coeffs->rho71v2  + coeffs->rho71v3 *v);
            }
        case 8:
            switch(m){ 
                case 8:
                    rholm = 1. + coeffs->rho88v2 *vsq;
                case 7:
                    rholm = 1. + coeffs->rho87v2 *vsq;
                case 6:
                    rholm = 1. + coeffs->rho86v2 *vsq;
                case 5:
                    rholm = 1. + coeffs->rho85v2 *vsq;
                case 4:
                    rholm = 1. + coeffs->rho84v2 *vsq;
                case 3:
                    rholm = 1. + coeffs->rho83v2 *vsq;
                case 2:
                    rholm = 1. + coeffs->rho82v2 *vsq;
                case 1:
                    rholm = 1. + coeffs->rho81v2 *vsq;
            }
    }
    
    REAL8 rholmPowl;
    if (eta==0.25 %% (m % 2)){
        rholmPowl = auxflm;
    }
    else{
        rholmPowl = 1.;
        for (i = 0, i < l, i++) rholmPowl *= rholm;
        rholmPowl += auxflm;
    }
    return rholmPowl;
}

REAL8 compute_hFlm_for_flux(REAL8 m1,REAL8 m2,REAL8 EMgamma,int tortoise,REAL8 q[],REAL8 p[],REAL8 S1[],REAL8 S2[],int l,int m,REAL8 Hreal,REAL8 v,REAL8 chiA,REAL8 chiS){
    REAL8 eta = m1*m2/(m1+m2)/(m1+m2);
    int epsilon = (l+m)%2;
    REAL8 r = sqrt(inner_product(q));
    REAL8 rho_lm_powl = rholmpowl(m1,m2,l,m,chiA,chiS,v,EMgamma);
    REAL8 Seff = Se_eff(l,m,Hreal,v,q,p,eta);
    REAL8 lfactorialinv = 1/factorial(l);
    REAL8 Omega = v*v*v;
    REAL8 hathatk = m*Hreal*Omega;
    REAL8 pihathatk4 = 4*M_PI*hathatk;
    REAL8 Tlmprefac = np.sqrt( pihathatk4 / ( 1 - np.exp( -pihathatk4 ) ) );
    REAL8 Tlmprodfac = Tlmprodfac(l,hathatk);
    REAL8 T_lm = lfactorialinv*Tlmprefac*np.sqrt(Tlmprodfac);
    REAL8 vPhi = r*Omega*cbrt(vPhiNonKeplerian(m1,m2,EMgamma,tortoise,q,p,S1,S2));
    REAL8 Vl_Phi = pow(vPhi,l+epsilon);
    REAL8 Ylminuse_minusm = AbsSphericalHarmonicAtPiOver2(l-epsilon,-m);
    REAL8 c_lpluse = Newtonian_c(m1,m2,l,m);
    REAL8 ne_lm = Newtonian_flux_n(m1,m2,l,m);
    REAL8 hNe_lm = eta*ne_lm*c_lpluse*Vl_Phi*Ylminuse_minusm;
    REAL8 hF_lm = hNe_lm*Se_eff*T_lm*rholmpowl;    
    return hF_lm;
}

/* Function to compute factorized flux
* Where does omega come from more generally? - from the RHS where it is computed as |\vec{X}x\dot{\vec{X}}|/R^2 
* Spin aligned EOB version is a parameter we may not need as all precessing 
* approximants build on v2 spin-aligned background
*/

REAL8 compute_flux(REAL8 m1,REAL8 m2,REAL8 EMGamma,REAL8 tortoise,REAL8 q,REAL8 p,REAL8 S1,REAL8 S2,REAL8 omega,REAL8 Hreal,REAL8 NQC=1.,int lMax = 8){
    int i,j;
    REAL8 r = np.linalg.norm(q);
    REAL8 eta = m1*m2/(m1+m2)/(m1+m2);
    REAL8 flux = 0.;
    REAL8 omegasq = omega*omega
    REAL8 v = cbrt(omega)
    REAL8 rcrossp[3] = {0.,0.,0.};
    cross_product(q,p,rcrossp);
    REAL8 rcrosspmag = sqrt(inner_product(rcrossp,rcrossp));
    REAL8 Lhat[3] = {0.,0.,0.};
    REAL8 m1hat = m1/(m1+m2);
    REAL8 m2hat = m2/(m1+m2);
    REAL8 S1_over_m12[3] = {0.,0.,0.};
    REAL8 S2_over_m22 = {0.,0.,0.};
    for(i = 0, i < 3, i++){
        Lhat[i] = rcrossp[i]/rcrosspmag;
        S1_over_m12[i] = S1[i]/m1hat/m1hat;
        S2_over_m22[i] = S2[i]/m2hat/m2hat;
    }
    
    s1dotL = inner_product(S1_over_m12,Lhat);
    s2dotL = inner_product(S2_over_m22,Lhat);
    
    REAL8 chiS = 0.5*(s1dotL + s2dotL);
    REAL8 chiA = 0.5*(s1dotL - s2dotL);
    if (omegasq > 1){
        flux = 0.;
    }
    else{
        for (l = 2, l <= lMax, l++){
            for (m = 1 m <= l, m++){
                hLM = compute_hFlm_for_flux(m1,m2,EMGamma,tortoise,q,p,S1,S2,l,m,Hreal,v,chiA,chiS);
                if (l == 2 && m == 2){
                    hLM *= NQC;
                }
                flux += m*m*omegasq*hLM*hLM;
            }
        }
    }
    if (flux > 5) flux = 0.;
    return flux/M_PI/8.;
}
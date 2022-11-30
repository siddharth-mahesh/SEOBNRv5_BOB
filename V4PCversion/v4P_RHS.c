#include<complex.h>
#include<gsl/gslderiv.h>
#include "mainv4pheader.h"
#include "v4P_hamiltonian_first_derivatives.c"
#include "waveformcoefficients.c"
#include "flux.c"

/* question for developer:
* in the appendix of Pan2010, the flux term is multiplied by pstar but in lal it is multiplied by p. Why so?
*/

int SEOBNR_RHS(REAL8 UNUSED t, const REAL8 values[], REAL8 dvalues[], REAL8 m1, REAL8 m2, REAL8 EMGamma, int tortoise){
    int i,j,k,L;
    REAL8 eta = m1*m2/(m1+m2)/(m1+m2);
    REAL8 q[3], p[3], S1[3], S2[3];
    REAL8 dHreal[12];
    REAL8 dxi[4];
    REAL8 Hreal;
    memset(q,0,3 * sizeof(REAL8));
    memset(p,0,3 * sizeof(REAL8));
    memset(S1,0,3 * sizeof(REAL8));
    memset(S2,0,3 * sizeof(REAL8));
    memset(dHreal,0,12 * sizeof(REAL8));
    memset(dxi,0,4 * sizeof(REAL8));
    
    q[0] = values[0]; q[1] = values[1]; q[2] = values[2];
    p[0] = values[3]; p[1] = values[4]; p[2] = values[5];
    S1[0] = values[6]; S1[1] = values[7]; S1[2] = values[8];
    S2[0] = values[9]; S2[1] = values[10]; S2[2] = values[11];
    
    Hreal = compute_v4P_Hreal(m1,m2,EMGamma,tortoise,values[0],values[1],values[2],values[3],values[4],values[5],values[6],values[7],values[8],values[9],values[10],values[11]);
    
    ham_first_derivs(dHreal,dxi,m1,m2,tortoise,values[0],values[1],values[2],values[3],values[4],values[5],values[6],values[7],values[8],values[9],values[10],values[11]);
    //dHdX = np.array([dHreal[0],dHreal[1],dHreal[2]])/eta
    //dHdP = np.array([dHreal[3],dHreal[4],dHreal[5]])/eta
    //dHdS1 = np.array([dHreal[6],dHreal[7],dHreal[8]])
    //dHdS2 = np.array([dHreal[9],dHreal[10],dHreal[11]])
    
    r = sqrt(inner_product(q,q));
    u = 1/r;
    u2 = u*u;
    
    xi_a_minusone = dxi[0] - 1;
    xi_a_inv = 1/dxi[0];
    
    //dxi_a = np.array([dHreal[12],dHreal[13],dHreal[14]])
    
    REAL8 T[3][3] = {{0.}};
    REAL8 Tinv[3][3] = {{0.}};
    REAL8 dT[3][3][3] = {{{0.}}};
    
    for (i = 0, i < 3, i++){
        for (j = 0, j < 3, j++){
            T[i][j] = kron_delta(i,j) + q[i]*q[j]*xi_a_minusone*u2;
            Tinv[i][j] = kron_delta(i,j) - q[i]*q[j]*xi_a_minusone*u2*xi_a_inv;
        }
    }        
    
    for (i = 0, i < 3, i++){
        for (j = 0, j < 3, j++){
            for (k = 0, k < 3, k++){
                dT[i][j][k] = ((kron_delta(i,k)*q[j] + q[i]*kron_delta(j,k) - 2.*q[i]*q[j]*q[k]*u2)*xi_a_minusone + q[i]*q[j]*dxi[k+1])*u2;
            }
        }
    }
    
    memset(dvalues, 0, 12 * sizeof(REAL8));
    REAL8 dx[3] = {0.,0.,0.};
    
    for (i = 0, i < 3, i++){ 
        for (j = 0, j < 3, j++){
            dx[i] += T[i][j] * dHreal[i + 3];
        }
        dvalues[i] = dx[i];
    }
    
    REAL8 q_cross_dx[3];
    cross_product(q,dx,q_cross_dx);
    REAL8 Omega = sqrt(inner_product(q_cross_dx,q_cross_dx))*u2;
    REAL8 Lvec;
    cross_product(q,p,Lvec);
    REAL8 Lmag = sqrt(inner_product(Lvec,Lvec));
    REAL8 dEdt = compute_flux(m1,m2,EMGamma,tortoise,q,p,S1,S2,Omega,Hreal);
    REAL8 dHdS1[3];
    REAL8 dS1[3];
    REAL8 dHdS2[3];
    REAL8 dS2[3];
    
    dHdS1[0] = dHreal[6]; dHdS1[1] = dHreal[7]; dHdS1[2] = dHreal[8];
    cross_product(dHdS1,S1,dS1);
    dHdS2[0] = dHreal[9]; dHdS2[1] = dHreal[10]; dHdS2[2] = dHreal[11];
    cross_product(dHdS2,S2,dS2);
    
    dvalues[6] = dS1[0]; dvalues[7] = dS1[1]; dvalues[8] = dS1[2]; 
    dvalues[9] = dS2[0]; dvalues[10] = dS2[1]; dvalues[11] = dS2[2];
    
    REAL8 Ham_deriv_term_1[3] = {0.,0.,0.};
    
    for (i = 0, i < 3, i++){
        Ham_deriv_term_1[i] = 0.; 
        for (j = 0, j < 3, j++){
            Ham_deriv_term_1[i] += T[i][j]*(-1.*dHreal[j]);
        }
    }
    
    REAL8 Ham_deriv_term_2[3] = {0.,0.,0.};
    REAL8 Ham_deriv_term_2_T11[3][3][3] = {{{0.}}};
    REAL8 Ham_deriv_term_2_T12[3][3] = {{0.}}
    REAL8 Ham_deriv_term_2_T2[3] = {0.,0.,0.};
    
    for (i = 0, < 3, i++){
        for (j = 0, < 3, j++){
            for (l = 0, l < 3, l++){
                for (k = 0, k < 3, k++){
                    Ham_deriv_term_2_T11[i][j][l] += dT[i][k][j]*Tinv[k][l];
                }
            }
        }
    }
    
    for (i = 0, < 3, i++){
        for (j = 0, < 3, j++){
            for (k = 0, k < 3, k++){
                Ham_deriv_term_2_T12[i][j] += Ham_deriv_term_2_T11[i][j][k]*p[k];
            }
        }
    }
    
    for (i = 0, < 3, i++){
        for (j = 0, j < 3, j++){
            Ham_deriv_term_2_T2[i] += dHdP[j]*T[i][j];
        }
    }
    
    for (i = 0, i < 3, i++){
        for (j = 0, j < 3, j++){
            Ham_deriv_term_2[i] += Ham_deriv_term_2_T12[i][j]*Ham_deriv_term_2_T2[j]
        }
    }
    
    REAL8 pstar[3] = {0.,0.,0.};
    
    for (i = 0, i < 3, i++){
        for (j = 0, j < 3, j++){
            pstar[i] += T[i][j]*p[j];
        }
    }
    
    REAL8 flux_term[3];
    for (i = 0, i < 3, i++) flux_term[i] =  -p[i]*dEdt/Omega/Lmag/eta;
    
    for (i = 0, i < 3, i++) dvalues[3+i] = Ham_deriv_term_1[i] + Ham_deriv_term_2[i] + flux_term[i];
    
    return 1;
}
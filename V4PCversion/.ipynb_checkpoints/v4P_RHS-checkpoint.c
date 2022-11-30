#include<complex.h>
#include<gsl/gslderiv.h>
#include "mainv4pheader.h"
#include "v4P_hamiltonian_first_derivatives.c"
#include "waveformcoefficients.c"
#include "flux.c"

/* question for developer:
* in the appendix of Pan2010, the flux term is multiplied by pstar but in lal it is *multiplied by p. Why so?
*/

int SEOBNR_RHS(REAL8 UNUSED t, const REAL8 values[], REAL8 dvalues[], REAL8 m1, REAL8 m2, REAL8 EMGamma, int tortoise){
    
    
    eta = m1*m2/(m1+m2)/(m1+m2)
    q = np.array([var[0], var[1], var[2]])
    p = np.array([var[3], var[4], var[5]])
    S1 = np.array([var[6], var[7], var[8]])
    S2 = np.array([var[9], var[10], var[11]])
    
    ham_and_potentials = hr.compute_Hreal_and_csi(m1,m2,EMGamma,tortoise,var[0],var[1],var[2],var[3],var[4],var[5],var[6],var[7],var[8],var[9],var[10],var[11])
    Hreal = ham_and_potentials[0]
    
    dHreal = dh.ham_first_derivs(m1,m2,tortoise,var[0],var[1],var[2],var[3],var[4],var[5],var[6],var[7],var[8],var[9],var[10],var[11])
    dHdX = np.array([dHreal[0],dHreal[1],dHreal[2]])/eta
    dHdP = np.array([dHreal[3],dHreal[4],dHreal[5]])/eta
    dHdS1 = np.array([dHreal[6],dHreal[7],dHreal[8]])
    dHdS2 = np.array([dHreal[9],dHreal[10],dHreal[11]])
    
    R = np.linalg.norm(q)
    Rinv = 1/R
    R2inv = Rinv*Rinv
    
    xi_a = ham_and_potentials[1]
    xi_a_minusone = xi_a - 1
    xi_a_inv = 1/xi_a
    dxi_a = np.array([dHreal[12],dHreal[13],dHreal[14]])
    
    kron_delta = lambda x,y: 1 if x== y else 0 
    
    T = np.zeros([3,3])
    Tinv = np.zeros([3,3])
    dT = np.zeros([3,3,3])
    
    for i in range(3):
        for j in range(3):
            T[i,j] = kron_delta(i,j) + q[i]*q[j]*xi_a_minusone*R2inv
            Tinv[i,j] = kron_delta(i,j) - q[i]*q[j]*xi_a_minusone*R2inv*xi_a_inv
    
    for i in range(3):
        for j in range(3):
            for k in range(3):
                dT[i,j,k] = ((kron_delta(i,k)*q[j] + q[i]*kron_delta(j,k) - 2*q[i]*q[j]*q[k]*R2inv)*xi_a_minusone + q[i]*q[j]*dxi_a[k])*R2inv 
    
    dX = np.dot(T,dHdP)
    
    Omega = np.linalg.norm(np.cross(q,dX))*R2inv
    Lvec = np.cross(q,p)
    L = np.linalg.norm(Lvec)
    dEdt = fl.compute_flux(m1,m2,EMGamma,tortoise,q,p,S1,S2,Omega,Hreal)
    dS1 = np.cross(dHdS1,S1)
    dS2 = np.cross(dHdS2,S2)
    Ham_deriv_term_1 = np.zeros(3)
    for i in range(3):
        Ham_deriv_term_1[i] = 0 
        for j in range(3):
            Ham_deriv_term_1[i] += T[i,j]*(-dHdX[j])
    
    Ham_deriv_term_2 = np.zeros(3)
    Ham_deriv_term_2_T11 = np.zeros([3,3,3])
    Ham_deriv_term_2_T12 = np.zeros([3,3])
    Ham_deriv_term_2_T2 = np.zeros([3])
    
    for i in range(3):
        for j in range(3):
            for l in range(3):
                for k in range(3):
                    Ham_deriv_term_2_T11[i,j,l] += dT[i,k,j]*Tinv[k,l]
    
    for i in range(3):
        for j in range(3):
            for k in range(3):
                Ham_deriv_term_2_T12[i,j] += Ham_deriv_term_2_T11[i,j,k]*p[k]
    
    for i in range(3):
        for j in range(3):
            Ham_deriv_term_2_T2[i] += dHdP[j]*T[i,j]
    
    for i in range(3):
        for j in range(3):
            Ham_deriv_term_2[i] += Ham_deriv_term_2_T12[i,j]*Ham_deriv_term_2_T2[j]
                
    pstar = np.zeros(3)
    
    for i in range(3):
        for j in range(3):
            pstar[i] += T[i,j]*p[j]
    
    flux_term = -p*dEdt/Omega/L/eta
    
    dP = Ham_deriv_term_1 + Ham_deriv_term_2 + flux_term
    
    return np.array([dX[0],dX[1],dX[2],dP[0],dP[1],dP[2],dS1[0],dS1[1],dS1[2],dS2[0],dS2[1],dS2[2]])
}
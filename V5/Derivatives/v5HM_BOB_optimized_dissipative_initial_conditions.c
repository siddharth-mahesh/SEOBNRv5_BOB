#include "SEOBNRBOB.h"
double  v5HM_BOB_optimized_dissipative_initial_conditions(double x, void *params){
double m1 = ((struct ic_diss_params *) params)->m1;
double m2 = ((struct ic_diss_params *) params)->m2;
double chi1 = ((struct ic_diss_params *) params)->chi1;
double chi2 = ((struct ic_diss_params *) params)->chi2;
double a6 = ((struct ic_diss_params *) params)->a6;
double dSO = ((struct ic_diss_params *) params)->dSO;
double r = ((struct ic_diss_params *) params)->r;
double pphi = ((struct ic_diss_params *) params)->pphi;
double prstar = x;
double dvalues[4] = {0.,0.,0.,0.};
double ddvalues[2] = {0.,0.};
double hamiltonian_and_xi[2] = {0.,0.};
double flux, omega;
int gsl_status;
gsl_status = v5HM_BOB_optimized_hamiltonian_first_derivatives(dvalues,hamiltonian_and_xi,m1,m2,r,prstar,pphi,chi1,chi2,a6,dSO);
gsl_status = v5HM_BOB_optimized_hamiltonian_second_derivatives(ddvalues,m1,m2,r,prstar,pphi,chi1,chi2,a6,dSO);
gsl_status = v5HM_BOB_optimized_omega(&omega,m1,m2,r,pphi,chi1,chi2,a6,dSO);
gsl_status = v5HM_BOB_optimized_flux(&flux,m1,m2,r,prstar,chi1,chi2,hamiltonian_and_xi[0],dvalues[3],omega);
double drdt = xi*dvalues[2];
double dLdr_inv = ddvalues[1]/ddvalues[0];
double prstareqn = drdt + dLdr_inv*flux;
return prstareqn;
}

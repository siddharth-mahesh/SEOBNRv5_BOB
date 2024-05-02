#include "SEOBNRBOB.h"
int v5HM_BOB_optimized_merger_ringdown(double *h, double *phi, double t, double t0, double hNR, double omegaNR, double omegaQNM, double tau){
const double tmp0 = ((omegaNR)*(omegaNR));
const double tmp1 = (1.0/(omegaQNM));
const double tmp5 = (1.0/16.0)*((omegaNR)*(omegaNR)*(omegaNR)*(omegaNR));
const double tmp16 = omegaQNM/omegaNR;
const double tmp2 = omegaNR*tmp1;
const double tmp3 = 2*log(tmp2);
const double tmp12 = (1.0/2.0)*tmp0*tmp1;
const double tmp13 = tau*tmp0*tmp1;
const double tmp4 = (t - t0 + tau*tmp3)/tau;
const double tmp6 = tanh(tmp3);
const double tmp7 = tmp5 + ((1.0/16.0)*((omegaQNM)*(omegaQNM)*(omegaQNM)*(omegaQNM)) - tmp5)*(-tmp6 + tanh(tmp4))/(1 - tmp6);
const double tmp10 = pow(tmp7, 0.25);
const double tmp15 = 2*tmp1*tmp10;
const double tmp17 = 2*omegaQNM*tmp10/tmp0;
*h = (1.0/4.0)*hNR*tmp0*(1.0/sqrt(tmp7))*cosh(tmp3)/cosh(tmp4);
*phi = omegaQNM*tau*(-atan2((1.0/2.0)*omegaNR, (1.0/2.0)*omegaQNM) + atan2(tmp10, (1.0/2.0)*omegaQNM)) + 0.5*omegaQNM*tau*log((1 - tmp2)*(tmp15 + 1)/((1 - tmp15)*(tmp2 + 1))) - tmp13*(-atan2((1.0/2.0)*omegaNR, tmp12) + atan2(tmp10, tmp12)) - 0.5*tmp13*log((1 - tmp16)*(tmp17 + 1)/((1 - tmp17)*(tmp16 + 1)));
return GSL_SUCCESS;
}

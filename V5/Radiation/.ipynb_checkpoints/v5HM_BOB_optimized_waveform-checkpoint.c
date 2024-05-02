#include "SEOBNRBOB.h"
int v5HM_BOB_optimized_waveform(double *h22,double m1, double m2, double r, double phi, double prstar, double pphi, double chi1, double chi2, double Hreal, double Omega, double Omega_circ){
const double tmp0 = ((Omega)*(Omega));
const double tmp1 = m1 + m2;
const double tmp5 = pow(Omega, 8.0/3.0);
const double tmp6 = log(4*cbrt(Omega));
const double tmp7 = pow(Omega, 10.0/3.0);
const double tmp14 = pow(Omega, 5.0/3.0);
const double tmp20 = ((M_PI)*(M_PI));
const double tmp23 = Hreal*Omega;
const double tmp10 = (1.0/2.0)*chi1 + (1.0/2.0)*chi2;
const double tmp11 = (1.0/2.0)*chi1 - 1.0/2.0*chi2;
const double tmp12 = (m1 - m2)/tmp1;
const double tmp15 = ((m1)*(m1))*((m2)*(m2))/((tmp1)*(tmp1)*(tmp1)*(tmp1));
const double tmp21 = ((m1)*(m1)*(m1))*((m2)*(m2)*(m2))/pow(tmp1, 6);
const double tmp24 = 4*I*tmp23;
const double tmp3 = (1.0/((tmp1)*(tmp1)));
const double tmp13 = tmp11*tmp12;
const double tmp16 = ((tmp10)*(tmp10));
const double tmp17 = ((tmp11)*(tmp11));
const double tmp4 = m1*m2*tmp3;
const double EulerGamma = 0.57721566490153286;
*h22 = -1.0092530088080638*I*M_PI*tmp0*tmp4*(1 + (1.0/2.0)*((tmp1)*(tmp1))*(((Hreal)*(Hreal)) - 1)/(m1*m2))*((pow(Omega, 7.0/3.0)*(((tmp10)*(tmp10)*(tmp10))*(tmp4 + 1.0/3.0) + tmp10*tmp17*(-4*tmp15 - 3*tmp4 + 1) + tmp10*(-245717.0/63504.0*tmp15 + (50803.0/63504.0)*tmp21 + (74749.0/5292.0)*tmp4 + 18733.0/15876.0) + ((tmp11)*(tmp11)*(tmp11))*tmp12*(1.0/3.0 - 4.0/3.0*tmp4) + tmp13*tmp16*(2*tmp4 + 1) + tmp13*((97865.0/63504.0)*tmp15 + (50140.0/3969.0)*tmp4 + 18733.0/15876.0)) + pow(Omega, 4.0/3.0)*(tmp10*tmp13 + (19583.0/42336.0)*tmp15 + (1.0/2.0)*tmp16 + tmp17*(1.0/2.0 - 2*tmp4) - 33025.0/21168.0*tmp4 - 20555.0/10584.0) + pow(Omega, 2.0/3.0)*((55.0/84.0)*tmp4 - 43.0/42.0) + Omega*(-0.33333333333333331*chi1 - 0.33333333333333331*chi2 + 0.66666666666666663*m1*m2*tmp10*tmp3 - 0.66666666666666663*tmp13) + tmp0*(tmp10*tmp13*(89.0/126.0 - 781.0/252.0*tmp4) - 6292061.0/3259872.0*tmp15 + tmp16*((10.0/9.0)*tmp15 - 1817.0/504.0*tmp4 + 89.0/252.0) + tmp17*(-27.0/14.0*tmp15 - 457.0/504.0*tmp4 + 89.0/252.0) + (41.0/192.0)*tmp20*tmp4 + (10620745.0/39118464.0)*tmp21 - 48993925.0/9779616.0*tmp4 - 428.0/105.0*tmp6 - 428.0/105.0*EulerGamma + 1556919113.0/122245200.0) + tmp14*(tmp10*((209.0/126.0)*tmp15 + (49.0/18.0)*tmp4 - 34.0/21.0) + tmp13*(-19.0/42.0*tmp4 - 34.0/21.0)) + tmp4*(21.199999999999999*tmp5 - 411*tmp7) + tmp5*((9202.0/2205.0)*tmp6 - 387216563023.0/160190110080.0 + (9202.0/2205.0)*EulerGamma) + tmp7*((439877.0/55566.0)*tmp6 - 16094530514677.0/533967033600.0 + (439877.0/55566.0)*EulerGamma) + 1)*(pow(Omega, 7.0/3.0)*(((tmp10)*(tmp10)*(tmp10))*(tmp4 + 1.0/3.0) + tmp10*tmp17*(-4*tmp15 - 3*tmp4 + 1) + tmp10*(-245717.0/63504.0*tmp15 + (50803.0/63504.0)*tmp21 + (74749.0/5292.0)*tmp4 + 18733.0/15876.0) + ((tmp11)*(tmp11)*(tmp11))*tmp12*(1.0/3.0 - 4.0/3.0*tmp4) + tmp13*tmp16*(2*tmp4 + 1) + tmp13*((97865.0/63504.0)*tmp15 + (50140.0/3969.0)*tmp4 + 18733.0/15876.0)) + pow(Omega, 4.0/3.0)*(tmp10*tmp13 + (19583.0/42336.0)*tmp15 + (1.0/2.0)*tmp16 + tmp17*(1.0/2.0 - 2*tmp4) - 33025.0/21168.0*tmp4 - 20555.0/10584.0) + pow(Omega, 2.0/3.0)*((55.0/84.0)*tmp4 - 43.0/42.0) + Omega*(-0.33333333333333331*chi1 - 0.33333333333333331*chi2 + 0.66666666666666663*m1*m2*tmp10*tmp3 - 0.66666666666666663*tmp13) + tmp0*(tmp10*tmp13*(89.0/126.0 - 781.0/252.0*tmp4) - 6292061.0/3259872.0*tmp15 + tmp16*((10.0/9.0)*tmp15 - 1817.0/504.0*tmp4 + 89.0/252.0) + tmp17*(-27.0/14.0*tmp15 - 457.0/504.0*tmp4 + 89.0/252.0) + (41.0/192.0)*tmp20*tmp4 + (10620745.0/39118464.0)*tmp21 - 48993925.0/9779616.0*tmp4 - 428.0/105.0*tmp6 - 428.0/105.0*EulerGamma + 1556919113.0/122245200.0) + tmp14*(tmp10*((209.0/126.0)*tmp15 + (49.0/18.0)*tmp4 - 34.0/21.0) + tmp13*(-19.0/42.0*tmp4 - 34.0/21.0)) + tmp4*(21.199999999999999*tmp5 - 411*tmp7) + tmp5*((9202.0/2205.0)*tmp6 - 387216563023.0/160190110080.0 + (9202.0/2205.0)*EulerGamma) + tmp7*((439877.0/55566.0)*tmp6 - 16094530514677.0/533967033600.0 + (439877.0/55566.0)*EulerGamma) + 1))*exp(-2*I*phi)*exp(1.0*I*(((Hreal)*(Hreal)*(Hreal))*((Omega)*(Omega)*(Omega))*((1712.0/315.0)*tmp20 - 2203.0/81.0) + ((Hreal)*(Hreal))*tmp0*(tmp10*((8.0/3.0)*tmp4 - 4.0/3.0) - 4.0/3.0*tmp13 + (428.0/105.0)*M_PI) - 24*tmp14*tmp4 + (7.0/3.0)*tmp23))*exp(2*M_PI*tmp23)*exp(tmp24*log(8*Omega*exp(-1.0/2.0)))*tgamma(3 - tmp24)/pow(Omega_circ, 4.0/3.0);
return GSL_SUCCESS
}
/* 
* Header file for sBOB waveform functions 
* Author: Siddharth Mahesh 
* Standard order of input for all functions is as follows: 
* time, QNM frequency, QNM damping time, reference frequency, peak strain 
* This module was generated using Nrpy+ 
*/
#include <stdio.h> 
 #include<math.h> 
 #include<stdlib.h> 
 #include<string.h>
#ifndef SBOB_FUNCS_H
#define SBOB_FUNCS_H
void get_sBOB_strainamplitude(const double t,const double Omega_QNM,const double tau,const double Omega_0,const double h_peak, double *strain_amplitude);
void get_sBOB_phase(const double t,const double Omega_QNM,const double tau,const double Omega_0,const double h_peak, double *phase);
#endif
/*
Output sBOB strain amplitude
 */
void get_sBOB_strainamplitude(const double t,const double Omega_QNM,const double tau,const double Omega_0,const double h_peak, double *strain_amplitude) {

   /*
    *  Original SymPy expression:
    *  "*strain_amplitude = h_peak*sqrt(Omega_0**4/2 + Omega_QNM**4/2 + (-Omega_0**4/2 + Omega_QNM**4/2)*tanh(log(Omega_0/Omega_QNM)))*cosh(log(Omega_0/Omega_QNM))/(sqrt(Omega_0**4/2 + Omega_QNM**4/2 + (-Omega_0**4/2 + Omega_QNM**4/2)*tanh((t + tau*log(Omega_0/Omega_QNM))/tau))*cosh((t + tau*log(Omega_0/Omega_QNM))/tau))"
    */
      const double tmp_0 = log(Omega_0/Omega_QNM);
      const double tmp_1 = (t + tau*tmp_0)/tau;
      const double tmp_2 = (1.0/2.0)*((Omega_0)*(Omega_0)*(Omega_0)*(Omega_0));
      const double tmp_3 = (1.0/2.0)*((Omega_QNM)*(Omega_QNM)*(Omega_QNM)*(Omega_QNM));
      const double tmp_4 = -tmp_2 + tmp_3;
      *strain_amplitude = h_peak*sqrt(tmp_2 + tmp_3 + tmp_4*tanh(tmp_0))*cosh(tmp_0)/(sqrt(tmp_2 + tmp_3 + tmp_4*tanh(tmp_1))*cosh(tmp_1));
}
/*
Output sBOB waveform phase
 */
void get_sBOB_phase(const double t,const double Omega_QNM,const double tau,const double Omega_0,const double h_peak, double *phase) {

   /*
    *  Original SymPy expression:
    *  "*phase = -2*Omega_0*tau*(-atan((Omega_0**4/2 + Omega_QNM**4/2)**(1/4)/Omega_0) + atan((Omega_0**4/2 + Omega_QNM**4/2 + (-Omega_0**4/2 + Omega_QNM**4/2)*tanh((t + tau*log(Omega_0/Omega_QNM))/tau))**(1/4)/Omega_0)) - Omega_0*tau*log((1 - (Omega_0**4/2 + Omega_QNM**4/2)**(1/4)/Omega_0)*(1 + (Omega_0**4/2 + Omega_QNM**4/2 + (-Omega_0**4/2 + Omega_QNM**4/2)*tanh((t + tau*log(Omega_0/Omega_QNM))/tau))**(1/4)/Omega_0)/((1 + (Omega_0**4/2 + Omega_QNM**4/2)**(1/4)/Omega_0)*(1 - (Omega_0**4/2 + Omega_QNM**4/2 + (-Omega_0**4/2 + Omega_QNM**4/2)*tanh((t + tau*log(Omega_0/Omega_QNM))/tau))**(1/4)/Omega_0))) + 2*Omega_QNM*tau*(-atan((Omega_0**4/2 + Omega_QNM**4/2)**(1/4)/Omega_QNM) + atan((Omega_0**4/2 + Omega_QNM**4/2 + (-Omega_0**4/2 + Omega_QNM**4/2)*tanh((t + tau*log(Omega_0/Omega_QNM))/tau))**(1/4)/Omega_QNM)) + Omega_QNM*tau*log((1 - (Omega_0**4/2 + Omega_QNM**4/2)**(1/4)/Omega_QNM)*(1 + (Omega_0**4/2 + Omega_QNM**4/2 + (-Omega_0**4/2 + Omega_QNM**4/2)*tanh((t + tau*log(Omega_0/Omega_QNM))/tau))**(1/4)/Omega_QNM)/((1 + (Omega_0**4/2 + Omega_QNM**4/2)**(1/4)/Omega_QNM)*(1 - (Omega_0**4/2 + Omega_QNM**4/2 + (-Omega_0**4/2 + Omega_QNM**4/2)*tanh((t + tau*log(Omega_0/Omega_QNM))/tau))**(1/4)/Omega_QNM)))"
    */
      const double tmp_0 = (1.0/(Omega_0));
      const double tmp_1 = (1.0/2.0)*((Omega_0)*(Omega_0)*(Omega_0)*(Omega_0));
      const double tmp_2 = (1.0/2.0)*((Omega_QNM)*(Omega_QNM)*(Omega_QNM)*(Omega_QNM));
      const double tmp_4 = pow(tmp_1 + tmp_2, 1.0/4.0);
      const double tmp_5 = tmp_0*tmp_4;
      const double tmp_6 = (1.0/(Omega_QNM));
      const double tmp_7 = pow(tmp_1 + tmp_2 + (-tmp_1 + tmp_2)*tanh((t + tau*log(Omega_0*tmp_6))/tau), 1.0/4.0);
      const double tmp_8 = tmp_0*tmp_7;
      const double tmp_10 = tmp_4*tmp_6;
      const double tmp_11 = tmp_6*tmp_7;
      *phase = -2*Omega_0*tau*(-atan(tmp_5) + atan(tmp_8)) - Omega_0*tau*log((-tmp_5 + 1)*(tmp_8 + 1)/((tmp_5 + 1)*(-tmp_8 + 1))) + 2*Omega_QNM*tau*(-atan(tmp_10) + atan(tmp_11)) + Omega_QNM*tau*log((-tmp_10 + 1)*(tmp_11 + 1)/((tmp_10 + 1)*(-tmp_11 + 1)));
}

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
void get_sBOB_strainamplitude(const double t, double *strain_amplitude);
void get_sBOB_phase(const double t, double *phase);
#endif
void get_sBOB_strainamplitude(const double t, double *strain_amplitude){
    double *strain_amplitude = (1.0/4.0)*pow((1.0/2.0)*Mf + 1.0/2.0, 2)/(pow(Mf, 2)*sqrt(1.3163592146527707e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + (-1.3163592146527707e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + 0.16906095035683122*pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4))*tanh((1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))*(Mf*(1.4187000000000001*pow(1 - af, -0.499) + 0.69999999999999996)*log(2*Mf*(0.097592999999999999 - 0.091933000000000001*af)/((1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))*(pow(af, 2) - 2.4228000000000001*af + 1.4366000000000001)))/(1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001)) + t)/(Mf*(1.4187000000000001*pow(1 - af, -0.499) + 0.69999999999999996))) + 0.16906095035683122*pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4))*cosh((1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))*(Mf*(1.4187000000000001*pow(1 - af, -0.499) + 0.69999999999999996)*log(2*Mf*(0.097592999999999999 - 0.091933000000000001*af)/((1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))*(pow(af, 2) - 2.4228000000000001*af + 1.4366000000000001)))/(1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001)) + t)/(Mf*(1.4187000000000001*pow(1 - af, -0.499) + 0.69999999999999996))));
}
void get_sBOB_phase(const double t, double *phase){
    double *phase = -2*Mf*(0.097592999999999999 - 0.091933000000000001*af)*(1.4187000000000001*pow(1 - af, -0.499) + 0.69999999999999996)*(-atan(0.64122556145171994*pow(7.7862996266989856e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4), 1.0/4.0)*(pow(af, 2) - 2.4228000000000001*af + 1.4366000000000001)/(0.097592999999999999 - 0.091933000000000001*af)) + atan((pow(af, 2) - 2.4228000000000001*af + 1.4366000000000001)*pow(1.3163592146527707e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + (-1.3163592146527707e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + 0.16906095035683122*pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4))*tanh((1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))*(Mf*(1.4187000000000001*pow(1 - af, -0.499) + 0.69999999999999996)*log(2*Mf*(0.097592999999999999 - 0.091933000000000001*af)/((1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))*(pow(af, 2) - 2.4228000000000001*af + 1.4366000000000001)))/(1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001)) + t)/(Mf*(1.4187000000000001*pow(1 - af, -0.499) + 0.69999999999999996))) + 0.16906095035683122*pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4), 1.0/4.0)/(0.097592999999999999 - 0.091933000000000001*af)))/((1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))*(pow(af, 2) - 2.4228000000000001*af + 1.4366000000000001)) - Mf*(0.097592999999999999 - 0.091933000000000001*af)*(1.4187000000000001*pow(1 - af, -0.499) + 0.69999999999999996)*log((1 - 0.64122556145171994*pow(7.7862996266989856e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4), 1.0/4.0)*(pow(af, 2) - 2.4228000000000001*af + 1.4366000000000001)/(0.097592999999999999 - 0.091933000000000001*af))*(1 + (pow(af, 2) - 2.4228000000000001*af + 1.4366000000000001)*pow(1.3163592146527707e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + (-1.3163592146527707e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + 0.16906095035683122*pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4))*tanh((1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))*(Mf*(1.4187000000000001*pow(1 - af, -0.499) + 0.69999999999999996)*log(2*Mf*(0.097592999999999999 - 0.091933000000000001*af)/((1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))*(pow(af, 2) - 2.4228000000000001*af + 1.4366000000000001)))/(1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001)) + t)/(Mf*(1.4187000000000001*pow(1 - af, -0.499) + 0.69999999999999996))) + 0.16906095035683122*pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4), 1.0/4.0)/(0.097592999999999999 - 0.091933000000000001*af))/((1 + 0.64122556145171994*pow(7.7862996266989856e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4), 1.0/4.0)*(pow(af, 2) - 2.4228000000000001*af + 1.4366000000000001)/(0.097592999999999999 - 0.091933000000000001*af))*(1 - (pow(af, 2) - 2.4228000000000001*af + 1.4366000000000001)*pow(1.3163592146527707e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + (-1.3163592146527707e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + 0.16906095035683122*pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4))*tanh((1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))*(Mf*(1.4187000000000001*pow(1 - af, -0.499) + 0.69999999999999996)*log(2*Mf*(0.097592999999999999 - 0.091933000000000001*af)/((1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))*(pow(af, 2) - 2.4228000000000001*af + 1.4366000000000001)))/(1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001)) + t)/(Mf*(1.4187000000000001*pow(1 - af, -0.499) + 0.69999999999999996))) + 0.16906095035683122*pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4), 1.0/4.0)/(0.097592999999999999 - 0.091933000000000001*af))))/((1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))*(pow(af, 2) - 2.4228000000000001*af + 1.4366000000000001)) + 2*(0.35467500000000002*pow(1 - af, -0.499) + 0.17499999999999999)*log((-1.2824511229034399*Mf*pow(7.7862996266989856e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4), 1.0/4.0)/(1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001)) + 1)*(2*Mf*pow(1.3163592146527707e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + (-1.3163592146527707e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + 0.16906095035683122*pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4))*tanh((1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))*(Mf*(1.4187000000000001*pow(1 - af, -0.499) + 0.69999999999999996)*log(2*Mf*(0.097592999999999999 - 0.091933000000000001*af)/((1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))*(pow(af, 2) - 2.4228000000000001*af + 1.4366000000000001)))/(1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001)) + t)/(Mf*(1.4187000000000001*pow(1 - af, -0.499) + 0.69999999999999996))) + 0.16906095035683122*pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4), 1.0/4.0)/(1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001)) + 1)/((1.2824511229034399*Mf*pow(7.7862996266989856e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4), 1.0/4.0)/(1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001)) + 1)*(-2*Mf*pow(1.3163592146527707e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + (-1.3163592146527707e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + 0.16906095035683122*pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4))*tanh((1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))*(Mf*(1.4187000000000001*pow(1 - af, -0.499) + 0.69999999999999996)*log(2*Mf*(0.097592999999999999 - 0.091933000000000001*af)/((1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))*(pow(af, 2) - 2.4228000000000001*af + 1.4366000000000001)))/(1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001)) + t)/(Mf*(1.4187000000000001*pow(1 - af, -0.499) + 0.69999999999999996))) + 0.16906095035683122*pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4), 1.0/4.0)/(1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001)) + 1))) + 2*(0.70935000000000004*pow(1 - af, -0.499) + 0.34999999999999998)*(-atan(1.2824511229034399*Mf*pow(7.7862996266989856e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4), 1.0/4.0)/(1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))) + atan(2*Mf*pow(1.3163592146527707e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + (-1.3163592146527707e-6*pow(1 - 0.94200403717479742*af, 4)/pow(0.41274558362225522*pow(af, 2) - af + 0.59295030543173188, 4) + 0.16906095035683122*pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4))*tanh((1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))*(Mf*(1.4187000000000001*pow(1 - af, -0.499) + 0.69999999999999996)*log(2*Mf*(0.097592999999999999 - 0.091933000000000001*af)/((1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))*(pow(af, 2) - 2.4228000000000001*af + 1.4366000000000001)))/(1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001)) + t)/(Mf*(1.4187000000000001*pow(1 - af, -0.499) + 0.69999999999999996))) + 0.16906095035683122*pow(1 - 0.75850763884335459*pow(1 - af, 0.12920000000000001), 4)/pow(Mf, 4), 1.0/4.0)/(1.5250999999999999 - 1.1568000000000001*pow(1 - af, 0.12920000000000001))));
}
#include "SEOBNRBOB.h"
int  (const gsl_vector * u, void *params, gsl_vector * f){
double m1 = ((struct ic_cons_params *) params)->m1;
double m2 = ((struct ic_cons_params *) params)->m2;
double chi1 = ((struct ic_cons_params *) params)->chi1;
double chi2 = ((struct ic_cons_params *) params)->chi2;
double a6 = ((struct ic_cons_params *) params)->a6;
double dSO = ((struct ic_cons_params *) params)->dSO;
double Omega = ((struct ic_cons_params *) params)->Omega;
double r = gsl_vector_get (u, 0);
double pphi = gsl_vector_get (u, 1);
const double tmp0 = ((r)*(r)*(r));
const double tmp1 = m1 + m2;
const double tmp16 = ((r)*(r));
const double tmp20 = m1 - m2;
const double tmp25 = (1.0/(r));
const double tmp29 = ((pphi)*(pphi));
const double tmp31 = ((m1)*(m1));
const double tmp32 = ((m2)*(m2));
const double tmp35 = ((r)*(r)*(r)*(r));
const double tmp37 = ((pphi)*(pphi)*(pphi)*(pphi));
const double tmp57 = ((r)*(r)*(r)*(r)*(r));
const double tmp62 = log(r);
const double tmp64 = a6*(3.0*r + 59.0);
const double tmp65 = ((m1)*(m1)*(m1));
const double tmp66 = ((m2)*(m2)*(m2));
const double tmp69 = 8.0*r;
const double tmp94 = 53760.0*r;
const double tmp95 = 907881.73955968802*r;
const double tmp97 = 1140531.4845898801*r;
const double tmp105 = 53760.0*a6 + 20550748.608964998;
const double tmp107 = 2.22557561555965e-7*r;
const double tmp2 = (1.0/((tmp1)*(tmp1)));
const double tmp5 = chi1*m1 + chi2*m2;
const double tmp10 = (1.0/(tmp0));
const double tmp13 = (1.0/((tmp1)*(tmp1)*(tmp1)));
const double tmp17 = (1.0/(tmp16));
const double tmp19 = (1.0/((tmp1)*(tmp1)*(tmp1)*(tmp1)));
const double tmp21 = chi1*m1 - chi2*m2;
const double tmp36 = (1.0/(tmp35));
const double tmp44 = 2*tmp25 + 1;
const double tmp58 = (1.0/(tmp57));
const double tmp67 = pow(tmp1, -6);
const double tmp70 = 2.0*tmp0 + 4.0*tmp16 + tmp69 + 16.0;
const double tmp74 = ((m1)*(m1)*(m1)*(m1))*((m2)*(m2)*(m2)*(m2))/pow(tmp1, 8);
const double tmp93 = 1692004.4940084701*r + 545881.12350211805*tmp0 + 1091762.2470042401*tmp16 - 877100.56175105902*tmp35 - 4316471.01198305;
const double tmp96 = 226970.43488992201*tmp0 + 453940.86977984401*tmp16 - 113485.217444961*tmp35 + tmp95 + 1815763.47911938;
const double tmp98 = 285132.87114747101*tmp0 + 570265.74229494203*tmp16 + 51667.379039702697*tmp35 + tmp97 + 5388804.0029947804;
const double tmp6 = ((tmp5)*(tmp5));
const double tmp11 = pphi*tmp10;
const double tmp18 = pphi*tmp17;
const double tmp22 = tmp20*tmp21;
const double tmp26 = m1*m2*tmp2;
const double tmp30 = tmp17*tmp29;
const double tmp33 = tmp19*tmp31*tmp32;
const double tmp34 = tmp10*tmp29;
const double tmp39 = tmp5/tmp1;
const double tmp55 = tmp2*((tmp21)*(tmp21));
const double tmp68 = tmp65*tmp66*tmp67;
const double tmp71 = a6*(tmp35 + tmp70);
const double tmp108 = ((pphi)*(pphi)*(pphi))*tmp36;
const double tmp7 = tmp2*tmp6;
const double tmp15 = dSO*m1*m2*tmp13*tmp5;
const double tmp23 = -tmp13*((tmp5)*(tmp5)*(tmp5)) + tmp19*tmp22*tmp6;
const double tmp41 = tmp2*tmp22;
const double tmp56 = tmp13*tmp22*tmp5;
const double tmp61 = -21.0/16.0*tmp19*tmp31*tmp32 + (81.0/64.0)*tmp26;
const double tmp63 = 756.0*tmp26;
const double tmp73 = tmp26*tmp62;
const double tmp75 = tmp33*((tmp62)*(tmp62));
const double tmp81 = -45.0/16.0*tmp26 - 0.9375;
const double tmp83 = (75.0/32.0)*tmp26 + (345.0/64.0)*tmp33 + 1.640625;
const double tmp85 = (9.0/16.0)*tmp26 - 0.9375;
const double tmp91 = 688128.0*r + 1548288.0*tmp26 + 833536.0;
const double tmp99 = tmp16*(1128960.0*tmp26 + 2071680.0);
const double tmp100 = tmp0*(725760.0*tmp26 + 1035840.0);
const double tmp101 = tmp35*(362880.0*tmp26 + 517920.0);
const double tmp102 = r*(-1791904.9465588899*tmp26 - 725760.0*tmp33 + 4143360.0);
const double tmp106 = r*(6572428.80109422*m1*m2*tmp2 - 1300341.64885297*tmp33 - 1720320.0);
const double tmp8 = tmp0 + tmp7*(r + 2);
const double tmp42 = tmp39*(tmp17*((109.0/192.0)*tmp19*tmp31*tmp32 - 177.0/32.0*tmp26 - 5.0/64.0) + tmp25*((23.0/32.0)*tmp26 - 3.0/32.0) + tmp30*(-45.0/32.0*tmp26 - 15.0/32.0) + tmp34*(-267.0/128.0*tmp26 - 1591.0/768.0*tmp33 + 59.0/256.0) + tmp36*tmp37*((75.0/128.0)*tmp26 + (345.0/256.0)*tmp33 + 105.0/256.0) + 7.0/4.0) + tmp41*(tmp17*(-1.0/32.0*tmp26 + (103.0/192.0)*tmp33 + 5.0/64.0) + tmp25*((11.0/32.0)*tmp26 + 3.0/32.0) + tmp30*(15.0/32.0 - 9.0/32.0*tmp26) + tmp34*(-35.0/128.0*tmp26 - 613.0/768.0*tmp33 - 59.0/256.0) + tmp36*tmp37*((75.0/256.0)*tmp19*tmp31*tmp32 - 45.0/128.0*tmp26 - 105.0/256.0) + 1.0/4.0);
const double tmp72 = tmp63 + 1079.0;
const double tmp76 = 13212057600.0*tmp57 + 241555486248.80701*tmp74 + 67645734912.0*tmp75;
const double tmp86 = (75.0/64.0)*tmp19*tmp31*tmp32 - 45.0/32.0*tmp26 - 1.640625;
const double tmp103 = 32768.0*tmp100 + 32768.0*tmp101 + 32768.0*tmp102 - 1272782568716.72*tmp26 - 61684325962.210297*tmp33 + 5284823040.0*tmp57 + 32768.0*tmp99 + 440653578240.0;
const double tmp9 = (1.0/(tmp8));
const double tmp43 = pphi*tmp42 + tmp11*tmp15 + (1.0/4.0)*tmp18*tmp23;
const double tmp46 = tmp16 + tmp44*tmp7;
const double tmp51 = tmp17*tmp44*tmp7 + 1;
const double tmp77 = tmp17*tmp7 + 7680.0*tmp35*(2048.0*m1*m2*tmp2*tmp62*(336.0*r + tmp63 + 407.0) + 28.0*m1*m2*tmp2*(1920.0*a6 + 733955.30746303697) - 7.0*r*(938918.40015631705*m1*m2*tmp2 - 185763.09269328101*tmp33 - 245760.0) - 5416406.5954118604*tmp33 - 3440640.0)/(53760.0*tmp26*(113485.217444961*r*(-tmp35 + tmp70) + 148.04406601634*r*(7704.0*r + 1926.0*tmp0 + 3852.0*tmp16 + 349.0*tmp35 + 36400.0) + 128.0*r*(13218.7851094412*r + 4264.6962773603*tmp0 + 8529.3925547206109*tmp16 - 6852.3481386801504*tmp35 - 33722.4297811176) + 7680.0*tmp71) + 7.0*tmp33*(745857848.11560404*a6 + 122635399361.987*r - 6178501407.5375795*tmp0 + 2064783811.32587*tmp16 - 3089250703.7687898*tmp35 + 1426660551.8843999*tmp57 - 39321600.0*tmp64 + 276057889687.01099) + 1120.0*tmp68*(-163683964.82255101*r - 17833256.898554999*tmp16 - 1188987459.03162) + 32768.0*tmp73*(240.0*r*(-7466.2706106620599*tmp26 - 3024.0*tmp33 + 17264.0) + 960.0*tmp0*tmp72 + 1920.0*tmp16*(588.0*tmp26 + 1079.0) - 38842241.476950698*tmp26 - 1882456.2366397199*tmp33 + 480.0*tmp35*tmp72 + 161280.0*tmp57 + 13447680.0) + tmp76) + tmp36*(tmp55*((1.0/2.0)*tmp26 + 1.0/8.0) - 5.0/4.0*tmp56 + (9.0/8.0)*tmp7) + tmp58*(tmp55*(-tmp61 - 9.0/64.0) + tmp56*(117.0/32.0 - 39.0/16.0*tmp26) + tmp7*(-175.0/64.0*tmp26 - 225.0/64.0));
const double tmp104 = (1.0/(tmp103*tmp73 + tmp26*(412876800.0*tmp71 + tmp93*tmp94 + tmp94*tmp96 + tmp94*tmp98) + tmp33*(5221004936.8092299*a6 + 858447795533.90906*r - 43249509852.7631*tmp0 + 14453486679.281099*tmp16 - 21624754926.3815*tmp35 + 9986623863.1907997*tmp57 - 275251200.0*tmp64 + 1932405227809.0801) + tmp68*(-183326040601.25699*r - 19973247726.381599*tmp16 - 1331665954115.4099) + tmp76));
const double tmp47 = (1.0/(tmp46));
const double tmp52 = (1.0/(tmp51));
const double tmp49 = -tmp44*tmp47;
const double tmp50 = tmp30*tmp49*tmp7 + tmp30 + 1;
const double tmp79 = sqrt(tmp50*tmp52*tmp77);
const double tmp80 = (1.0/sqrt(tmp26*(2*tmp43*tmp9 + 2*tmp79 - 2) + 1));
const double Hreal_prmr = tmp80*(tmp43*(-3*tmp16 - tmp7)/((tmp8)*(tmp8)) + tmp9*(-3*pphi*tmp15*tmp36 + pphi*(tmp39*(-tmp10*((109.0/96.0)*tmp19*tmp31*tmp32 - 177.0/16.0*tmp26 - 0.15625) - tmp17*((23.0/32.0)*tmp26 - 0.09375) - tmp29*tmp36*(-801.0/128.0*tmp26 - 1591.0/256.0*tmp33 + 0.69140625) - tmp34*tmp81 - tmp37*tmp58*tmp83) + tmp41*(tmp10*tmp29*tmp85 - tmp10*(-1.0/16.0*tmp26 + (103.0/96.0)*tmp33 + 0.15625) - tmp17*((11.0/32.0)*tmp26 + 0.09375) - tmp29*tmp36*(-105.0/128.0*tmp26 - 613.0/256.0*tmp33 - 0.69140625) - tmp37*tmp58*tmp86)) - 1.0/2.0*tmp11*tmp23) + tmp51*tmp79*((1.0/2.0)*tmp50*(tmp52*(tmp0*tmp104*(30720.0*m1*m2*tmp105*tmp2 + 30720.0*m1*m2*tmp2*tmp62*tmp91 - 30720.0*tmp106 - 166392010611.052*tmp33 - 105696460800.0) - 2*tmp10*tmp7 + tmp104*tmp35*(7680.0*tmp25*tmp26*tmp91 - 50476253192.403603*tmp26 + 9986623863.1907902*tmp33 + 5284823040.0*tmp73 + 13212057600.0) + tmp35*(1.3162167359092599e-19*m1*m2*tmp105*tmp2 + 1.3162167359092599e-19*m1*m2*tmp2*tmp62*tmp91 - 1.3162167359092599e-19*tmp106 - 7.1291650093703898e-13*tmp33 - 4.5286279502388399e-13)*(-tmp103*tmp25*tmp26 - 135291469824.0*tmp25*tmp33*tmp62 - tmp26*(412876800.0*a6*(4*tmp0 + 6.0*tmp16 + tmp69 + 8.0) + 201084856528.177*r + 56877242932.044098*tmp0 + 113754485864.088*tmp16 - 50476253192.403603*tmp35 + tmp94*(2183524.4940084699*r - 3508402.2470042398*tmp0 + 1637643.37050636*tmp16 + 1692004.4940084701) + tmp94*(-453940.86977984401*tmp0 + 680911.30466976599*tmp16 + tmp95 + 907881.73955968802) + tmp94*(206669.51615881099*tmp0 + 855398.613442412*tmp16 + tmp97 + 1140531.4845898801) + 155264066234.24799) - tmp33*(-825753600.0*a6 + 28906973358.562199*r - 86499019705.526398*tmp0 - 129748529558.289*tmp16 + 49933119315.954002*tmp35 + 858447795533.90906) - 66060288000.0*tmp35 + tmp65*tmp66*tmp67*(39946495452.763199*r + 183326040601.25699) - tmp73*(32768.0*r*(2257920.0*tmp26 + 4143360.0) + 32768.0*tmp0*(1451520.0*tmp26 + 2071680.0) + 32768.0*tmp16*(2177280.0*tmp26 + 3107520.0) - 58717141288.841904*tmp26 - 23781703680.0*tmp33 + 26424115200.0*tmp35 + 135769620480.0))/((tmp26*(tmp107*tmp93 + tmp107*tmp96 + tmp107*tmp98 + 0.00170924207274981*tmp71) + tmp33*(0.021614102076040202*a6 + 3.5538327399018002*r - 0.17904586032964401*tmp0 + 0.059835058618348702*tmp16 - 0.089522930164821907*tmp35 + 0.041342980936911501*tmp57 - 0.0011394947151665399*tmp64 + 7.9998399449253697) + 0.054695746327994003*tmp57 + tmp68*(-0.75893966826498704*r - 0.082685961873822905*tmp16 - 5.5128781167229297) + tmp73*(1.3565413275792199e-7*tmp100 + 1.3565413275792199e-7*tmp101 + 1.3565413275792199e-7*tmp102 - 5.2691105819295299*tmp26 - 0.25536296823610299*tmp33 + 0.021878298531197701*tmp57 + 1.3565413275792199e-7*tmp99 + 1.8242333680060501) + tmp74 + 0.28004222119932898*tmp75)*(tmp26*(tmp107*tmp93 + tmp107*tmp96 + tmp107*tmp98 + 0.00170924207274981*tmp71) + tmp33*(0.021614102076040202*a6 + 3.5538327399018002*r - 0.17904586032964401*tmp0 + 0.059835058618348702*tmp16 - 0.089522930164821907*tmp35 + 0.041342980936911501*tmp57 - 0.0011394947151665399*tmp64 + 7.9998399449253697) + 0.054695746327994003*tmp57 + tmp68*(-0.75893966826498704*r - 0.082685961873822905*tmp16 - 5.5128781167229297) + tmp73*(1.3565413275792199e-7*tmp100 + 1.3565413275792199e-7*tmp101 + 1.3565413275792199e-7*tmp102 - 5.2691105819295299*tmp26 - 0.25536296823610299*tmp33 + 0.021878298531197701*tmp57 + 1.3565413275792199e-7*tmp99 + 1.8242333680060501) + tmp74 + 0.28004222119932898*tmp75)) + tmp58*(5*tmp13*tmp20*tmp21*tmp5 - 4*tmp55*((1.0/2.0)*tmp26 + 0.125) - 9.0/2.0*tmp7) - (5*tmp55*(-tmp61 - 0.140625) + 5*tmp56*(3.65625 - 39.0/16.0*tmp26) + 5*tmp7*(-175.0/64.0*tmp26 - 3.515625))/pow(r, 6)) + tmp77*(tmp10*tmp7*(4*tmp25 + 2) + 2*tmp36*tmp7)/((tmp51)*(tmp51))) + (1.0/2.0)*tmp52*tmp77*(tmp17*tmp2*tmp29*tmp6*(2*tmp17*tmp47 - tmp44*(-2*r + 2*tmp17*tmp2*tmp6)/((tmp46)*(tmp46))) - 2*tmp34*tmp49*tmp7 - 2*tmp34))/(tmp50*tmp77));
const double Hreal_prmpphi = tmp80*(tmp79*(2*pphi*tmp17*tmp49*tmp7 + 2*pphi*tmp17)/(2*tmp30*tmp49*tmp7 + 2*tmp30 + 2) + tmp9*(pphi*(tmp39*(tmp108*tmp83 + tmp11*(-267.0/64.0*tmp26 - 1591.0/384.0*tmp33 + 0.4609375) + tmp18*tmp81) + tmp41*(tmp108*tmp86 + tmp11*(-35.0/64.0*tmp26 - 613.0/384.0*tmp33 - 0.4609375) - tmp18*tmp85)) + tmp10*tmp15 + (1.0/4.0)*tmp17*tmp23 + tmp42));
gsl_vector_set (f, 0, Hreal_prmr);
gsl_vector_set (f, 1, Omega - Hreal_prmpphi);
return GSL_SUCCESS;
}
#include "SEOBNRBOB.h"
int v5HM_BOB_optimized_hamiltonian_first_derivatives(double dvalues[], double hamiltonian_and_xi[] double m1, double m2, double r, double prstar, double pphi, double chi1, double chi2, double a6, double dSO){
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
const double tmp45 = ((prstar)*(prstar));
const double tmp49 = log(r);
const double tmp57 = ((m1)*(m1)*(m1));
const double tmp58 = ((m2)*(m2)*(m2));
const double tmp80 = ((r)*(r)*(r)*(r)*(r));
const double tmp81 = a6*(3.0*r + 59.0);
const double tmp82 = 8.0*r;
const double tmp110 = ((prstar)*(prstar)*(prstar)*(prstar));
const double tmp116 = pow(prstar, 6);
const double tmp120 = pow(prstar, 8);
const double tmp125 = pow(r, -4.5);
const double tmp128 = pow(r, -2.5);
const double tmp131 = pow(r, -3.5);
const double tmp176 = pow(r, -6);
const double tmp181 = 53760.0*r;
const double tmp182 = 907881.73955968802*r;
const double tmp184 = 1140531.4845898801*r;
const double tmp192 = 53760.0*a6 + 20550748.608964998;
const double tmp194 = 2.22557561555965e-7*r;
const double tmp227 = ((prstar)*(prstar)*(prstar));
const double tmp228 = ((prstar)*(prstar)*(prstar)*(prstar)*(prstar));
const double tmp229 = pow(prstar, 7);
const double tmp2 = (1.0/((tmp1)*(tmp1)));
const double tmp5 = chi1*m1 + chi2*m2;
const double tmp10 = (1.0/(tmp0));
const double tmp13 = (1.0/((tmp1)*(tmp1)*(tmp1)));
const double tmp17 = (1.0/(tmp16));
const double tmp19 = (1.0/((tmp1)*(tmp1)*(tmp1)*(tmp1)));
const double tmp21 = chi1*m1 - chi2*m2;
const double tmp36 = (1.0/(tmp35));
const double tmp55 = ((tmp49)*(tmp49));
const double tmp59 = pow(tmp1, -6);
const double tmp75 = 879923036160.0*tmp0;
const double tmp83 = 2.0*tmp0 + 4.0*tmp16 + tmp82 + 16.0;
const double tmp88 = ((m1)*(m1)*(m1)*(m1))*((m2)*(m2)*(m2)*(m2))/pow(tmp1, 8);
const double tmp134 = 2*tmp25 + 1;
const double tmp149 = (1.0/(tmp80));
const double tmp180 = 1692004.4940084701*r + 545881.12350211805*tmp0 + 1091762.2470042401*tmp16 - 877100.56175105902*tmp35 - 4316471.01198305;
const double tmp183 = 226970.43488992201*tmp0 + 453940.86977984401*tmp16 + tmp182 - 113485.217444961*tmp35 + 1815763.47911938;
const double tmp185 = 285132.87114747101*tmp0 + 570265.74229494203*tmp16 + tmp184 + 51667.379039702697*tmp35 + 5388804.0029947804;
const double tmp6 = ((tmp5)*(tmp5));
const double tmp11 = pphi*tmp10;
const double tmp18 = pphi*tmp17;
const double tmp22 = tmp20*tmp21;
const double tmp26 = m1*m2*tmp2;
const double tmp30 = tmp17*tmp29;
const double tmp33 = tmp19*tmp31*tmp32;
const double tmp34 = tmp10*tmp29;
const double tmp39 = tmp5/tmp1;
const double tmp56 = 5787938193408.0*tmp55;
const double tmp60 = tmp57*tmp58*tmp59;
const double tmp84 = a6*(tmp35 + tmp83);
const double tmp98 = tmp2*((tmp21)*(tmp21));
const double tmp111 = tmp110*tmp36;
const double tmp232 = ((pphi)*(pphi)*(pphi))*tmp36;
const double tmp7 = tmp2*tmp6;
const double tmp15 = dSO*m1*m2*tmp13*tmp5;
const double tmp23 = -tmp13*((tmp5)*(tmp5)*(tmp5)) + tmp19*tmp22*tmp6;
const double tmp41 = tmp2*tmp22;
const double tmp50 = 14700.0*tmp26 + 42911.0;
const double tmp52 = 5041721180160.0*tmp33 - 104186110149937.0;
const double tmp61 = 133772083200.0*tmp33;
const double tmp62 = 1822680546449.21*tmp33;
const double tmp63 = r*tmp26;
const double tmp64 = 2589101062873.8101*tmp26;
const double tmp69 = tmp16*(-630116198.87329912*tmp26 - 197773496.79353401*tmp33 + 5805304367.8791304);
const double tmp70 = tmp0*(-2675575.66847905*tmp26 - 138240.0*tmp33 - 5278341.3229328999);
const double tmp71 = tmp26*(-2510664218.2812872*tmp26 - 42636451.603233099*tmp33 + 14515200.0*tmp60 + 1002013764.01019);
const double tmp72 = ((1 - 0.49694878161693501*tmp26)*(1 - 0.49694878161693501*tmp26));
const double tmp85 = 756.0*tmp26;
const double tmp87 = tmp26*tmp49;
const double tmp99 = tmp13*tmp22*tmp5;
const double tmp103 = (115.0/64.0)*tmp26 + (1.0/16.0)*tmp33;
const double tmp108 = tmp26*(452.54216699669303 - 51.695238095238103*tmp49);
const double tmp109 = tmp33*(118.40000000000001*tmp49 - 1796.13660498019);
const double tmp114 = tmp10*tmp110*(92.711044284954397*tmp26 - 131.0*tmp33 + 10.0*tmp60);
const double tmp115 = tmp110*tmp17*(8.0*m1*m2*tmp2 - 6.0*tmp33);
const double tmp118 = tmp10*tmp116*(-33.978212217043598*tmp26 - 89.529832736123396*tmp33 + 188.0*tmp57*tmp58*tmp59 - 14.0*tmp88);
const double tmp119 = tmp116*tmp17*(-2.7830076369500598*tmp26 - 5.4000000000000004*tmp33 + 6.0*tmp57*tmp58*tmp59);
const double tmp122 = tmp120*tmp17*(1.3897775099612799*tmp26 + 3.3384202364832198*tmp33 + 3.4285714285714302*tmp60 - 6.0*tmp88);
const double tmp123 = tmp25*tmp26;
const double tmp147 = (165.0/32.0)*tmp26 - 5*tmp33;
const double tmp148 = -75.0/32.0*m1*m2*tmp2 + (15.0/8.0)*tmp33;
const double tmp157 = -21.0/16.0*tmp19*tmp31*tmp32 + (81.0/64.0)*tmp26;
const double tmp165 = -45.0/16.0*tmp26 - 0.9375;
const double tmp167 = (75.0/32.0)*tmp26 + (345.0/64.0)*tmp33 + 1.640625;
const double tmp169 = (9.0/16.0)*tmp26 - 0.9375;
const double tmp175 = 4*tmp98;
const double tmp179 = 688128.0*r + 1548288.0*tmp26 + 833536.0;
const double tmp186 = tmp16*(1128960.0*tmp26 + 2071680.0);
const double tmp187 = tmp0*(725760.0*tmp26 + 1035840.0);
const double tmp188 = tmp35*(362880.0*tmp26 + 517920.0);
const double tmp189 = r*(-1791904.9465588899*tmp26 - 725760.0*tmp33 + 4143360.0);
const double tmp193 = r*(6572428.80109422*m1*m2*tmp2 - 1300341.64885297*tmp33 - 1720320.0);
const double tmp209 = 4161798144000.0*tmp26 + 12148770078720.0;
const double tmp212 = tmp26*(-19691250647040.0*tmp26 - 118648618352640.0);
const double tmp8 = tmp0 + tmp7*(r + 2);
const double tmp42 = tmp39*(tmp17*((109.0/192.0)*tmp19*tmp31*tmp32 - 177.0/32.0*tmp26 - 5.0/64.0) + tmp25*((23.0/32.0)*tmp26 - 3.0/32.0) + tmp30*(-45.0/32.0*tmp26 - 15.0/32.0) + tmp34*(-267.0/128.0*tmp26 - 1591.0/768.0*tmp33 + 59.0/256.0) + tmp36*tmp37*((75.0/128.0)*tmp26 + (345.0/256.0)*tmp33 + 105.0/256.0) + 7.0/4.0) + tmp41*(tmp17*(-1.0/32.0*tmp26 + (103.0/192.0)*tmp33 + 5.0/64.0) + tmp25*((11.0/32.0)*tmp26 + 3.0/32.0) + tmp30*(15.0/32.0 - 9.0/32.0*tmp26) + tmp34*(-35.0/128.0*tmp26 - 613.0/768.0*tmp33 - 59.0/256.0) + tmp36*tmp37*((75.0/256.0)*tmp19*tmp31*tmp32 - 45.0/128.0*tmp26 - 105.0/256.0) + 1.0/4.0);
const double tmp46 = tmp17*tmp7;
const double tmp51 = r*tmp50;
const double tmp53 = 879923036160.0*tmp16 + 25392914995744.301*tmp26 + tmp52;
const double tmp65 = r*tmp62 - 39476764256925.602*r + tmp16*tmp61 + tmp16*tmp64 + 5107745331375.71*tmp16 - 12049908701745.199*tmp26 + 80059249540278.203*tmp33 + tmp56 + 6730497718123.0195*tmp60 + 10611661054566.199*tmp63 + 275059053208689.0;
const double tmp73 = 5927865218923.0195*tmp60 + 86618264430493.297*tmp72 + 188440788778196.0;
const double tmp76 = tmp49*(49152.0*r*(409207698.13607502*tmp26 + 102574080.0*tmp33 - 2119671837.3603799) + 283115520.0*tmp16*tmp50 - 1698693120.0*tmp26*(11592.0*tmp26 + 69847.0) + tmp75);
const double tmp86 = tmp85 + 1079.0;
const double tmp90 = 67645734912.0*tmp33*tmp55 + 13212057600.0*tmp80 + 241555486248.80701*tmp88;
const double tmp100 = tmp10*(tmp7*(3*tmp26 + 45.0/16.0) + tmp98*((3.0/4.0)*tmp26 - 3.0/16.0) - 21.0/8.0*tmp99);
const double tmp104 = tmp36*(tmp7*(-1171.0/64.0*tmp26 - 861.0/64.0) + tmp98*(tmp103 - 37.0/64.0) + tmp99*((13.0/16.0)*tmp26 + 449.0/32.0));
const double tmp112 = tmp111*(tmp108 + tmp109 + 602.31854041656402*tmp60);
const double tmp127 = tmp110*tmp125*tmp26;
const double tmp130 = tmp120*tmp128*tmp26;
const double tmp133 = tmp116*tmp131*tmp26;
const double tmp135 = tmp134*tmp7 + tmp16;
const double tmp151 = tmp110*tmp149*(tmp7*(tmp147 + 25.0/32.0) + tmp98*(-tmp148 - 15.0/32.0) + tmp99*((45.0/8.0)*tmp26 - 5.0/16.0));
const double tmp170 = (75.0/64.0)*tmp19*tmp31*tmp32 - 45.0/32.0*tmp26 - 1.640625;
const double tmp190 = 32768.0*tmp186 + 32768.0*tmp187 + 32768.0*tmp188 + 32768.0*tmp189 - 1272782568716.72*tmp26 - 61684325962.210297*tmp33 + 5284823040.0*tmp80 + 440653578240.0;
const double tmp199 = 4*tmp108 + 4*tmp109 + 2409.2741616662602*tmp60;
const double tmp202 = (45.0/8.0)*tmp26 - 0.3125;
const double tmp203 = tmp147 + 0.78125;
const double tmp205 = -tmp148 - 0.46875;
const double tmp211 = tmp16*tmp209;
const double tmp213 = r*(20113376778784.398*tmp26 + tmp52);
const double tmp9 = (1.0/(tmp8));
const double tmp43 = pphi*tmp42 + tmp11*tmp15 + (1.0/4.0)*tmp18*tmp23;
const double tmp47 = tmp46 + 1;
const double tmp54 = tmp49*(-283115520.0*tmp51 - tmp53);
const double tmp74 = r*(43393301259014.688*tmp26 + 43133561885859.406*tmp33 + tmp73);
const double tmp91 = (2048.0*m1*m2*tmp2*tmp49*(336.0*r + tmp85 + 407.0) + 28.0*m1*m2*tmp2*(1920.0*a6 + 733955.30746303697) - 7.0*r*(938918.40015631705*m1*m2*tmp2 - 185763.09269328101*tmp33 - 245760.0) - 5416406.5954118604*tmp33 - 3440640.0)/(53760.0*tmp26*(113485.217444961*r*(-tmp35 + tmp83) + 148.04406601634*r*(7704.0*r + 1926.0*tmp0 + 3852.0*tmp16 + 349.0*tmp35 + 36400.0) + 128.0*r*(13218.7851094412*r + 4264.6962773603*tmp0 + 8529.3925547206109*tmp16 - 6852.3481386801504*tmp35 - 33722.4297811176) + 7680.0*tmp84) + 7.0*tmp33*(745857848.11560404*a6 + 122635399361.987*r - 6178501407.5375795*tmp0 + 2064783811.32587*tmp16 - 3089250703.7687898*tmp35 + 1426660551.8843999*tmp80 - 39321600.0*tmp81 + 276057889687.01099) + 1120.0*tmp60*(-163683964.82255101*r - 17833256.898554999*tmp16 - 1188987459.03162) + 32768.0*tmp87*(240.0*r*(-7466.2706106620599*tmp26 - 3024.0*tmp33 + 17264.0) + 960.0*tmp0*tmp86 + 1920.0*tmp16*(588.0*tmp26 + 1079.0) - 38842241.476950698*tmp26 - 1882456.2366397199*tmp33 + 480.0*tmp35*tmp86 + 161280.0*tmp80 + 13447680.0) + tmp90);
const double tmp136 = (1.0/(tmp135));
const double tmp153 = tmp134*tmp46 + 1;
const double tmp177 = 2*tmp10*tmp7;
const double tmp191 = (1.0/(tmp190*tmp87 + tmp26*(tmp180*tmp181 + tmp181*tmp183 + tmp181*tmp185 + 412876800.0*tmp84) + tmp33*(5221004936.8092299*a6 + 858447795533.90906*r - 43249509852.7631*tmp0 + 14453486679.281099*tmp16 - 21624754926.3815*tmp35 + 9986623863.1907997*tmp80 - 275251200.0*tmp81 + 1932405227809.0801) + tmp60*(-183326040601.25699*r - 19973247726.381599*tmp16 - 1331665954115.4099) + tmp90));
const double tmp210 = -r*tmp209 - tmp53;
const double tmp214 = (1.0/(r*tmp56 + r*(43393301259015.0*tmp26 + 43133561885859.0*tmp33 + tmp73) + tmp0*(tmp61 + tmp64 + 5107745331375.71) - tmp16*(-5807150888816.293*tmp26 - tmp62 + 53501685054374.102) + tmp26*(-138829688614082.19*tmp26 - 2357625227852.3799*tmp33 + 802632499200.0*tmp60 + 55407353094707.5) + tmp49*(tmp211 + tmp212 + tmp213 + tmp75)));
const double tmp218 = tmp175*tmp205 + 4*tmp202*tmp99 + 4*tmp203*tmp7;
const double tmp221 = 2*tmp46;
const double tmp48 = ((tmp47)*(tmp47));
const double tmp66 = -326837426.24148601*tmp51 - tmp54 + tmp65;
const double tmp77 = r*tmp56 - 9216.0*tmp69 - 967680.0*tmp70 + 55296.0*tmp71 + tmp74 + tmp76;
const double tmp138 = -tmp134*tmp136;
const double tmp140 = ((tmp47)*(tmp47)*(tmp47)*(tmp47));
const double tmp142 = ((r*tmp55 - 1.5922768509339455e-9*tmp69 - 1.6718906934806428e-7*tmp70 + 9.5536611056036724e-9*tmp71 + 1.7277309580446458e-13*tmp74 + 1.7277309580446458e-13*tmp76)*(r*tmp55 - 1.5922768509339455e-9*tmp69 - 1.6718906934806428e-7*tmp70 + 9.5536611056036724e-9*tmp71 + 1.7277309580446458e-13*tmp74 + 1.7277309580446458e-13*tmp76));
const double tmp144 = (1.0/((0.0066265062908739489*r*tmp33 - 0.14352105046684044*r + 0.0094128916415248589*tmp16*tmp26 + 0.00048633950287942822*tmp16*tmp33 + 0.018569631763766858*tmp16 - 0.04380844244600398*tmp26 + 0.29106204142837921*tmp33 - 1.18824456940711e-6*tmp51 - 3.6355829351353647e-15*tmp54 + 0.021042529325572337*tmp55 + 0.024469282648975563*tmp60 + 0.038579573843421422*tmp63 + 1)*(0.0066265062908739489*r*tmp33 - 0.14352105046684044*r + 0.0094128916415248589*tmp16*tmp26 + 0.00048633950287942822*tmp16*tmp33 + 0.018569631763766858*tmp16 - 0.04380844244600398*tmp26 + 0.29106204142837921*tmp33 - 1.18824456940711e-6*tmp51 - 3.6355829351353647e-15*tmp54 + 0.021042529325572337*tmp55 + 0.024469282648975563*tmp60 + 0.038579573843421422*tmp63 + 1)));
const double tmp154 = (1.0/(tmp153));
const double tmp196 = tmp0*tmp191*(30720.0*m1*m2*tmp179*tmp2*tmp49 + 30720.0*m1*m2*tmp192*tmp2 - 30720.0*tmp193 - 166392010611.052*tmp33 - 105696460800.0) + tmp191*tmp35*(7680.0*tmp123*tmp179 - 50476253192.403603*tmp26 + 9986623863.1907902*tmp33 + 5284823040.0*tmp87 + 13212057600.0) + tmp35*(1.3162167359092599e-19*m1*m2*tmp179*tmp2*tmp49 + 1.3162167359092599e-19*m1*m2*tmp192*tmp2 - 1.3162167359092599e-19*tmp193 - 7.1291650093703898e-13*tmp33 - 4.5286279502388399e-13)*(-tmp123*tmp190 - 135291469824.0*tmp25*tmp33*tmp49 - tmp26*(412876800.0*a6*(4*tmp0 + 6.0*tmp16 + tmp82 + 8.0) + 201084856528.177*r + 56877242932.044098*tmp0 + 113754485864.088*tmp16 + tmp181*(2183524.4940084699*r - 3508402.2470042398*tmp0 + 1637643.37050636*tmp16 + 1692004.4940084701) + tmp181*(-453940.86977984401*tmp0 + 680911.30466976599*tmp16 + tmp182 + 907881.73955968802) + tmp181*(206669.51615881099*tmp0 + 855398.613442412*tmp16 + tmp184 + 1140531.4845898801) - 50476253192.403603*tmp35 + 155264066234.24799) - tmp33*(-825753600.0*a6 + 28906973358.562199*r - 86499019705.526398*tmp0 - 129748529558.289*tmp16 + 49933119315.954002*tmp35 + 858447795533.90906) - 66060288000.0*tmp35 + tmp57*tmp58*tmp59*(39946495452.763199*r + 183326040601.25699) - tmp87*(32768.0*r*(2257920.0*tmp26 + 4143360.0) + 32768.0*tmp0*(1451520.0*tmp26 + 2071680.0) + 32768.0*tmp16*(2177280.0*tmp26 + 3107520.0) - 58717141288.841904*tmp26 - 23781703680.0*tmp33 + 26424115200.0*tmp35 + 135769620480.0))/((tmp26*(tmp180*tmp194 + tmp183*tmp194 + tmp185*tmp194 + 0.00170924207274981*tmp84) + 0.28004222119932898*tmp33*tmp55 + tmp33*(0.021614102076040202*a6 + 3.5538327399018002*r - 0.17904586032964401*tmp0 + 0.059835058618348702*tmp16 - 0.089522930164821907*tmp35 + 0.041342980936911501*tmp80 - 0.0011394947151665399*tmp81 + 7.9998399449253697) + tmp60*(-0.75893966826498704*r - 0.082685961873822905*tmp16 - 5.5128781167229297) + 0.054695746327994003*tmp80 + tmp87*(1.3565413275792199e-7*tmp186 + 1.3565413275792199e-7*tmp187 + 1.3565413275792199e-7*tmp188 + 1.3565413275792199e-7*tmp189 - 5.2691105819295299*tmp26 - 0.25536296823610299*tmp33 + 0.021878298531197701*tmp80 + 1.8242333680060501) + tmp88)*(tmp26*(tmp180*tmp194 + tmp183*tmp194 + tmp185*tmp194 + 0.00170924207274981*tmp84) + 0.28004222119932898*tmp33*tmp55 + tmp33*(0.021614102076040202*a6 + 3.5538327399018002*r - 0.17904586032964401*tmp0 + 0.059835058618348702*tmp16 - 0.089522930164821907*tmp35 + 0.041342980936911501*tmp80 - 0.0011394947151665399*tmp81 + 7.9998399449253697) + tmp60*(-0.75893966826498704*r - 0.082685961873822905*tmp16 - 5.5128781167229297) + 0.054695746327994003*tmp80 + tmp87*(1.3565413275792199e-7*tmp186 + 1.3565413275792199e-7*tmp187 + 1.3565413275792199e-7*tmp188 + 1.3565413275792199e-7*tmp189 - 5.2691105819295299*tmp26 - 0.25536296823610299*tmp33 + 0.021878298531197701*tmp80 + 1.8242333680060501) + tmp88));
const double tmp215 = r*(-4804510165749.8398*tmp26 - 14024920797448.4) - tmp210*tmp49 + tmp65;
const double tmp67 = (1.0/(tmp66));
const double tmp93 = tmp35*tmp91 + 0.00013020833333333333*tmp46;
const double tmp105 = tmp66/tmp77;
const double tmp145 = tmp140*tmp142*tmp144;
const double tmp159 = 7680.0*tmp35*tmp91;
const double tmp197 = -tmp177 + tmp196;
const double tmp216 = r*tmp214*(267544166400.0*r*tmp33 + 10215490662751.4*r - tmp210*tmp25 + 11575876386816.0*tmp25*tmp49 + 5807150888816.3496*tmp26 - tmp49*(-1759846072320.0*r - tmp209) + tmp62 + 5178202125747.6201*tmp63 - 53501685054374.0) + r*tmp215*(2.9850542633858698e-26*r*(-11614301777632.688*tmp26 - 3645361092898.4199*tmp33 + 107003370108748.0) + 2.9850542633858698e-26*tmp16*(-7767303188621.4199*tmp26 - 401316249600.0*tmp33 - 15323235994127.1) - tmp25*(2.6266180105408469e-14*tmp0 + 2.9850542633858698e-26*tmp211 + 2.9850542633858698e-26*tmp212 + 2.9850542633858698e-26*tmp213) - 1.2953135892561021e-12*tmp26 - 1.2875602280240172e-12*tmp33 - tmp49*(2.9850542633858698e-26*r*(8323596288000.0*tmp26 + 24297540157440.0) + 7.8798540316225408e-14*tmp16 + 6.0039521104596725e-13*tmp26 + 1.5049811303639448e-13*tmp33 - 3.110011922886593e-12) - 3.4554619160892902e-13*tmp49 - 1.7277309580446501e-13*tmp55 - 1.7694999344523e-13*tmp60 - 2.5856021952532902e-12*tmp72 - 5.6250597993815013e-12)/((r*tmp55 + r*(7.4971949956958071*tmp26 + 7.4523190200933342*tmp33 + 1.0241756253849399*tmp60 + 14.9653056988661*tmp72 + 32.557498453044211) + tmp0*(0.44732700598333702*tmp26 + 0.023112216946676398*tmp33 + 0.88248097348258103) - tmp16*(-1.0033194368644476*tmp26 - 0.31491016067260402*tmp33 + 9.2436517575996895) + tmp26*(-23.986035091424775*tmp26 - 0.40733420936276099*tmp33 + 0.13867330168005801*tmp60 + 9.5728999245036892) + tmp49*(0.152027027027027*tmp0 + 1.7277309580446501e-13*tmp211 + 1.7277309580446501e-13*tmp212 + 1.7277309580446501e-13*tmp213))*(r*tmp55 + r*(7.4971949956958071*tmp26 + 7.4523190200933342*tmp33 + 1.0241756253849399*tmp60 + 14.9653056988661*tmp72 + 32.557498453044211) + tmp0*(0.44732700598333702*tmp26 + 0.023112216946676398*tmp33 + 0.88248097348258103) - tmp16*(-1.0033194368644476*tmp26 - 0.31491016067260402*tmp33 + 9.2436517575996895) + tmp26*(-23.986035091424775*tmp26 - 0.40733420936276099*tmp33 + 0.13867330168005801*tmp60 + 9.5728999245036892) + tmp49*(0.152027027027027*tmp0 + 1.7277309580446501e-13*tmp211 + 1.7277309580446501e-13*tmp212 + 1.7277309580446501e-13*tmp213))) + tmp214*tmp215;
const double tmp94 = (1.0/((tmp93)*(tmp93)));
const double tmp106 = tmp105*tmp80*tmp91;
const double tmp139 = (1.0/((tmp93)*(tmp93)*(tmp93)*(tmp93)));
const double tmp160 = tmp159 + tmp46;
const double tmp217 = r*tmp105;
const double tmp79 = tmp25*tmp67*tmp77;
const double tmp95 = 1.6954210069444444e-8*tmp94;
const double tmp107 = tmp100 + tmp104 + 7680.0*tmp106 + tmp46;
const double tmp152 = tmp112 + tmp114 + tmp115 + tmp118 + tmp119 + 0.12195486878044901*tmp120*tmp123 + tmp122 + 147.44375299014601*tmp127 + 1.4827534202436501*tmp130 - 11.317508579186301*tmp133 + tmp138*tmp30*tmp7 + 1.2727731413908503e-19*tmp139*tmp145*tmp151 + tmp30 + 1;
const double tmp161 = tmp149*(tmp7*(-175.0/64.0*tmp26 - 225.0/64.0) + tmp98*(-tmp157 - 9.0/64.0) + tmp99*(117.0/32.0 - 39.0/16.0*tmp26)) + tmp160 + tmp36*((9.0/8.0)*tmp7 + tmp98*((1.0/2.0)*tmp26 + 1.0/8.0) - 5.0/4.0*tmp99);
const double tmp219 = sqrt(tmp217);
const double tmp224 = 2*tmp100 + 2*tmp104 + 15360.0*tmp106 + tmp221;
const double tmp225 = (1.0/(tmp160));
const double tmp162 = tmp154*tmp161;
const double tmp171 = tmp48*tmp79*tmp95;
const double tmp220 = tmp219/tmp47;
const double tmp222 = (1.0/(tmp219));
const double tmp163 = sqrt(tmp162*(tmp107*tmp45*tmp48*tmp79*tmp95 + tmp152));
const double tmp172 = tmp107*tmp171*tmp45 + tmp152;
const double tmp223 = tmp160*tmp177*tmp219/tmp48 + tmp160*tmp216*tmp222/(tmp221 + 2) + tmp197*tmp220;
const double tmp164 = (1.0/sqrt(tmp26*(2*tmp163 + 2*tmp43*tmp9 - 2) + 1));
const double tmp173 = sqrt(tmp162*tmp172);
const double tmp231 = tmp173/(3.3908420138888889e-8*tmp107*tmp45*tmp48*tmp79*tmp94 + 2*tmp112 + 2*tmp114 + 2*tmp115 + 2*tmp118 + 2*tmp119 + 0.24390973756089801*tmp120*tmp123 + 2*tmp122 + 294.88750598029202*tmp127 + 2.9655068404873002*tmp130 - 22.635017158372602*tmp133 + 2*tmp138*tmp30*tmp7 + 2.5455462827817007e-19*tmp139*tmp145*tmp151 + 2*tmp30 + 2);
const double Hreal_prmr = tmp164*(tmp153*tmp173*((1.0/2.0)*tmp162*(39.611280027152098*m1*m2*tmp116*tmp125*tmp2 - 663.49688845565697*pow(r, -5.5)*tmp110*tmp26 - tmp10*tmp110*(16.0*m1*m2*tmp2 - 12.0*tmp33) - tmp10*tmp116*(-5.5660152739001196*tmp26 - 10.800000000000001*tmp33 + 12.0*tmp57*tmp58*tmp59) - tmp10*tmp120*(2.7795550199225598*tmp26 + 6.6768404729664397*tmp33 + 6.8571428571428603*tmp60 - 12.0*tmp88) + 1.2727731413908503e-19*tmp110*tmp139*tmp140*tmp142*tmp144*tmp176*(-3*tmp202*tmp99 - 3*tmp203*tmp7 - 3*tmp205*tmp98) - tmp110*tmp149*tmp199 + tmp110*tmp36*(-51.695238095238103*tmp123 + 118.40000000000001*tmp25*tmp33) - tmp111*(278.13313285486299*tmp26 - 393.0*tmp33 + 30.0*tmp60) - 3.7427765505058777e-20*tmp111*tmp218*tmp223*((tmp47)*(tmp47)*(tmp47)*(tmp47)*(tmp47))*tmp67*tmp77/(pow(tmp217, 3.0/2.0)*((tmp93)*(tmp93)*(tmp93)*(tmp93)*(tmp93))) - tmp116*tmp36*(-101.934636651131*tmp26 - 268.58949820837*tmp33 + 564.0*tmp57*tmp58*tmp59 - 42.0*tmp88) - 3.7068835506091302*tmp120*tmp131*tmp26 - 0.12195486878044901*tmp120*tmp17*tmp26 - 2*tmp138*tmp34*tmp7 + tmp17*tmp2*tmp29*tmp6*(-tmp134*(-2*r + 2*tmp17*tmp2*tmp6)/((tmp135)*(tmp135)) + 2*tmp136*tmp17) - tmp222*tmp223*tmp224*tmp225*tmp45*((tmp47)*(tmp47)*(tmp47))*tmp79*tmp95 + 1.6954210069444444e-8*tmp25*tmp45*tmp48*tmp67*tmp77*tmp94*(tmp149*(-tmp175*(tmp103 - 0.578125) - 4*tmp7*(-1171.0/64.0*tmp26 - 13.453125) - 4*tmp99*((13.0/16.0)*tmp26 + 14.03125)) + tmp159*tmp216 - tmp177 + tmp196*tmp217 + tmp36*((63.0/8.0)*tmp13*tmp20*tmp21*tmp5 - 3*tmp7*(3*tmp26 + 2.8125) - 3*tmp98*((3.0/4.0)*tmp26 - 0.1875))) - 2*tmp34) + (1.0/2.0)*tmp172*(tmp154*(tmp149*(5*tmp13*tmp20*tmp21*tmp5 - tmp175*((1.0/2.0)*tmp26 + 0.125) - 9.0/2.0*tmp7) - tmp176*(5*tmp7*(-175.0/64.0*tmp26 - 3.515625) + 5*tmp98*(-tmp157 - 0.140625) + 5*tmp99*(3.65625 - 39.0/16.0*tmp26)) + tmp197) + tmp161*(tmp10*tmp7*(4*tmp25 + 2) + 2*tmp36*tmp7)/((tmp153)*(tmp153))))/(tmp161*tmp172) + tmp43*(-3*tmp16 - tmp7)/((tmp8)*(tmp8)) + tmp9*(-3*pphi*tmp15*tmp36 + pphi*(tmp39*(-tmp10*((109.0/96.0)*tmp19*tmp31*tmp32 - 177.0/16.0*tmp26 - 0.15625) - tmp149*tmp167*tmp37 - tmp165*tmp34 - tmp17*((23.0/32.0)*tmp26 - 0.09375) - tmp29*tmp36*(-801.0/128.0*tmp26 - 1591.0/256.0*tmp33 + 0.69140625)) + tmp41*(tmp10*tmp169*tmp29 - tmp10*(-1.0/16.0*tmp26 + (103.0/96.0)*tmp33 + 0.15625) - tmp149*tmp170*tmp37 - tmp17*((11.0/32.0)*tmp26 + 0.09375) - tmp29*tmp36*(-105.0/128.0*tmp26 - 613.0/256.0*tmp33 - 0.69140625))) - 1.0/2.0*tmp11*tmp23));
const double Hreal_prmprstar = tmp164*tmp231*(prstar*tmp171*tmp224 + tmp10*tmp227*(370.84417713981799*tmp26 - 524.0*tmp33 + 40.0*tmp60) + tmp10*tmp228*(-203.869273302262*tmp26 - 537.17899641674001*tmp33 + 1128.0*tmp57*tmp58*tmp59 - 84.0*tmp88) + 0.97563895024359204*tmp123*tmp229 + 589.77501196058404*tmp125*tmp227*tmp26 + 11.862027361949201*tmp128*tmp229*tmp26 - 67.905051475117801*tmp131*tmp228*tmp26 + 9.7748977258817303e-16*tmp145*tmp149*tmp218*tmp225*tmp227/((tmp93)*(tmp93)*(tmp93)) + tmp17*tmp227*(32.0*m1*m2*tmp2 - 24.0*tmp33) + tmp17*tmp228*(-16.698045821700401*tmp26 - 32.399999999999999*tmp33 + 36.0*tmp57*tmp58*tmp59) + tmp17*tmp229*(11.1182200796902*tmp26 + 26.707361891865801*tmp33 + 27.428571428571399*tmp60 - 48.0*tmp88) + tmp199*tmp227*tmp36);
const double Hreal_prmpphi = tmp164*(tmp231*(pphi*tmp138*tmp221 + 2*pphi*tmp17) + tmp9*(pphi*(tmp39*(tmp11*(-267.0/64.0*tmp26 - 1591.0/384.0*tmp33 + 0.4609375) + tmp165*tmp18 + tmp167*tmp232) + tmp41*(tmp11*(-35.0/64.0*tmp26 - 613.0/384.0*tmp33 - 0.4609375) - tmp169*tmp18 + tmp170*tmp232)) + tmp10*tmp15 + (1.0/4.0)*tmp17*tmp23 + tmp42));
const double Hreal = sqrt(2*tmp26*(tmp163 + tmp43*tmp9 - 1) + 1);
const double xi = tmp160*tmp220;
dvalues[0] = Hreal_prmr;
dvalues[1] = 0.;
dvalues[2] = Hreal_prmprstar;
dvalues[3] = Hreal_prmpphi;
hamiltonian_and_xi[0] = Hreal;
hamiltonian_and_xi[1] = xi;
return GSL_SUCCESS;
}

import numpy as np
def v5HM_optimized_hamiltonian( m1 , m2 , r, prstar , pphi , chi1 , chi2 , a6 , dSO):
    tmp0 = ((r)*(r)*(r))
    tmp1 = m1 + m2
    tmp12 = ((r)*(r))
    tmp17 = (1.0/(r))
    tmp21 = ((m1)*(m1))
    tmp22 = ((m2)*(m2))
    tmp25 = ((r)*(r)*(r)*(r))
    tmp30 = ((r)*(r)*(r)*(r)*(r))
    tmp35 = np.log(r)
    tmp37 = ((m1)*(m1)*(m1))
    tmp38 = ((m2)*(m2)*(m2))
    tmp46 = ((prstar)*(prstar)*(prstar)*(prstar))
    tmp47 = np.power(prstar, 6)
    tmp48 = np.power(prstar, 8)
    tmp2 = (1.0/((tmp1)*(tmp1)))
    tmp5 = chi1*m1 + chi2*m2
    tmp8 = (1.0/(tmp0))
    tmp10 = (1.0/((tmp1)*(tmp1)*(tmp1)))
    tmp13 = (1.0/(tmp12))
    tmp14 = (1.0/((tmp1)*(tmp1)*(tmp1)*(tmp1)))
    tmp15 = chi1*m1 - chi2*m2
    tmp26 = (1.0/(tmp25))
    tmp28 = 2*tmp17 + 1
    tmp31 = (1.0/(tmp30))
    tmp34 = ((m1)*(m1)*(m1)*(m1))*((m2)*(m2)*(m2)*(m2))/np.power(tmp1, 8)
    tmp36 = ((tmp35)*(tmp35))
    tmp39 = np.power(tmp1, -6)
    tmp41 = 8*r + 2*tmp0 + 4*tmp12 + 16.0
    tmp16 = tmp15*(m1 - m2)
    tmp18 = m1*m2*tmp2
    tmp20 = ((pphi)*(pphi))*tmp13
    tmp23 = tmp14*tmp21*tmp22
    tmp24 = ((pphi)*(pphi))*tmp8
    tmp27 = ((pphi)*(pphi)*(pphi)*(pphi))*tmp26
    tmp33 = ((tmp15)*(tmp15))*tmp2
    tmp40 = tmp37*tmp38*tmp39
    tmp7 = tmp2*((tmp5)*(tmp5))
    tmp32 = tmp10*tmp16*tmp5
    tmp53 = 14700*tmp18 + 42911.0
    tmp58 = tmp12*(-7876452485916241.0/12500000.0*tmp18 - 98886748396767.0/500000.0*tmp23 + 580530436787913.0/100000.0)
    tmp59 = tmp0*(-53511513369581.0/20000000.0*tmp18 - 138240*tmp23 - 52783413229329.0/10000000.0)
    tmp60 = tmp18*(-31383302728516087.0/12500000.0*tmp18 - 426364516032331.0/10000000.0*tmp23 + 14515200*tmp40 + 100201376401019.0/100000.0)
    tmp61 = r*(43393301259014.688*tmp18 + (3369809522332764779.0/78125.0)*tmp23 + (296393260946151.0/50.0)*tmp40 + (866182644304933.0/10.0)*((1 - 0.49694878161693501*tmp18)*(1 - 0.49694878161693501*tmp18)) + 188440788778196.0)
    tmp29 = tmp13*tmp7
    tmp43 = 756*tmp18 + 1079.0
    tmp54 = r*tmp53
    tmp62 = tmp35*(49152*r*((16368307925443.0/40000.0)*tmp18 + 102574080*tmp23 - 105983591868019.0/50000.0) + 879923036160*tmp0 + 283115520*tmp12*tmp53 - 1698693120*tmp18*(11592*tmp18 + 69847.0))
    tmp44 = 7680*(2048.0*m1*m2*tmp2*tmp35*(336*r + 756*tmp18 + 407.0) + 28*m1*m2*tmp2*(1920*a6 + 733955307463037.0/1000000000.0) - 7.0*r*((938918400156317.0/1000000000.0)*m1*m2*tmp2 - 185763092693281.0/1000000000.0*tmp23 - 245760.0) - 270820329770593.0/50000000.0*tmp23 - 3440640.0)/(32768.0*tmp18*tmp35*(240*r*(-373313530533103.0/50000000000.0*tmp18 - 3024*tmp23 + 17264.0) + 960*tmp0*tmp43 + 1920*tmp12*(588*tmp18 + 1079.0) - 388422414769507.0/10000000.0*tmp18 - 47061405915993.0/25000000.0*tmp23 + 480*tmp25*tmp43 + 161280*tmp30 + 13447680.0) + 53760*tmp18*(7680*a6*(tmp25 + tmp41) + (113485217444961.0/1000000000.0)*r*(-tmp25 + tmp41) + (7402203300817.0/50000000000.0)*r*(7704*r + 1926*tmp0 + 3852*tmp12 + 349*tmp25 + 36400.0) + 128.0*r*((33046962773603.0/2500000000.0)*r + (42646962773603.0/10000000000.0)*tmp0 + (852939255472061.0/100000000000.0)*tmp12 - 137046962773603.0/20000000000.0*tmp25 - 42153037226397.0/1250000000.0)) + 67645734912*tmp23*tmp36 + 7.0*tmp23*(-39321600*a6*(3*r + 59.0) + (186464462028901.0/250000.0)*a6 + (122635399361987.0/1000.0)*r - 308925070376879.0/50000.0*tmp0 + (206478381132587.0/100000.0)*tmp12 - 308925070376879.0/100000.0*tmp25 + (3566651379711.0/2500.0)*tmp30 + 276057889687011.0/1000.0) + 13212057600*tmp30 + (241555486248807.0/1000.0)*tmp34 + 1120*tmp40*(-163683964822551.0/1000000.0*r - 3566651379711.0/200000.0*tmp12 - 59449372951581.0/50000.0))
    tmp50 = tmp29 + 1
    tmp57 = tmp35*(-879923036160*tmp12 - 253929149957443.0/10.0*tmp18 - 5041721180160*tmp23 - 283115520*tmp54 + 104186110149937.0)
    tmp65 = 5787938193408*r*tmp36 - 9216.0*tmp58 - 967680.0*tmp59 + 55296*tmp60 + tmp61 + tmp62
    tmp45 = tmp25*tmp44 + tmp29
    tmp64 = (53058305272831.0/5.0)*r*tmp18 + (182268054644921.0/100.0)*r*tmp23 - 197383821284628.0/5.0*r + (258910106287381.0/100.0)*tmp12*tmp18 + 133772083200*tmp12*tmp23 + (510774533137571.0/100.0)*tmp12 - 60249543508726.0/5.0*tmp18 + (400296247701391.0/5.0)*tmp23 + 5787938193408*tmp36 + (336524885906151.0/50.0)*tmp40 - 163418713120743.0/500000.0*tmp54 - tmp57 + 275059053208689.0
    tmp66 = tmp64/tmp65
    Hreal = np.sqrt(2*tmp18*(np.sqrt((tmp26*(-5.0/4.0*tmp32 + tmp33*((1.0/2.0)*tmp18 + 1.0/8.0) + (9.0/8.0)*tmp7) + tmp31*(tmp32*(117.0/32.0 - 39.0/16.0*tmp18) + tmp33*((21.0/16.0)*tmp14*tmp21*tmp22 - 81.0/64.0*tmp18 - 9.0/64.0) + tmp7*(-175.0/64.0*tmp18 - 225.0/64.0)) + tmp45)*(((prstar)*(prstar))*tmp17*((tmp50)*(tmp50))*tmp65*(tmp26*(tmp32*((13.0/16.0)*tmp18 + 449.0/32.0) + tmp33*((115.0/64.0)*tmp18 + (1.0/16.0)*tmp23 - 37.0/64.0) + tmp7*(-1171.0/64.0*tmp18 - 861.0/64.0)) + tmp29 + tmp30*tmp44*tmp66 + tmp8*(-21.0/8.0*tmp32 + tmp33*((3.0/4.0)*tmp18 - 3.0/16.0) + tmp7*(3*tmp18 + 45.0/16.0)))/(((tmp45)*(tmp45))*tmp64) + tmp13*tmp46*(8*m1*m2*tmp2 - 6*tmp23) + tmp13*tmp47*(-139150381847503.0/50000000000000.0*tmp18 - 27.0/5.0*tmp23 + 6*tmp37*tmp38*tmp39) + tmp13*tmp48*((4343054718629.0/3125000000000.0)*tmp18 + (166921011824161.0/50000000000000.0)*tmp23 - 6*tmp34 + (342857142857143.0/100000000000000.0)*tmp40) + (121954868780449.0/1000000000000000.0)*tmp17*tmp18*tmp48 - tmp20*tmp28*tmp7/(tmp12 + tmp28*tmp7) + tmp20 + tmp26*tmp46*(tmp18*(452542166996693.0/1000000000000.0 - 516952380952381.0/10000000000000.0*tmp35) + tmp23*((592.0/5.0)*tmp35 - 179613660498019.0/100000000000.0) + (150579635104141.0/250000000000.0)*tmp40) + 0.0004427880404175718*tmp31*tmp46*((tmp50)*(tmp50)*(tmp50)*(tmp50))*(tmp32*((45.0/8.0)*tmp18 - 5.0/16.0) + tmp33*((75.0/32.0)*m1*m2*tmp2 - 15.0/8.0*tmp23 - 15.0/32.0) + tmp7*((165.0/32.0)*tmp18 - 5*tmp23 + 25.0/32.0))*((r*tmp36 - 1.5922768509339455e-9*tmp58 - 1.6718906934806428e-7*tmp59 + (3.0/314015744.0)*tmp60 + (1.0/5787938193408.0)*tmp61 + (1.0/5787938193408.0)*tmp62)*(r*tmp36 - 1.5922768509339455e-9*tmp58 - 1.6718906934806428e-7*tmp59 + (3.0/314015744.0)*tmp60 + (1.0/5787938193408.0)*tmp61 + (1.0/5787938193408.0)*tmp62))/(((tmp45)*(tmp45)*(tmp45)*(tmp45))*((0.038579573843421422*r*tmp18 + 0.0066265062908739489*r*tmp23 - 0.14352105046684044*r + 0.0094128916415248589*tmp12*tmp18 + 0.00048633950287942822*tmp12*tmp23 + 0.018569631763766858*tmp12 - 0.04380844244600398*tmp18 + 0.29106204142837921*tmp23 + 0.021042529325572337*tmp36 + 0.024469282648975563*tmp40 - 1.18824456940711e-6*tmp54 - 3.6355829351353647e-15*tmp57 + 1)*(0.038579573843421422*r*tmp18 + 0.0066265062908739489*r*tmp23 - 0.14352105046684044*r + 0.0094128916415248589*tmp12*tmp18 + 0.00048633950287942822*tmp12*tmp23 + 0.018569631763766858*tmp12 - 0.04380844244600398*tmp18 + 0.29106204142837921*tmp23 + 0.021042529325572337*tmp36 + 0.024469282648975563*tmp40 - 1.18824456940711e-6*tmp54 - 3.6355829351353647e-15*tmp57 + 1))) + tmp46*tmp8*((115888805356193.0/1250000000000.0)*tmp18 - 131*tmp23 + 10*tmp40) + tmp47*tmp8*(-84945530542609.0/2500000000000.0*tmp18 - 447649163680617.0/5000000000000.0*tmp23 - 14*tmp34 + 188*tmp37*tmp38*tmp39) + 1 + (29655068404873.0/20000000000000.0)*tmp18*tmp48/np.power(r, 5.0/2.0) - 113175085791863.0/10000000000000.0*tmp18*tmp47/np.power(r, 7.0/2.0) + (73721876495073.0/500000000000.0)*tmp18*tmp46/np.power(r, 9.0/2.0))/(tmp28*tmp29 + 1)) - 1 + (dSO*m1*m2*pphi*tmp10*tmp5*tmp8 + (1.0/4.0)*pphi*tmp13*(-tmp10*((tmp5)*(tmp5)*(tmp5)) + tmp14*tmp16*((tmp5)*(tmp5))) + pphi*(tmp16*tmp2*(tmp13*(-1.0/32.0*tmp18 + (103.0/192.0)*tmp23 + 5.0/64.0) + tmp17*((11.0/32.0)*tmp18 + 3.0/32.0) + tmp20*(15.0/32.0 - 9.0/32.0*tmp18) + tmp24*(-35.0/128.0*tmp18 - 613.0/768.0*tmp23 - 59.0/256.0) + tmp27*((75.0/256.0)*tmp14*tmp21*tmp22 - 45.0/128.0*tmp18 - 105.0/256.0) + 1.0/4.0) + tmp5*(tmp13*((109.0/192.0)*tmp14*tmp21*tmp22 - 177.0/32.0*tmp18 - 5.0/64.0) + tmp17*((23.0/32.0)*tmp18 - 3.0/32.0) + tmp20*(-45.0/32.0*tmp18 - 15.0/32.0) + tmp24*(-267.0/128.0*tmp18 - 1591.0/768.0*tmp23 + 59.0/256.0) + tmp27*((75.0/128.0)*tmp18 + (345.0/256.0)*tmp23 + 105.0/256.0) + 7.0/4.0)/tmp1))/(tmp0 + tmp7*(r + 2))) + 1)
    xi = tmp45*np.sqrt(r*tmp66)/tmp50
    return Hreal , xi

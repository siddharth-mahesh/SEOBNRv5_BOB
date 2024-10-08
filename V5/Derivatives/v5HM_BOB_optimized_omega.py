import numpy as np
def v5HM_BOB_optimized_omega(m1, m2, r, pphi, chi1, chi2, a6, dSO):
    tmp0 = ((r)*(r)*(r))
    tmp1 = m1 + m2
    tmp15 = ((r)*(r))
    tmp22 = (1.0/(r))
    tmp27 = ((m1)*(m1))
    tmp28 = ((m2)*(m2))
    tmp31 = ((r)*(r)*(r)*(r))
    tmp44 = ((r)*(r)*(r)*(r)*(r))
    tmp45 = np.log(r)
    tmp2 = (1.0/((tmp1)*(tmp1)))
    tmp5 = chi1*m1 + chi2*m2
    tmp9 = (1.0/(tmp0))
    tmp12 = (1.0/((tmp1)*(tmp1)*(tmp1)))
    tmp16 = (1.0/(tmp15))
    tmp18 = (1.0/((tmp1)*(tmp1)*(tmp1)*(tmp1)))
    tmp19 = chi1*m1 - chi2*m2
    tmp32 = (1.0/(tmp31))
    tmp37 = 2*tmp22 + 1
    tmp47 = 8.0*r + 2.0*tmp0 + 4.0*tmp15 + 16.0
    tmp10 = pphi*tmp9
    tmp17 = pphi*tmp16
    tmp20 = tmp19*(m1 - m2)
    tmp23 = m1*m2*tmp2
    tmp26 = ((pphi)*(pphi))*tmp16
    tmp29 = tmp18*tmp27*tmp28
    tmp30 = ((pphi)*(pphi))*tmp9
    tmp33 = ((pphi)*(pphi)*(pphi)*(pphi))*tmp32
    tmp34 = tmp5/tmp1
    tmp42 = ((tmp19)*(tmp19))*tmp2
    tmp50 = ((pphi)*(pphi)*(pphi))*tmp32
    tmp7 = tmp2*((tmp5)*(tmp5))
    tmp14 = dSO*m1*m2*tmp12*tmp5
    tmp21 = -1.0/4.0*tmp12*((tmp5)*(tmp5)*(tmp5)) + (1.0/4.0)*tmp18*tmp20*((tmp5)*(tmp5))
    tmp43 = tmp12*tmp20*tmp5
    tmp46 = 756.0*tmp23
    tmp8 = (1.0/(tmp0 + tmp7*(r + 2)))
    tmp36 = tmp2*tmp20*(tmp16*(-1.0/32.0*tmp23 + (103.0/192.0)*tmp29 + 5.0/64.0) + tmp22*((11.0/32.0)*tmp23 + 3.0/32.0) + tmp26*(15.0/32.0 - 9.0/32.0*tmp23) + tmp30*(-35.0/128.0*tmp23 - 613.0/768.0*tmp29 - 59.0/256.0) + tmp33*((75.0/256.0)*tmp18*tmp27*tmp28 - 45.0/128.0*tmp23 - 105.0/256.0) + 1.0/4.0) + tmp34*(tmp16*((109.0/192.0)*tmp18*tmp27*tmp28 - 177.0/32.0*tmp23 - 5.0/64.0) + tmp22*((23.0/32.0)*tmp23 - 3.0/32.0) + tmp26*(-45.0/32.0*tmp23 - 15.0/32.0) + tmp30*(-267.0/128.0*tmp23 - 1591.0/768.0*tmp29 + 59.0/256.0) + tmp33*((75.0/128.0)*tmp23 + (345.0/256.0)*tmp29 + 105.0/256.0) + 7.0/4.0)
    tmp48 = tmp46 + 1079.0
    tmp39 = -tmp37/(tmp15 + tmp37*tmp7)
    tmp49 = np.sqrt((tmp26*tmp39*tmp7 + tmp26 + 1)*(tmp16*tmp7 + 7680.0*tmp31*(2048.0*m1*m2*tmp2*tmp45*(336.0*r + tmp46 + 407.0) + 28.0*m1*m2*tmp2*(1920.0*a6 + 733955.30746303697) - 7.0*r*(938918.40015631705*m1*m2*tmp2 - 185763.09269328101*tmp29 - 245760.0) - 5416406.5954118604*tmp29 - 3440640.0)/(241555486248.80701*((m1)*(m1)*(m1)*(m1))*((m2)*(m2)*(m2)*(m2))/np.power(tmp1, 8) + 1120.0*((m1)*(m1)*(m1))*((m2)*(m2)*(m2))*(-163683964.82255101*r - 17833256.898554999*tmp15 - 1188987459.03162)/np.power(tmp1, 6) + 32768.0*tmp23*tmp45*(240.0*r*(-7466.2706106620599*tmp23 - 3024.0*tmp29 + 17264.0) + 960.0*tmp0*tmp48 + 1920.0*tmp15*(588.0*tmp23 + 1079.0) - 38842241.476950698*tmp23 - 1882456.2366397199*tmp29 + 480.0*tmp31*tmp48 + 161280.0*tmp44 + 13447680.0) + 53760.0*tmp23*(7680.0*a6*(tmp31 + tmp47) + 113485.217444961*r*(-tmp31 + tmp47) + 148.04406601634*r*(7704.0*r + 1926.0*tmp0 + 3852.0*tmp15 + 349.0*tmp31 + 36400.0) + 128.0*r*(13218.7851094412*r + 4264.6962773603*tmp0 + 8529.3925547206109*tmp15 - 6852.3481386801504*tmp31 - 33722.4297811176)) + 67645734912.0*tmp29*((tmp45)*(tmp45)) + 7.0*tmp29*(-39321600.0*a6*(3.0*r + 59.0) + 745857848.11560404*a6 + 122635399361.987*r - 6178501407.5375795*tmp0 + 2064783811.32587*tmp15 - 3089250703.7687898*tmp31 + 1426660551.8843999*tmp44 + 276057889687.01099) + 13212057600.0*tmp44) + tmp32*(tmp42*((1.0/2.0)*tmp23 + 1.0/8.0) - 5.0/4.0*tmp43 + (9.0/8.0)*tmp7) + (tmp42*((21.0/16.0)*tmp18*tmp27*tmp28 - 81.0/64.0*tmp23 - 9.0/64.0) + tmp43*(117.0/32.0 - 39.0/16.0*tmp23) + tmp7*(-175.0/64.0*tmp23 - 225.0/64.0))/tmp44)/(tmp16*tmp37*tmp7 + 1))
    omega = (tmp49*(2*pphi*tmp16*tmp39*tmp7 + 2*tmp17)/(2*tmp26*tmp39*tmp7 + 2*tmp26 + 2) + tmp8*(pphi*(tmp2*tmp20*(tmp10*(-35.0/64.0*tmp23 - 613.0/384.0*tmp29 - 0.4609375) + tmp17*(0.9375 - 9.0/16.0*tmp23) + tmp50*((75.0/64.0)*tmp18*tmp27*tmp28 - 45.0/32.0*tmp23 - 1.640625)) + tmp34*(tmp10*(-267.0/64.0*tmp23 - 1591.0/384.0*tmp29 + 0.4609375) + tmp17*(-45.0/16.0*tmp23 - 0.9375) + tmp50*((75.0/32.0)*tmp23 + (345.0/64.0)*tmp29 + 1.640625))) + tmp14*tmp9 + tmp16*tmp21 + tmp36))/np.sqrt(tmp23*(2*tmp49 + 2*tmp8*(pphi*tmp36 + tmp10*tmp14 + tmp17*tmp21) - 2) + 1)
    return omega

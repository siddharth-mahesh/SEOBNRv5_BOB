import numpy as np
def v5HM_BOB_optimized_nqc_terms(t,t0,hNR,omegaNR,omegaQNM,tau):
    damp = np.zeros(3)
    domega = np.zeros(2)
    tmp0 = (1.0/(tau))
    tmp1 = 2*np.log(omegaNR/omegaQNM)
    tmp5 = (1.0/16.0)*((omegaNR)*(omegaNR)*(omegaNR)*(omegaNR))
    tmp24 = (1.0/((tau)*(tau)))
    tmp2 = tmp0*(t - t0 + tau*tmp1)
    tmp6 = np.tanh(tmp1)
    tmp9 = (1.0/16.0)*((omegaQNM)*(omegaQNM)*(omegaQNM)*(omegaQNM)) - tmp5
    tmp12 = hNR*((omegaNR)*(omegaNR))*np.cosh(tmp1)
    tmp3 = np.cosh(tmp2)
    tmp7 = np.tanh(tmp2)
    tmp8 = 1 - tmp6
    tmp16 = np.sinh(tmp2)
    tmp4 = (1.0/(tmp3))
    tmp10 = tmp9/tmp8
    tmp17 = tmp16/((tmp3)*(tmp3))
    tmp18 = 1 - ((tmp7)*(tmp7))
    tmp11 = tmp10*(-tmp6 + tmp7) + tmp5
    tmp13 = (1.0/np.sqrt(tmp11))*tmp12
    tmp20 = tmp0*tmp10*tmp18
    tmp21 = np.power(tmp11, -1.5)
    tmp25 = 0.25*tmp10*tmp18*tmp24
    tmp23 = tmp12*tmp21*tmp4
    tmp15 = (1.0/4.0)*tmp13*tmp4
    damp[0] = tmp15
    damp[1] = -1.0/4.0*tmp0*tmp13*tmp17 - 0.125*tmp20*tmp23
    damp[2] = 0.1875*np.power(tmp11, -2.5)*tmp12*((tmp18)*(tmp18))*tmp24*tmp4*((tmp9)*(tmp9))/((tmp8)*(tmp8)) + tmp12*tmp17*tmp21*tmp25 + (1.0/2.0)*tmp13*((tmp16)*(tmp16))*tmp24/((tmp3)*(tmp3)*(tmp3)) - tmp15*tmp24 + tmp23*tmp25*tmp7
    domega[0] = 2*np.power(tmp11, 0.25)
    domega[1] = 0.5*np.power(tmp11, -0.75)*tmp20
    return damp, domega
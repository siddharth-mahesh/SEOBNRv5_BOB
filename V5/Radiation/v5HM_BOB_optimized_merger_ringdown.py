import numpy as np
def v5HM_BOB_optimized_merger_ringdown(t,t0,hNR,omegaNR,omegaQNM,tau):
    tmp0 = ((omegaNR)*(omegaNR))
    tmp1 = (1.0/(omegaQNM))
    tmp5 = (1.0/16.0)*((omegaNR)*(omegaNR)*(omegaNR)*(omegaNR))
    tmp16 = omegaQNM/omegaNR
    tmp2 = omegaNR*tmp1
    tmp3 = 2*np.log(tmp2)
    tmp12 = (1.0/2.0)*tmp0*tmp1
    tmp13 = tau*tmp0*tmp1
    tmp4 = (t - t0 + tau*tmp3)/tau
    tmp6 = np.tanh(tmp3)
    tmp7 = tmp5 + ((1.0/16.0)*((omegaQNM)*(omegaQNM)*(omegaQNM)*(omegaQNM)) - tmp5)*(-tmp6 + np.tanh(tmp4))/(1 - tmp6)
    tmp10 = np.power(tmp7, 0.25)
    tmp15 = 2*tmp1*tmp10
    tmp17 = 2*omegaQNM*tmp10/tmp0
    h = (1.0/4.0)*hNR*tmp0*(1.0/np.sqrt(tmp7))*np.cosh(tmp3)/np.cosh(tmp4)
    phi = omegaQNM*tau*(-np.arctan2((1.0/2.0)*omegaNR, (1.0/2.0)*omegaQNM) + np.arctan2(tmp10, (1.0/2.0)*omegaQNM)) + 0.5*omegaQNM*tau*np.log((1 - tmp2)*(tmp15 + 1)/((1 - tmp15)*(tmp2 + 1))) - tmp13*(-np.arctan2((1.0/2.0)*omegaNR, tmp12) + np.arctan2(tmp10, tmp12)) - 0.5*tmp13*np.log((1 - tmp16)*(tmp17 + 1)/((1 - tmp17)*(tmp16 + 1)))
    return h,phi
import numpy as np
def v5HM_BOB_optimized_nqc_terms(t,t0,hNR,omegaNR,omegaQNM,tau):
    damp = np.zeros(3)
    domega = np.zeros(2)

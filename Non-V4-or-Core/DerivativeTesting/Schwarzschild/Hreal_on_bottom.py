import numpy as np
def compute_Hreal(x=2.129681018601393e+01, y=0.000000000000000e+00, z=0.000000000000000e+00, p1=0.000000000000000e+00, p2=2.335391115580442e-01, p3=-4.235164736271502e-22):
    r = np.sqrt(x*x + y*y + z*z)
    u = 1/r
    rhat3 = z*u
    rhat2 = y*u
    rhat1 = x*u
    rhoinv = 1/np.sqrt(x*x + y*y)
    phihat3 = 0
    phihat2 = x*rhoinv
    phihat1 = -y*rhoinv
    thetahat3 = phihat1*rhat2 - phihat2*rhat1
    thetahat2 = -phihat1*rhat3 + phihat3*rhat1
    thetahat1 = phihat2*rhat3 - phihat3*rhat2
    costheta = z*u
    sin2theta = 1 - costheta*costheta
    A = 1 - 2*u
    B = 1/A
    pphi = (phihat1*p1 + phihat2*p2 + phihat3*p3)*r
    pphi2 = pphi*pphi*sin2theta
    pr = rhat1*p1 + rhat2*p2 + rhat3*p3
    ptheta = (thetahat1*p1 + thetahat2*p2 + thetahat3*p3)*r
    gammappsum = pr*pr/B + ptheta*ptheta/r/r + pphi2/r/r/sin2theta
    Hnsradicand = 1 + gammappsum
    alpha = np.sqrt(A)
    Hns = alpha*np.sqrt(Hnsradicand)
    Heff = Hns
    return Heff
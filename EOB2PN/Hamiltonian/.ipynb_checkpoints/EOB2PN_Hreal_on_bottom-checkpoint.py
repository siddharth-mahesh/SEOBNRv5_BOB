import numpy as np
def compute_EOB2PN_Hreal(m1=23., m2=10., q1=2.1242072413581923e+01, q2=0., q3=0., p1=0., p2=2.1696072000958128e-01, p3=1.0000000000000000e-03):
    M = m1 + m2
    mu=m1*m2/M
    eta=mu/M
    q=np.sqrt(q1*q1+q2*q2+q3*q3)
    u=1/q
    psq=p1*p1+p2*p2+p3*p3
    n3=q3*u
    n2 = q2*u
    n1 = q1*u
    ndotp=n1*p1+n2*p2+n3*p3
    qprime3=q3
    qprime2 = q2
    qprime1 = q1
    pprime3=p3
    pprime2 = p2 
    pprime1 = p1 
    qprimesq=qprime1*qprime1+qprime2*qprime2+qprime3*qprime3
    qprime = np.sqrt(qprimesq)
    uprime=1/qprime
    nprime3=qprime3*uprime
    nprime2 = qprime2*uprime
    nprime1 = qprime1*uprime
    pprimedotnprime=pprime1*nprime1+pprime2*nprime2+pprime3*nprime3
    nprimedotpprimesq = pprimedotnprime*pprimedotnprime
    pprimesq=pprime1*pprime1+pprime2*pprime2+pprime3*pprime3
    B=1+2*uprime+(4-6*eta)*uprime*uprime
    Binv = 1/B
    A=1-2*uprime+2*eta*uprime*uprime*uprime
    Heff=np.sqrt(A*(1+pprimesq+nprimedotpprimesq*(Binv-1)))
    Hreal=(np.sqrt(1+2*eta*(Heff-1))-1)/eta
    return Hreal
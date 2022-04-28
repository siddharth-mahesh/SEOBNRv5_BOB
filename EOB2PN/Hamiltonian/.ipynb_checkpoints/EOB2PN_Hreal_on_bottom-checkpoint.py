import numpy as np
def compute_EOB2PN_Hreal(m1=23., m2=10., q1=2.1242072413581923e+01, q2=0., q3=0., p1=0., p2=2.1696072000958128e-01, p3=1.0000000000000000e-03):
    M = m1 + m2
    mu=m1*m2/M
    eta=mu/M
    q=np.sqrt(q1*q1+q2*q2+q3*q3)
    u=1/q
    psq=p1*p1+p2*p2+p3*p3
    n3=q3/q
    n2 = q2/q
    n1 = q1/q
    ndotp=n1*p1+n2*p2+n3*p3
    qprime2PN3=(1+0.125*eta)*eta*ndotp*ndotp*n3+0.25*(5-0.5*eta)*eta*psq*n3+1.5*(1-0.5*eta)*eta*ndotp*p3+0.25*(1-7*eta+eta*eta)*u*n3+0.125*(1-eta)*eta*q*psq*psq*n3+0.5*(1+eta)*eta*q*ndotp*psq*p3
    qprime2PN2 = (1 + 0.125*eta)*eta*ndotp*ndotp*n2 + 0.25*(5 - 0.5*eta)*eta*psq*n2 + 1.5*(1 - 0.5*eta)*eta*ndotp*p2  + 0.25*(1 - 7*eta + eta*eta)*u*n2 + 0.125*(1 - eta)*eta*q*psq*psq*n2 + 0.5*(1 + eta)*eta*q*ndotp*psq*p2
    qprime2PN1 = (1 + 0.125*eta)*eta*ndotp*ndotp*n1 + 0.25*(5 - 0.5*eta)*eta*psq*n1 + 1.5*(1 - 0.5*eta)*eta*ndotp*p1  + 0.25*(1 - 7*eta + eta*eta)*u*n1 + 0.125*(1 - eta)*eta*q*psq*psq*n1 + 0.5*(1 + eta)*eta*q*ndotp*psq*p1
    qprime1PN3=(1+0.5*eta)*n3+0.5*eta*psq*q3+eta*q*ndotp*p3
    qprime1PN2 = (1 + 0.5*eta)*n2 + 0.5*eta*psq*q2 + eta*q*ndotp*p2
    qprime1PN1 = (1 + 0.5*eta)*n1 + 0.5*eta*psq*q1 + eta*q*ndotp*p1
    qprime3=q3+qprime1PN3+qprime2PN3
    qprime2 = q2 + qprime1PN2 + qprime2PN2
    qprime1 = q1 + qprime1PN1 + qprime2PN1
    pprime2PN3=0.125*eta*psq*psq*(-1+3*eta)*p3+0.25*(3+11*eta)*u*u*p3-0.75*(3+0.5*eta)*ndotp*u*u*n3+0.25*(-2-18*eta+eta*eta)*u*u*n3+0.125*(10-eta)*eta*u*ndotp*psq*n3-0.125*(16+5*eta)*eta*u*ndotp*ndotp*p3+0.375*eta*u*ndotp*ndotp*ndotp*n3
    pprime2PN2 = 0.125*eta*psq*psq*(-1 + 3*eta)*p2 + 0.25*(3 + 11*eta)*u*u*p2 - 0.75*(3 + 0.5*eta)*ndotp*u*u*n2 + 0.25*(-2 -18*eta+ eta*eta)*u*u*n2 + 0.125*(10 - eta)*eta*u*ndotp*psq*n2 - 0.125*(16+5*eta)*eta*u*ndotp*ndotp*p2 + 0.375*eta*u*ndotp*ndotp*ndotp*n2
    pprime2PN1 = 0.125*eta*psq*psq*(-1 + 3*eta)*p1 + 0.25*(3 + 11*eta)*u*u*p1 - 0.75*(3 + 0.5*eta)*ndotp*u*u*n1 + 0.25*(-2 -18*eta+ eta*eta)*u*u*n1 + 0.125*(10 - eta)*eta*u*ndotp*psq*n1 - 0.125*(16+5*eta)*eta*u*ndotp*ndotp*p1 + 0.375*eta*u*ndotp*ndotp*ndotp*n1
    pprime1PN3=-(1+0.5*eta)*u*p3+0.5*eta*psq*p3+(1+0.5*eta)*u*ndotp*n3
    pprime1PN2 = -(1 + 0.5*eta)*u*p2 + 0.5*eta*psq*p2 + (1 + 0.5*eta)*u*ndotp*n2
    pprime1PN1 = -(1 + 0.5*eta)*u*p1 + 0.5*eta*psq*p1 + (1 + 0.5*eta)*u*ndotp*n1
    pprime3=p3+pprime1PN3+pprime2PN3
    pprime2 = p2 + pprime1PN2 + pprime2PN2
    pprime1 = p1 + pprime1PN1 + pprime2PN1
    qprimesq=qprime1*qprime1+qprime2*qprime2+qprime3*qprime3
    qprime = np.sqrt(qprimesq)
    uprime=1/qprime
    nprime3=qprime3/qprime
    nprime2 = qprime2/qprime
    nprime1 = qprime1/qprime
    pprimedotnprime=pprime1*nprime1+pprime2*nprime2+pprime3*nprime3
    nprimedotpprimesq = pprimedotnprime*pprimedotnprime
    pprimesq=pprime1*pprime1+pprime2*pprime2+pprime3*pprime3
    B=1+2*uprime+(4-6*eta)*uprime*uprime
    Binv = 1/B
    A=1-2*uprime+2*eta*uprime*uprime*uprime
    Heff=np.sqrt(A*(1+pprimesq+nprimedotpprimesq*(Binv-1)))
    Hreal=np.sqrt(1+2*eta*(Heff-1))/eta
    return Hreal
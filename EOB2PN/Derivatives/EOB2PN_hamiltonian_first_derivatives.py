from __future__ import division
import numpy as np
def ham_first_derivs(m1, m2, q1, q2, q3, p1, p2, p3):
    M = m1+m2
    mu = m1*m2/M
    eta = mu/M
    q = np.sqrt(q1*q1+q2*q2+q3*q3)
    u = 1/q
    psq = p1*p1+p2*p2+p3*p3
    n3 = q3/q
    n2 = q2/q
    n1 = q1/q
    ndotp = n1*p1+n2*p2+n3*p3
    qprime2PN3 = (1+0.125*eta)*eta*ndotp*ndotp*n3+0.25*(5-0.5*eta)*eta*psq*n3+1.5*(1-0.5*eta)*eta*ndotp*p3+0.25*(1-7*eta+eta*eta)*u*n3+0.125*(1-eta)*eta*q*psq*psq*n3+0.5*(1+eta)*eta*q*ndotp*psq*p3
    qprime2PN2 = (1+0.125*eta)*eta*ndotp*ndotp*n2+0.25*(5-0.5*eta)*eta*psq*n2+1.5*(1-0.5*eta)*eta*ndotp*p2+0.25*(1-7*eta+eta*eta)*u*n2+0.125*(1-eta)*eta*q*psq*psq*n2+0.5*(1+eta)*eta*q*ndotp*psq*p2
    qprime2PN1 = (1+0.125*eta)*eta*ndotp*ndotp*n1+0.25*(5-0.5*eta)*eta*psq*n1+1.5*(1-0.5*eta)*eta*ndotp*p1+0.25*(1-7*eta+eta*eta)*u*n1+0.125*(1-eta)*eta*q*psq*psq*n1+0.5*(1+eta)*eta*q*ndotp*psq*p1
    qprime1PN3 = (1+0.5*eta)*n3+0.5*eta*psq*q3+eta*q*ndotp*p3
    qprime1PN2 = (1+0.5*eta)*n2+0.5*eta*psq*q2+eta*q*ndotp*p2
    qprime1PN1 = (1+0.5*eta)*n1+0.5*eta*psq*q1+eta*q*ndotp*p1
    qprime3 = q3+qprime1PN3+qprime2PN3
    qprime2 = q2+qprime1PN2+qprime2PN2
    qprime1 = q1+qprime1PN1+qprime2PN1
    pprime2PN3 = 0.125*eta*psq*psq*(-1+3*eta)*p3+0.25*(3+11*eta)*u*u*p3-0.75*(3+0.5*eta)*ndotp*u*u*n3+0.25*(-2-18*eta+eta*eta)*u*u*n3+0.125*(10-eta)*eta*u*ndotp*psq*n3-0.125*(16+5*eta)*eta*u*ndotp*ndotp*p3+0.375*eta*u*ndotp*ndotp*ndotp*n3
    pprime2PN2 = 0.125*eta*psq*psq*(-1+3*eta)*p2+0.25*(3+11*eta)*u*u*p2-0.75*(3+0.5*eta)*ndotp*u*u*n2+0.25*(-2-18*eta+eta*eta)*u*u*n2+0.125*(10-eta)*eta*u*ndotp*psq*n2-0.125*(16+5*eta)*eta*u*ndotp*ndotp*p2+0.375*eta*u*ndotp*ndotp*ndotp*n2
    pprime2PN1 = 0.125*eta*psq*psq*(-1+3*eta)*p1+0.25*(3+11*eta)*u*u*p1-0.75*(3+0.5*eta)*ndotp*u*u*n1+0.25*(-2-18*eta+eta*eta)*u*u*n1+0.125*(10-eta)*eta*u*ndotp*psq*n1-0.125*(16+5*eta)*eta*u*ndotp*ndotp*p1+0.375*eta*u*ndotp*ndotp*ndotp*n1
    pprime1PN3 = -(1+0.5*eta)*u*p3+0.5*eta*psq*p3+(1+0.5*eta)*u*ndotp*n3
    pprime1PN2 = -(1+0.5*eta)*u*p2+0.5*eta*psq*p2+(1+0.5*eta)*u*ndotp*n2
    pprime1PN1 = -(1+0.5*eta)*u*p1+0.5*eta*psq*p1+(1+0.5*eta)*u*ndotp*n1
    pprime3 = p3+pprime1PN3+pprime2PN3
    pprime2 = p2+pprime1PN2+pprime2PN2
    pprime1 = p1+pprime1PN1+pprime2PN1
    qprimesq = qprime1*qprime1+qprime2*qprime2+qprime3*qprime3
    qprime = np.sqrt(qprimesq)
    uprime = 1/qprime
    nprime3 = qprime3/qprime
    nprime2 = qprime2/qprime
    nprime1 = qprime1/qprime
    pprimedotnprime = pprime1*nprime1+pprime2*nprime2+pprime3*nprime3
    nprimedotpprimesq = pprimedotnprime*pprimedotnprime
    pprimesq = pprime1*pprime1+pprime2*pprime2+pprime3*pprime3
    B = 1+2*uprime+(4-6*eta)*uprime*uprime
    Binv = 1/B
    A = 1-2*uprime+2*eta*uprime*uprime*uprime
    Heff = np.sqrt(A*(1+pprimesq+nprimedotpprimesq*(Binv-1)))
    Hreal = np.sqrt(1+2*eta*(Heff-1))/eta
    q_prmq1 = q1/np.sqrt(q1**2 + q2**2 + q3**2)
    q_prmq2 = q2/np.sqrt(q1**2 + q2**2 + q3**2)
    q_prmq3 = q3/np.sqrt(q1**2 + q2**2 + q3**2)
    u_prmq1 = -q_prmq1/q**2
    u_prmq2 = -q_prmq2/q**2
    u_prmq3 = -q_prmq3/q**2
    psq_prmp1 = 2*p1
    psq_prmp2 = 2*p2
    psq_prmp3 = 2*p3
    n3_prmq1 = -q3*q_prmq1/q**2
    n3_prmq2 = -q3*q_prmq2/q**2
    n3_prmq3 = 1/q - q3*q_prmq3/q**2
    n2_prmq1 = -q2*q_prmq1/q**2
    n2_prmq2 = 1/q - q2*q_prmq2/q**2
    n2_prmq3 = -q2*q_prmq3/q**2
    n1_prmq1 = 1/q - q1*q_prmq1/q**2
    n1_prmq2 = -q1*q_prmq2/q**2
    n1_prmq3 = -q1*q_prmq3/q**2
    ndotp_prmq1 = n1_prmq1*p1 + n2_prmq1*p2 + n3_prmq1*p3
    ndotp_prmq2 = n1_prmq2*p1 + n2_prmq2*p2 + n3_prmq2*p3
    ndotp_prmq3 = n1_prmq3*p1 + n2_prmq3*p2 + n3_prmq3*p3
    ndotp_prmp1 = n1
    ndotp_prmp2 = n2
    ndotp_prmp3 = n3
    qprime2PN3_prmq1 = 2*eta*n3*ndotp*ndotp_prmq1*(0.125*eta + 1) + eta*n3*psq**2*q_prmq1*(0.125 - 0.125*eta) + eta*n3_prmq1*ndotp**2*(0.125*eta + 1) + eta*n3_prmq1*psq**2*q*(0.125 - 0.125*eta) + eta*n3_prmq1*psq*(1.25 - 0.125*eta) + eta*ndotp*p3*psq*q_prmq1*(0.5*eta + 0.5) + eta*ndotp_prmq1*p3*psq*q*(0.5*eta + 0.5) + eta*ndotp_prmq1*p3*(1.5 - 0.75*eta) + n3*u_prmq1*(0.25*eta**2 - 1.75*eta + 0.25) + n3_prmq1*u*(0.25*eta**2 - 1.75*eta + 0.25)
    qprime2PN3_prmq2 = 2*eta*n3*ndotp*ndotp_prmq2*(0.125*eta + 1) + eta*n3*psq**2*q_prmq2*(0.125 - 0.125*eta) + eta*n3_prmq2*ndotp**2*(0.125*eta + 1) + eta*n3_prmq2*psq**2*q*(0.125 - 0.125*eta) + eta*n3_prmq2*psq*(1.25 - 0.125*eta) + eta*ndotp*p3*psq*q_prmq2*(0.5*eta + 0.5) + eta*ndotp_prmq2*p3*psq*q*(0.5*eta + 0.5) + eta*ndotp_prmq2*p3*(1.5 - 0.75*eta) + n3*u_prmq2*(0.25*eta**2 - 1.75*eta + 0.25) + n3_prmq2*u*(0.25*eta**2 - 1.75*eta + 0.25)
    qprime2PN3_prmq3 = 2*eta*n3*ndotp*ndotp_prmq3*(0.125*eta + 1) + eta*n3*psq**2*q_prmq3*(0.125 - 0.125*eta) + eta*n3_prmq3*ndotp**2*(0.125*eta + 1) + eta*n3_prmq3*psq**2*q*(0.125 - 0.125*eta) + eta*n3_prmq3*psq*(1.25 - 0.125*eta) + eta*ndotp*p3*psq*q_prmq3*(0.5*eta + 0.5) + eta*ndotp_prmq3*p3*psq*q*(0.5*eta + 0.5) + eta*ndotp_prmq3*p3*(1.5 - 0.75*eta) + n3*u_prmq3*(0.25*eta**2 - 1.75*eta + 0.25) + n3_prmq3*u*(0.25*eta**2 - 1.75*eta + 0.25)
    qprime2PN3_prmp1 = 2*eta*n3*ndotp*ndotp_prmp1*(0.125*eta + 1) + 2*eta*n3*psq*psq_prmp1*q*(0.125 - 0.125*eta) + eta*n3*psq_prmp1*(1.25 - 0.125*eta) + eta*ndotp*p3*psq_prmp1*q*(0.5*eta + 0.5) + eta*ndotp_prmp1*p3*psq*q*(0.5*eta + 0.5) + eta*ndotp_prmp1*p3*(1.5 - 0.75*eta)
    qprime2PN3_prmp2 = 2*eta*n3*ndotp*ndotp_prmp2*(0.125*eta + 1) + 2*eta*n3*psq*psq_prmp2*q*(0.125 - 0.125*eta) + eta*n3*psq_prmp2*(1.25 - 0.125*eta) + eta*ndotp*p3*psq_prmp2*q*(0.5*eta + 0.5) + eta*ndotp_prmp2*p3*psq*q*(0.5*eta + 0.5) + eta*ndotp_prmp2*p3*(1.5 - 0.75*eta)
    qprime2PN3_prmp3 = 2*eta*n3*ndotp*ndotp_prmp3*(0.125*eta + 1) + 2*eta*n3*psq*psq_prmp3*q*(0.125 - 0.125*eta) + eta*n3*psq_prmp3*(1.25 - 0.125*eta) + eta*ndotp*p3*psq_prmp3*q*(0.5*eta + 0.5) + eta*ndotp*psq*q*(0.5*eta + 0.5) + eta*ndotp*(1.5 - 0.75*eta) + eta*ndotp_prmp3*p3*psq*q*(0.5*eta + 0.5) + eta*ndotp_prmp3*p3*(1.5 - 0.75*eta)
    qprime2PN2_prmq1 = 2*eta*n2*ndotp*ndotp_prmq1*(0.125*eta + 1) + eta*n2*psq**2*q_prmq1*(0.125 - 0.125*eta) + eta*n2_prmq1*ndotp**2*(0.125*eta + 1) + eta*n2_prmq1*psq**2*q*(0.125 - 0.125*eta) + eta*n2_prmq1*psq*(1.25 - 0.125*eta) + eta*ndotp*p2*psq*q_prmq1*(0.5*eta + 0.5) + eta*ndotp_prmq1*p2*psq*q*(0.5*eta + 0.5) + eta*ndotp_prmq1*p2*(1.5 - 0.75*eta) + n2*u_prmq1*(0.25*eta**2 - 1.75*eta + 0.25) + n2_prmq1*u*(0.25*eta**2 - 1.75*eta + 0.25)
    qprime2PN2_prmq2 = 2*eta*n2*ndotp*ndotp_prmq2*(0.125*eta + 1) + eta*n2*psq**2*q_prmq2*(0.125 - 0.125*eta) + eta*n2_prmq2*ndotp**2*(0.125*eta + 1) + eta*n2_prmq2*psq**2*q*(0.125 - 0.125*eta) + eta*n2_prmq2*psq*(1.25 - 0.125*eta) + eta*ndotp*p2*psq*q_prmq2*(0.5*eta + 0.5) + eta*ndotp_prmq2*p2*psq*q*(0.5*eta + 0.5) + eta*ndotp_prmq2*p2*(1.5 - 0.75*eta) + n2*u_prmq2*(0.25*eta**2 - 1.75*eta + 0.25) + n2_prmq2*u*(0.25*eta**2 - 1.75*eta + 0.25)
    qprime2PN2_prmq3 = 2*eta*n2*ndotp*ndotp_prmq3*(0.125*eta + 1) + eta*n2*psq**2*q_prmq3*(0.125 - 0.125*eta) + eta*n2_prmq3*ndotp**2*(0.125*eta + 1) + eta*n2_prmq3*psq**2*q*(0.125 - 0.125*eta) + eta*n2_prmq3*psq*(1.25 - 0.125*eta) + eta*ndotp*p2*psq*q_prmq3*(0.5*eta + 0.5) + eta*ndotp_prmq3*p2*psq*q*(0.5*eta + 0.5) + eta*ndotp_prmq3*p2*(1.5 - 0.75*eta) + n2*u_prmq3*(0.25*eta**2 - 1.75*eta + 0.25) + n2_prmq3*u*(0.25*eta**2 - 1.75*eta + 0.25)
    qprime2PN2_prmp1 = 2*eta*n2*ndotp*ndotp_prmp1*(0.125*eta + 1) + 2*eta*n2*psq*psq_prmp1*q*(0.125 - 0.125*eta) + eta*n2*psq_prmp1*(1.25 - 0.125*eta) + eta*ndotp*p2*psq_prmp1*q*(0.5*eta + 0.5) + eta*ndotp_prmp1*p2*psq*q*(0.5*eta + 0.5) + eta*ndotp_prmp1*p2*(1.5 - 0.75*eta)
    qprime2PN2_prmp2 = 2*eta*n2*ndotp*ndotp_prmp2*(0.125*eta + 1) + 2*eta*n2*psq*psq_prmp2*q*(0.125 - 0.125*eta) + eta*n2*psq_prmp2*(1.25 - 0.125*eta) + eta*ndotp*p2*psq_prmp2*q*(0.5*eta + 0.5) + eta*ndotp*psq*q*(0.5*eta + 0.5) + eta*ndotp*(1.5 - 0.75*eta) + eta*ndotp_prmp2*p2*psq*q*(0.5*eta + 0.5) + eta*ndotp_prmp2*p2*(1.5 - 0.75*eta)
    qprime2PN2_prmp3 = 2*eta*n2*ndotp*ndotp_prmp3*(0.125*eta + 1) + 2*eta*n2*psq*psq_prmp3*q*(0.125 - 0.125*eta) + eta*n2*psq_prmp3*(1.25 - 0.125*eta) + eta*ndotp*p2*psq_prmp3*q*(0.5*eta + 0.5) + eta*ndotp_prmp3*p2*psq*q*(0.5*eta + 0.5) + eta*ndotp_prmp3*p2*(1.5 - 0.75*eta)
    qprime2PN1_prmq1 = 2*eta*n1*ndotp*ndotp_prmq1*(0.125*eta + 1) + eta*n1*psq**2*q_prmq1*(0.125 - 0.125*eta) + eta*n1_prmq1*ndotp**2*(0.125*eta + 1) + eta*n1_prmq1*psq**2*q*(0.125 - 0.125*eta) + eta*n1_prmq1*psq*(1.25 - 0.125*eta) + eta*ndotp*p1*psq*q_prmq1*(0.5*eta + 0.5) + eta*ndotp_prmq1*p1*psq*q*(0.5*eta + 0.5) + eta*ndotp_prmq1*p1*(1.5 - 0.75*eta) + n1*u_prmq1*(0.25*eta**2 - 1.75*eta + 0.25) + n1_prmq1*u*(0.25*eta**2 - 1.75*eta + 0.25)
    qprime2PN1_prmq2 = 2*eta*n1*ndotp*ndotp_prmq2*(0.125*eta + 1) + eta*n1*psq**2*q_prmq2*(0.125 - 0.125*eta) + eta*n1_prmq2*ndotp**2*(0.125*eta + 1) + eta*n1_prmq2*psq**2*q*(0.125 - 0.125*eta) + eta*n1_prmq2*psq*(1.25 - 0.125*eta) + eta*ndotp*p1*psq*q_prmq2*(0.5*eta + 0.5) + eta*ndotp_prmq2*p1*psq*q*(0.5*eta + 0.5) + eta*ndotp_prmq2*p1*(1.5 - 0.75*eta) + n1*u_prmq2*(0.25*eta**2 - 1.75*eta + 0.25) + n1_prmq2*u*(0.25*eta**2 - 1.75*eta + 0.25)
    qprime2PN1_prmq3 = 2*eta*n1*ndotp*ndotp_prmq3*(0.125*eta + 1) + eta*n1*psq**2*q_prmq3*(0.125 - 0.125*eta) + eta*n1_prmq3*ndotp**2*(0.125*eta + 1) + eta*n1_prmq3*psq**2*q*(0.125 - 0.125*eta) + eta*n1_prmq3*psq*(1.25 - 0.125*eta) + eta*ndotp*p1*psq*q_prmq3*(0.5*eta + 0.5) + eta*ndotp_prmq3*p1*psq*q*(0.5*eta + 0.5) + eta*ndotp_prmq3*p1*(1.5 - 0.75*eta) + n1*u_prmq3*(0.25*eta**2 - 1.75*eta + 0.25) + n1_prmq3*u*(0.25*eta**2 - 1.75*eta + 0.25)
    qprime2PN1_prmp1 = 2*eta*n1*ndotp*ndotp_prmp1*(0.125*eta + 1) + 2*eta*n1*psq*psq_prmp1*q*(0.125 - 0.125*eta) + eta*n1*psq_prmp1*(1.25 - 0.125*eta) + eta*ndotp*p1*psq_prmp1*q*(0.5*eta + 0.5) + eta*ndotp*psq*q*(0.5*eta + 0.5) + eta*ndotp*(1.5 - 0.75*eta) + eta*ndotp_prmp1*p1*psq*q*(0.5*eta + 0.5) + eta*ndotp_prmp1*p1*(1.5 - 0.75*eta)
    qprime2PN1_prmp2 = 2*eta*n1*ndotp*ndotp_prmp2*(0.125*eta + 1) + 2*eta*n1*psq*psq_prmp2*q*(0.125 - 0.125*eta) + eta*n1*psq_prmp2*(1.25 - 0.125*eta) + eta*ndotp*p1*psq_prmp2*q*(0.5*eta + 0.5) + eta*ndotp_prmp2*p1*psq*q*(0.5*eta + 0.5) + eta*ndotp_prmp2*p1*(1.5 - 0.75*eta)
    qprime2PN1_prmp3 = 2*eta*n1*ndotp*ndotp_prmp3*(0.125*eta + 1) + 2*eta*n1*psq*psq_prmp3*q*(0.125 - 0.125*eta) + eta*n1*psq_prmp3*(1.25 - 0.125*eta) + eta*ndotp*p1*psq_prmp3*q*(0.5*eta + 0.5) + eta*ndotp_prmp3*p1*psq*q*(0.5*eta + 0.5) + eta*ndotp_prmp3*p1*(1.5 - 0.75*eta)
    qprime1PN3_prmq1 = eta*ndotp*p3*q_prmq1 + eta*ndotp_prmq1*p3*q + n3_prmq1*(0.5*eta + 1)
    qprime1PN3_prmq2 = eta*ndotp*p3*q_prmq2 + eta*ndotp_prmq2*p3*q + n3_prmq2*(0.5*eta + 1)
    qprime1PN3_prmq3 = eta*ndotp*p3*q_prmq3 + eta*ndotp_prmq3*p3*q + 0.5*eta*psq + n3_prmq3*(0.5*eta + 1)
    qprime1PN3_prmp1 = eta*ndotp_prmp1*p3*q + 0.5*eta*psq_prmp1*q3
    qprime1PN3_prmp2 = eta*ndotp_prmp2*p3*q + 0.5*eta*psq_prmp2*q3
    qprime1PN3_prmp3 = eta*ndotp*q + eta*ndotp_prmp3*p3*q + 0.5*eta*psq_prmp3*q3
    qprime1PN2_prmq1 = eta*ndotp*p2*q_prmq1 + eta*ndotp_prmq1*p2*q + n2_prmq1*(0.5*eta + 1)
    qprime1PN2_prmq2 = eta*ndotp*p2*q_prmq2 + eta*ndotp_prmq2*p2*q + 0.5*eta*psq + n2_prmq2*(0.5*eta + 1)
    qprime1PN2_prmq3 = eta*ndotp*p2*q_prmq3 + eta*ndotp_prmq3*p2*q + n2_prmq3*(0.5*eta + 1)
    qprime1PN2_prmp1 = eta*ndotp_prmp1*p2*q + 0.5*eta*psq_prmp1*q2
    qprime1PN2_prmp2 = eta*ndotp*q + eta*ndotp_prmp2*p2*q + 0.5*eta*psq_prmp2*q2
    qprime1PN2_prmp3 = eta*ndotp_prmp3*p2*q + 0.5*eta*psq_prmp3*q2
    qprime1PN1_prmq1 = eta*ndotp*p1*q_prmq1 + eta*ndotp_prmq1*p1*q + 0.5*eta*psq + n1_prmq1*(0.5*eta + 1)
    qprime1PN1_prmq2 = eta*ndotp*p1*q_prmq2 + eta*ndotp_prmq2*p1*q + n1_prmq2*(0.5*eta + 1)
    qprime1PN1_prmq3 = eta*ndotp*p1*q_prmq3 + eta*ndotp_prmq3*p1*q + n1_prmq3*(0.5*eta + 1)
    qprime1PN1_prmp1 = eta*ndotp*q + eta*ndotp_prmp1*p1*q + 0.5*eta*psq_prmp1*q1
    qprime1PN1_prmp2 = eta*ndotp_prmp2*p1*q + 0.5*eta*psq_prmp2*q1
    qprime1PN1_prmp3 = eta*ndotp_prmp3*p1*q + 0.5*eta*psq_prmp3*q1
    qprime3_prmq1 = qprime1PN3_prmq1 + qprime2PN3_prmq1
    qprime3_prmq2 = qprime1PN3_prmq2 + qprime2PN3_prmq2
    qprime3_prmq3 = qprime1PN3_prmq3 + qprime2PN3_prmq3 + 1
    qprime3_prmp1 = qprime1PN3_prmp1 + qprime2PN3_prmp1
    qprime3_prmp2 = qprime1PN3_prmp2 + qprime2PN3_prmp2
    qprime3_prmp3 = qprime1PN3_prmp3 + qprime2PN3_prmp3
    qprime2_prmq1 = qprime1PN2_prmq1 + qprime2PN2_prmq1
    qprime2_prmq2 = qprime1PN2_prmq2 + qprime2PN2_prmq2 + 1
    qprime2_prmq3 = qprime1PN2_prmq3 + qprime2PN2_prmq3
    qprime2_prmp1 = qprime1PN2_prmp1 + qprime2PN2_prmp1
    qprime2_prmp2 = qprime1PN2_prmp2 + qprime2PN2_prmp2
    qprime2_prmp3 = qprime1PN2_prmp3 + qprime2PN2_prmp3
    qprime1_prmq1 = qprime1PN1_prmq1 + qprime2PN1_prmq1 + 1
    qprime1_prmq2 = qprime1PN1_prmq2 + qprime2PN1_prmq2
    qprime1_prmq3 = qprime1PN1_prmq3 + qprime2PN1_prmq3
    qprime1_prmp1 = qprime1PN1_prmp1 + qprime2PN1_prmp1
    qprime1_prmp2 = qprime1PN1_prmp2 + qprime2PN1_prmp2
    qprime1_prmp3 = qprime1PN1_prmp3 + qprime2PN1_prmp3
    pprime2PN3_prmq1 = 0.375*eta*n3*ndotp**3*u_prmq1 + 1.125*eta*n3*ndotp**2*ndotp_prmq1*u + eta*n3*ndotp*psq*u_prmq1*(1.25 - 0.125*eta) + eta*n3*ndotp_prmq1*psq*u*(1.25 - 0.125*eta) + 0.375*eta*n3_prmq1*ndotp**3*u + eta*n3_prmq1*ndotp*psq*u*(1.25 - 0.125*eta) - eta*ndotp**2*p3*u_prmq1*(0.625*eta + 2.0) - 2*eta*ndotp*ndotp_prmq1*p3*u*(0.625*eta + 2.0) - n3*ndotp*u*u_prmq1*(0.75*eta + 4.5) - n3*ndotp_prmq1*u**2*(0.375*eta + 2.25) + n3*u*u_prmq1*(0.5*eta**2 - 9.0*eta - 1.0) - n3_prmq1*ndotp*u**2*(0.375*eta + 2.25) + n3_prmq1*u**2*(0.25*eta**2 - 4.5*eta - 0.5) + 2*p3*u*u_prmq1*(2.75*eta + 0.75)
    pprime2PN3_prmq2 = 0.375*eta*n3*ndotp**3*u_prmq2 + 1.125*eta*n3*ndotp**2*ndotp_prmq2*u + eta*n3*ndotp*psq*u_prmq2*(1.25 - 0.125*eta) + eta*n3*ndotp_prmq2*psq*u*(1.25 - 0.125*eta) + 0.375*eta*n3_prmq2*ndotp**3*u + eta*n3_prmq2*ndotp*psq*u*(1.25 - 0.125*eta) - eta*ndotp**2*p3*u_prmq2*(0.625*eta + 2.0) - 2*eta*ndotp*ndotp_prmq2*p3*u*(0.625*eta + 2.0) - n3*ndotp*u*u_prmq2*(0.75*eta + 4.5) - n3*ndotp_prmq2*u**2*(0.375*eta + 2.25) + n3*u*u_prmq2*(0.5*eta**2 - 9.0*eta - 1.0) - n3_prmq2*ndotp*u**2*(0.375*eta + 2.25) + n3_prmq2*u**2*(0.25*eta**2 - 4.5*eta - 0.5) + 2*p3*u*u_prmq2*(2.75*eta + 0.75)
    pprime2PN3_prmq3 = 0.375*eta*n3*ndotp**3*u_prmq3 + 1.125*eta*n3*ndotp**2*ndotp_prmq3*u + eta*n3*ndotp*psq*u_prmq3*(1.25 - 0.125*eta) + eta*n3*ndotp_prmq3*psq*u*(1.25 - 0.125*eta) + 0.375*eta*n3_prmq3*ndotp**3*u + eta*n3_prmq3*ndotp*psq*u*(1.25 - 0.125*eta) - eta*ndotp**2*p3*u_prmq3*(0.625*eta + 2.0) - 2*eta*ndotp*ndotp_prmq3*p3*u*(0.625*eta + 2.0) - n3*ndotp*u*u_prmq3*(0.75*eta + 4.5) - n3*ndotp_prmq3*u**2*(0.375*eta + 2.25) + n3*u*u_prmq3*(0.5*eta**2 - 9.0*eta - 1.0) - n3_prmq3*ndotp*u**2*(0.375*eta + 2.25) + n3_prmq3*u**2*(0.25*eta**2 - 4.5*eta - 0.5) + 2*p3*u*u_prmq3*(2.75*eta + 0.75)
    pprime2PN3_prmp1 = 1.125*eta*n3*ndotp**2*ndotp_prmp1*u + eta*n3*ndotp*psq_prmp1*u*(1.25 - 0.125*eta) + eta*n3*ndotp_prmp1*psq*u*(1.25 - 0.125*eta) - 2*eta*ndotp*ndotp_prmp1*p3*u*(0.625*eta + 2.0) + 0.25*eta*p3*psq*psq_prmp1*(3*eta - 1) - n3*ndotp_prmp1*u**2*(0.375*eta + 2.25)
    pprime2PN3_prmp2 = 1.125*eta*n3*ndotp**2*ndotp_prmp2*u + eta*n3*ndotp*psq_prmp2*u*(1.25 - 0.125*eta) + eta*n3*ndotp_prmp2*psq*u*(1.25 - 0.125*eta) - 2*eta*ndotp*ndotp_prmp2*p3*u*(0.625*eta + 2.0) + 0.25*eta*p3*psq*psq_prmp2*(3*eta - 1) - n3*ndotp_prmp2*u**2*(0.375*eta + 2.25)
    pprime2PN3_prmp3 = 1.125*eta*n3*ndotp**2*ndotp_prmp3*u + eta*n3*ndotp*psq_prmp3*u*(1.25 - 0.125*eta) + eta*n3*ndotp_prmp3*psq*u*(1.25 - 0.125*eta) - eta*ndotp**2*u*(0.625*eta + 2.0) - 2*eta*ndotp*ndotp_prmp3*p3*u*(0.625*eta + 2.0) + 0.25*eta*p3*psq*psq_prmp3*(3*eta - 1) + 0.125*eta*psq**2*(3*eta - 1) - n3*ndotp_prmp3*u**2*(0.375*eta + 2.25) + u**2*(2.75*eta + 0.75)
    pprime2PN2_prmq1 = 0.375*eta*n2*ndotp**3*u_prmq1 + 1.125*eta*n2*ndotp**2*ndotp_prmq1*u + eta*n2*ndotp*psq*u_prmq1*(1.25 - 0.125*eta) + eta*n2*ndotp_prmq1*psq*u*(1.25 - 0.125*eta) + 0.375*eta*n2_prmq1*ndotp**3*u + eta*n2_prmq1*ndotp*psq*u*(1.25 - 0.125*eta) - eta*ndotp**2*p2*u_prmq1*(0.625*eta + 2.0) - 2*eta*ndotp*ndotp_prmq1*p2*u*(0.625*eta + 2.0) - n2*ndotp*u*u_prmq1*(0.75*eta + 4.5) - n2*ndotp_prmq1*u**2*(0.375*eta + 2.25) + n2*u*u_prmq1*(0.5*eta**2 - 9.0*eta - 1.0) - n2_prmq1*ndotp*u**2*(0.375*eta + 2.25) + n2_prmq1*u**2*(0.25*eta**2 - 4.5*eta - 0.5) + 2*p2*u*u_prmq1*(2.75*eta + 0.75)
    pprime2PN2_prmq2 = 0.375*eta*n2*ndotp**3*u_prmq2 + 1.125*eta*n2*ndotp**2*ndotp_prmq2*u + eta*n2*ndotp*psq*u_prmq2*(1.25 - 0.125*eta) + eta*n2*ndotp_prmq2*psq*u*(1.25 - 0.125*eta) + 0.375*eta*n2_prmq2*ndotp**3*u + eta*n2_prmq2*ndotp*psq*u*(1.25 - 0.125*eta) - eta*ndotp**2*p2*u_prmq2*(0.625*eta + 2.0) - 2*eta*ndotp*ndotp_prmq2*p2*u*(0.625*eta + 2.0) - n2*ndotp*u*u_prmq2*(0.75*eta + 4.5) - n2*ndotp_prmq2*u**2*(0.375*eta + 2.25) + n2*u*u_prmq2*(0.5*eta**2 - 9.0*eta - 1.0) - n2_prmq2*ndotp*u**2*(0.375*eta + 2.25) + n2_prmq2*u**2*(0.25*eta**2 - 4.5*eta - 0.5) + 2*p2*u*u_prmq2*(2.75*eta + 0.75)
    pprime2PN2_prmq3 = 0.375*eta*n2*ndotp**3*u_prmq3 + 1.125*eta*n2*ndotp**2*ndotp_prmq3*u + eta*n2*ndotp*psq*u_prmq3*(1.25 - 0.125*eta) + eta*n2*ndotp_prmq3*psq*u*(1.25 - 0.125*eta) + 0.375*eta*n2_prmq3*ndotp**3*u + eta*n2_prmq3*ndotp*psq*u*(1.25 - 0.125*eta) - eta*ndotp**2*p2*u_prmq3*(0.625*eta + 2.0) - 2*eta*ndotp*ndotp_prmq3*p2*u*(0.625*eta + 2.0) - n2*ndotp*u*u_prmq3*(0.75*eta + 4.5) - n2*ndotp_prmq3*u**2*(0.375*eta + 2.25) + n2*u*u_prmq3*(0.5*eta**2 - 9.0*eta - 1.0) - n2_prmq3*ndotp*u**2*(0.375*eta + 2.25) + n2_prmq3*u**2*(0.25*eta**2 - 4.5*eta - 0.5) + 2*p2*u*u_prmq3*(2.75*eta + 0.75)
    pprime2PN2_prmp1 = 1.125*eta*n2*ndotp**2*ndotp_prmp1*u + eta*n2*ndotp*psq_prmp1*u*(1.25 - 0.125*eta) + eta*n2*ndotp_prmp1*psq*u*(1.25 - 0.125*eta) - 2*eta*ndotp*ndotp_prmp1*p2*u*(0.625*eta + 2.0) + 0.25*eta*p2*psq*psq_prmp1*(3*eta - 1) - n2*ndotp_prmp1*u**2*(0.375*eta + 2.25)
    pprime2PN2_prmp2 = 1.125*eta*n2*ndotp**2*ndotp_prmp2*u + eta*n2*ndotp*psq_prmp2*u*(1.25 - 0.125*eta) + eta*n2*ndotp_prmp2*psq*u*(1.25 - 0.125*eta) - eta*ndotp**2*u*(0.625*eta + 2.0) - 2*eta*ndotp*ndotp_prmp2*p2*u*(0.625*eta + 2.0) + 0.25*eta*p2*psq*psq_prmp2*(3*eta - 1) + 0.125*eta*psq**2*(3*eta - 1) - n2*ndotp_prmp2*u**2*(0.375*eta + 2.25) + u**2*(2.75*eta + 0.75)
    pprime2PN2_prmp3 = 1.125*eta*n2*ndotp**2*ndotp_prmp3*u + eta*n2*ndotp*psq_prmp3*u*(1.25 - 0.125*eta) + eta*n2*ndotp_prmp3*psq*u*(1.25 - 0.125*eta) - 2*eta*ndotp*ndotp_prmp3*p2*u*(0.625*eta + 2.0) + 0.25*eta*p2*psq*psq_prmp3*(3*eta - 1) - n2*ndotp_prmp3*u**2*(0.375*eta + 2.25)
    pprime2PN1_prmq1 = 0.375*eta*n1*ndotp**3*u_prmq1 + 1.125*eta*n1*ndotp**2*ndotp_prmq1*u + eta*n1*ndotp*psq*u_prmq1*(1.25 - 0.125*eta) + eta*n1*ndotp_prmq1*psq*u*(1.25 - 0.125*eta) + 0.375*eta*n1_prmq1*ndotp**3*u + eta*n1_prmq1*ndotp*psq*u*(1.25 - 0.125*eta) - eta*ndotp**2*p1*u_prmq1*(0.625*eta + 2.0) - 2*eta*ndotp*ndotp_prmq1*p1*u*(0.625*eta + 2.0) - n1*ndotp*u*u_prmq1*(0.75*eta + 4.5) - n1*ndotp_prmq1*u**2*(0.375*eta + 2.25) + n1*u*u_prmq1*(0.5*eta**2 - 9.0*eta - 1.0) - n1_prmq1*ndotp*u**2*(0.375*eta + 2.25) + n1_prmq1*u**2*(0.25*eta**2 - 4.5*eta - 0.5) + 2*p1*u*u_prmq1*(2.75*eta + 0.75)
    pprime2PN1_prmq2 = 0.375*eta*n1*ndotp**3*u_prmq2 + 1.125*eta*n1*ndotp**2*ndotp_prmq2*u + eta*n1*ndotp*psq*u_prmq2*(1.25 - 0.125*eta) + eta*n1*ndotp_prmq2*psq*u*(1.25 - 0.125*eta) + 0.375*eta*n1_prmq2*ndotp**3*u + eta*n1_prmq2*ndotp*psq*u*(1.25 - 0.125*eta) - eta*ndotp**2*p1*u_prmq2*(0.625*eta + 2.0) - 2*eta*ndotp*ndotp_prmq2*p1*u*(0.625*eta + 2.0) - n1*ndotp*u*u_prmq2*(0.75*eta + 4.5) - n1*ndotp_prmq2*u**2*(0.375*eta + 2.25) + n1*u*u_prmq2*(0.5*eta**2 - 9.0*eta - 1.0) - n1_prmq2*ndotp*u**2*(0.375*eta + 2.25) + n1_prmq2*u**2*(0.25*eta**2 - 4.5*eta - 0.5) + 2*p1*u*u_prmq2*(2.75*eta + 0.75)
    pprime2PN1_prmq3 = 0.375*eta*n1*ndotp**3*u_prmq3 + 1.125*eta*n1*ndotp**2*ndotp_prmq3*u + eta*n1*ndotp*psq*u_prmq3*(1.25 - 0.125*eta) + eta*n1*ndotp_prmq3*psq*u*(1.25 - 0.125*eta) + 0.375*eta*n1_prmq3*ndotp**3*u + eta*n1_prmq3*ndotp*psq*u*(1.25 - 0.125*eta) - eta*ndotp**2*p1*u_prmq3*(0.625*eta + 2.0) - 2*eta*ndotp*ndotp_prmq3*p1*u*(0.625*eta + 2.0) - n1*ndotp*u*u_prmq3*(0.75*eta + 4.5) - n1*ndotp_prmq3*u**2*(0.375*eta + 2.25) + n1*u*u_prmq3*(0.5*eta**2 - 9.0*eta - 1.0) - n1_prmq3*ndotp*u**2*(0.375*eta + 2.25) + n1_prmq3*u**2*(0.25*eta**2 - 4.5*eta - 0.5) + 2*p1*u*u_prmq3*(2.75*eta + 0.75)
    pprime2PN1_prmp1 = 1.125*eta*n1*ndotp**2*ndotp_prmp1*u + eta*n1*ndotp*psq_prmp1*u*(1.25 - 0.125*eta) + eta*n1*ndotp_prmp1*psq*u*(1.25 - 0.125*eta) - eta*ndotp**2*u*(0.625*eta + 2.0) - 2*eta*ndotp*ndotp_prmp1*p1*u*(0.625*eta + 2.0) + 0.25*eta*p1*psq*psq_prmp1*(3*eta - 1) + 0.125*eta*psq**2*(3*eta - 1) - n1*ndotp_prmp1*u**2*(0.375*eta + 2.25) + u**2*(2.75*eta + 0.75)
    pprime2PN1_prmp2 = 1.125*eta*n1*ndotp**2*ndotp_prmp2*u + eta*n1*ndotp*psq_prmp2*u*(1.25 - 0.125*eta) + eta*n1*ndotp_prmp2*psq*u*(1.25 - 0.125*eta) - 2*eta*ndotp*ndotp_prmp2*p1*u*(0.625*eta + 2.0) + 0.25*eta*p1*psq*psq_prmp2*(3*eta - 1) - n1*ndotp_prmp2*u**2*(0.375*eta + 2.25)
    pprime2PN1_prmp3 = 1.125*eta*n1*ndotp**2*ndotp_prmp3*u + eta*n1*ndotp*psq_prmp3*u*(1.25 - 0.125*eta) + eta*n1*ndotp_prmp3*psq*u*(1.25 - 0.125*eta) - 2*eta*ndotp*ndotp_prmp3*p1*u*(0.625*eta + 2.0) + 0.25*eta*p1*psq*psq_prmp3*(3*eta - 1) - n1*ndotp_prmp3*u**2*(0.375*eta + 2.25)
    pprime1PN3_prmq1 = n3*ndotp*u_prmq1*(0.5*eta + 1) + n3*ndotp_prmq1*u*(0.5*eta + 1) + n3_prmq1*ndotp*u*(0.5*eta + 1) + p3*u_prmq1*(-0.5*eta - 1)
    pprime1PN3_prmq2 = n3*ndotp*u_prmq2*(0.5*eta + 1) + n3*ndotp_prmq2*u*(0.5*eta + 1) + n3_prmq2*ndotp*u*(0.5*eta + 1) + p3*u_prmq2*(-0.5*eta - 1)
    pprime1PN3_prmq3 = n3*ndotp*u_prmq3*(0.5*eta + 1) + n3*ndotp_prmq3*u*(0.5*eta + 1) + n3_prmq3*ndotp*u*(0.5*eta + 1) + p3*u_prmq3*(-0.5*eta - 1)
    pprime1PN3_prmp1 = 0.5*eta*p3*psq_prmp1 + n3*ndotp_prmp1*u*(0.5*eta + 1)
    pprime1PN3_prmp2 = 0.5*eta*p3*psq_prmp2 + n3*ndotp_prmp2*u*(0.5*eta + 1)
    pprime1PN3_prmp3 = 0.5*eta*p3*psq_prmp3 + 0.5*eta*psq + n3*ndotp_prmp3*u*(0.5*eta + 1) + u*(-0.5*eta - 1)
    pprime1PN2_prmq1 = n2*ndotp*u_prmq1*(0.5*eta + 1) + n2*ndotp_prmq1*u*(0.5*eta + 1) + n2_prmq1*ndotp*u*(0.5*eta + 1) + p2*u_prmq1*(-0.5*eta - 1)
    pprime1PN2_prmq2 = n2*ndotp*u_prmq2*(0.5*eta + 1) + n2*ndotp_prmq2*u*(0.5*eta + 1) + n2_prmq2*ndotp*u*(0.5*eta + 1) + p2*u_prmq2*(-0.5*eta - 1)
    pprime1PN2_prmq3 = n2*ndotp*u_prmq3*(0.5*eta + 1) + n2*ndotp_prmq3*u*(0.5*eta + 1) + n2_prmq3*ndotp*u*(0.5*eta + 1) + p2*u_prmq3*(-0.5*eta - 1)
    pprime1PN2_prmp1 = 0.5*eta*p2*psq_prmp1 + n2*ndotp_prmp1*u*(0.5*eta + 1)
    pprime1PN2_prmp2 = 0.5*eta*p2*psq_prmp2 + 0.5*eta*psq + n2*ndotp_prmp2*u*(0.5*eta + 1) + u*(-0.5*eta - 1)
    pprime1PN2_prmp3 = 0.5*eta*p2*psq_prmp3 + n2*ndotp_prmp3*u*(0.5*eta + 1)
    pprime1PN1_prmq1 = n1*ndotp*u_prmq1*(0.5*eta + 1) + n1*ndotp_prmq1*u*(0.5*eta + 1) + n1_prmq1*ndotp*u*(0.5*eta + 1) + p1*u_prmq1*(-0.5*eta - 1)
    pprime1PN1_prmq2 = n1*ndotp*u_prmq2*(0.5*eta + 1) + n1*ndotp_prmq2*u*(0.5*eta + 1) + n1_prmq2*ndotp*u*(0.5*eta + 1) + p1*u_prmq2*(-0.5*eta - 1)
    pprime1PN1_prmq3 = n1*ndotp*u_prmq3*(0.5*eta + 1) + n1*ndotp_prmq3*u*(0.5*eta + 1) + n1_prmq3*ndotp*u*(0.5*eta + 1) + p1*u_prmq3*(-0.5*eta - 1)
    pprime1PN1_prmp1 = 0.5*eta*p1*psq_prmp1 + 0.5*eta*psq + n1*ndotp_prmp1*u*(0.5*eta + 1) + u*(-0.5*eta - 1)
    pprime1PN1_prmp2 = 0.5*eta*p1*psq_prmp2 + n1*ndotp_prmp2*u*(0.5*eta + 1)
    pprime1PN1_prmp3 = 0.5*eta*p1*psq_prmp3 + n1*ndotp_prmp3*u*(0.5*eta + 1)
    pprime3_prmq1 = pprime1PN3_prmq1 + pprime2PN3_prmq1
    pprime3_prmq2 = pprime1PN3_prmq2 + pprime2PN3_prmq2
    pprime3_prmq3 = pprime1PN3_prmq3 + pprime2PN3_prmq3
    pprime3_prmp1 = pprime1PN3_prmp1 + pprime2PN3_prmp1
    pprime3_prmp2 = pprime1PN3_prmp2 + pprime2PN3_prmp2
    pprime3_prmp3 = pprime1PN3_prmp3 + pprime2PN3_prmp3 + 1
    pprime2_prmq1 = pprime1PN2_prmq1 + pprime2PN2_prmq1
    pprime2_prmq2 = pprime1PN2_prmq2 + pprime2PN2_prmq2
    pprime2_prmq3 = pprime1PN2_prmq3 + pprime2PN2_prmq3
    pprime2_prmp1 = pprime1PN2_prmp1 + pprime2PN2_prmp1
    pprime2_prmp2 = pprime1PN2_prmp2 + pprime2PN2_prmp2 + 1
    pprime2_prmp3 = pprime1PN2_prmp3 + pprime2PN2_prmp3
    pprime1_prmq1 = pprime1PN1_prmq1 + pprime2PN1_prmq1
    pprime1_prmq2 = pprime1PN1_prmq2 + pprime2PN1_prmq2
    pprime1_prmq3 = pprime1PN1_prmq3 + pprime2PN1_prmq3
    pprime1_prmp1 = pprime1PN1_prmp1 + pprime2PN1_prmp1 + 1
    pprime1_prmp2 = pprime1PN1_prmp2 + pprime2PN1_prmp2
    pprime1_prmp3 = pprime1PN1_prmp3 + pprime2PN1_prmp3
    qprimesq_prmq1 = 2*qprime1*qprime1_prmq1 + 2*qprime2*qprime2_prmq1 + 2*qprime3*qprime3_prmq1
    qprimesq_prmq2 = 2*qprime1*qprime1_prmq2 + 2*qprime2*qprime2_prmq2 + 2*qprime3*qprime3_prmq2
    qprimesq_prmq3 = 2*qprime1*qprime1_prmq3 + 2*qprime2*qprime2_prmq3 + 2*qprime3*qprime3_prmq3
    qprimesq_prmp1 = 2*qprime1*qprime1_prmp1 + 2*qprime2*qprime2_prmp1 + 2*qprime3*qprime3_prmp1
    qprimesq_prmp2 = 2*qprime1*qprime1_prmp2 + 2*qprime2*qprime2_prmp2 + 2*qprime3*qprime3_prmp2
    qprimesq_prmp3 = 2*qprime1*qprime1_prmp3 + 2*qprime2*qprime2_prmp3 + 2*qprime3*qprime3_prmp3
    qprime_prmq1 = qprimesq_prmq1/(2*np.sqrt(qprimesq))
    qprime_prmq2 = qprimesq_prmq2/(2*np.sqrt(qprimesq))
    qprime_prmq3 = qprimesq_prmq3/(2*np.sqrt(qprimesq))
    qprime_prmp1 = qprimesq_prmp1/(2*np.sqrt(qprimesq))
    qprime_prmp2 = qprimesq_prmp2/(2*np.sqrt(qprimesq))
    qprime_prmp3 = qprimesq_prmp3/(2*np.sqrt(qprimesq))
    uprime_prmq1 = -qprime_prmq1/qprime**2
    uprime_prmq2 = -qprime_prmq2/qprime**2
    uprime_prmq3 = -qprime_prmq3/qprime**2
    uprime_prmp1 = -qprime_prmp1/qprime**2
    uprime_prmp2 = -qprime_prmp2/qprime**2
    uprime_prmp3 = -qprime_prmp3/qprime**2
    nprime3_prmq1 = qprime3_prmq1/qprime - qprime3*qprime_prmq1/qprime**2
    nprime3_prmq2 = qprime3_prmq2/qprime - qprime3*qprime_prmq2/qprime**2
    nprime3_prmq3 = qprime3_prmq3/qprime - qprime3*qprime_prmq3/qprime**2
    nprime3_prmp1 = qprime3_prmp1/qprime - qprime3*qprime_prmp1/qprime**2
    nprime3_prmp2 = qprime3_prmp2/qprime - qprime3*qprime_prmp2/qprime**2
    nprime3_prmp3 = qprime3_prmp3/qprime - qprime3*qprime_prmp3/qprime**2
    nprime2_prmq1 = qprime2_prmq1/qprime - qprime2*qprime_prmq1/qprime**2
    nprime2_prmq2 = qprime2_prmq2/qprime - qprime2*qprime_prmq2/qprime**2
    nprime2_prmq3 = qprime2_prmq3/qprime - qprime2*qprime_prmq3/qprime**2
    nprime2_prmp1 = qprime2_prmp1/qprime - qprime2*qprime_prmp1/qprime**2
    nprime2_prmp2 = qprime2_prmp2/qprime - qprime2*qprime_prmp2/qprime**2
    nprime2_prmp3 = qprime2_prmp3/qprime - qprime2*qprime_prmp3/qprime**2
    nprime1_prmq1 = qprime1_prmq1/qprime - qprime1*qprime_prmq1/qprime**2
    nprime1_prmq2 = qprime1_prmq2/qprime - qprime1*qprime_prmq2/qprime**2
    nprime1_prmq3 = qprime1_prmq3/qprime - qprime1*qprime_prmq3/qprime**2
    nprime1_prmp1 = qprime1_prmp1/qprime - qprime1*qprime_prmp1/qprime**2
    nprime1_prmp2 = qprime1_prmp2/qprime - qprime1*qprime_prmp2/qprime**2
    nprime1_prmp3 = qprime1_prmp3/qprime - qprime1*qprime_prmp3/qprime**2
    pprimedotnprime_prmq1 = nprime1*pprime1_prmq1 + nprime1_prmq1*pprime1 + nprime2*pprime2_prmq1 + nprime2_prmq1*pprime2 + nprime3*pprime3_prmq1 + nprime3_prmq1*pprime3
    pprimedotnprime_prmq2 = nprime1*pprime1_prmq2 + nprime1_prmq2*pprime1 + nprime2*pprime2_prmq2 + nprime2_prmq2*pprime2 + nprime3*pprime3_prmq2 + nprime3_prmq2*pprime3
    pprimedotnprime_prmq3 = nprime1*pprime1_prmq3 + nprime1_prmq3*pprime1 + nprime2*pprime2_prmq3 + nprime2_prmq3*pprime2 + nprime3*pprime3_prmq3 + nprime3_prmq3*pprime3
    pprimedotnprime_prmp1 = nprime1*pprime1_prmp1 + nprime1_prmp1*pprime1 + nprime2*pprime2_prmp1 + nprime2_prmp1*pprime2 + nprime3*pprime3_prmp1 + nprime3_prmp1*pprime3
    pprimedotnprime_prmp2 = nprime1*pprime1_prmp2 + nprime1_prmp2*pprime1 + nprime2*pprime2_prmp2 + nprime2_prmp2*pprime2 + nprime3*pprime3_prmp2 + nprime3_prmp2*pprime3
    pprimedotnprime_prmp3 = nprime1*pprime1_prmp3 + nprime1_prmp3*pprime1 + nprime2*pprime2_prmp3 + nprime2_prmp3*pprime2 + nprime3*pprime3_prmp3 + nprime3_prmp3*pprime3
    nprimedotpprimesq_prmq1 = 2*pprimedotnprime*pprimedotnprime_prmq1
    nprimedotpprimesq_prmq2 = 2*pprimedotnprime*pprimedotnprime_prmq2
    nprimedotpprimesq_prmq3 = 2*pprimedotnprime*pprimedotnprime_prmq3
    nprimedotpprimesq_prmp1 = 2*pprimedotnprime*pprimedotnprime_prmp1
    nprimedotpprimesq_prmp2 = 2*pprimedotnprime*pprimedotnprime_prmp2
    nprimedotpprimesq_prmp3 = 2*pprimedotnprime*pprimedotnprime_prmp3
    pprimesq_prmq1 = 2*pprime1*pprime1_prmq1 + 2*pprime2*pprime2_prmq1 + 2*pprime3*pprime3_prmq1
    pprimesq_prmq2 = 2*pprime1*pprime1_prmq2 + 2*pprime2*pprime2_prmq2 + 2*pprime3*pprime3_prmq2
    pprimesq_prmq3 = 2*pprime1*pprime1_prmq3 + 2*pprime2*pprime2_prmq3 + 2*pprime3*pprime3_prmq3
    pprimesq_prmp1 = 2*pprime1*pprime1_prmp1 + 2*pprime2*pprime2_prmp1 + 2*pprime3*pprime3_prmp1
    pprimesq_prmp2 = 2*pprime1*pprime1_prmp2 + 2*pprime2*pprime2_prmp2 + 2*pprime3*pprime3_prmp2
    pprimesq_prmp3 = 2*pprime1*pprime1_prmp3 + 2*pprime2*pprime2_prmp3 + 2*pprime3*pprime3_prmp3
    B_prmq1 = uprime*uprime_prmq1*(8 - 12*eta) + 2*uprime_prmq1
    B_prmq2 = uprime*uprime_prmq2*(8 - 12*eta) + 2*uprime_prmq2
    B_prmq3 = uprime*uprime_prmq3*(8 - 12*eta) + 2*uprime_prmq3
    B_prmp1 = uprime*uprime_prmp1*(8 - 12*eta) + 2*uprime_prmp1
    B_prmp2 = uprime*uprime_prmp2*(8 - 12*eta) + 2*uprime_prmp2
    B_prmp3 = uprime*uprime_prmp3*(8 - 12*eta) + 2*uprime_prmp3
    Binv_prmq1 = -B_prmq1/B**2
    Binv_prmq2 = -B_prmq2/B**2
    Binv_prmq3 = -B_prmq3/B**2
    Binv_prmp1 = -B_prmp1/B**2
    Binv_prmp2 = -B_prmp2/B**2
    Binv_prmp3 = -B_prmp3/B**2
    A_prmq1 = 6*eta*uprime**2*uprime_prmq1 - 2*uprime_prmq1
    A_prmq2 = 6*eta*uprime**2*uprime_prmq2 - 2*uprime_prmq2
    A_prmq3 = 6*eta*uprime**2*uprime_prmq3 - 2*uprime_prmq3
    A_prmp1 = 6*eta*uprime**2*uprime_prmp1 - 2*uprime_prmp1
    A_prmp2 = 6*eta*uprime**2*uprime_prmp2 - 2*uprime_prmp2
    A_prmp3 = 6*eta*uprime**2*uprime_prmp3 - 2*uprime_prmp3
    Heff_prmq1 = np.sqrt(A*(nprimedotpprimesq*(Binv - 1) + pprimesq + 1))*(A*(Binv_prmq1*nprimedotpprimesq + nprimedotpprimesq_prmq1*(Binv - 1) + pprimesq_prmq1)/2 + A_prmq1*(nprimedotpprimesq*(Binv - 1) + pprimesq + 1)/2)/(A*(nprimedotpprimesq*(Binv - 1) + pprimesq + 1))
    Heff_prmq2 = np.sqrt(A*(nprimedotpprimesq*(Binv - 1) + pprimesq + 1))*(A*(Binv_prmq2*nprimedotpprimesq + nprimedotpprimesq_prmq2*(Binv - 1) + pprimesq_prmq2)/2 + A_prmq2*(nprimedotpprimesq*(Binv - 1) + pprimesq + 1)/2)/(A*(nprimedotpprimesq*(Binv - 1) + pprimesq + 1))
    Heff_prmq3 = np.sqrt(A*(nprimedotpprimesq*(Binv - 1) + pprimesq + 1))*(A*(Binv_prmq3*nprimedotpprimesq + nprimedotpprimesq_prmq3*(Binv - 1) + pprimesq_prmq3)/2 + A_prmq3*(nprimedotpprimesq*(Binv - 1) + pprimesq + 1)/2)/(A*(nprimedotpprimesq*(Binv - 1) + pprimesq + 1))
    Heff_prmp1 = np.sqrt(A*(nprimedotpprimesq*(Binv - 1) + pprimesq + 1))*(A*(Binv_prmp1*nprimedotpprimesq + nprimedotpprimesq_prmp1*(Binv - 1) + pprimesq_prmp1)/2 + A_prmp1*(nprimedotpprimesq*(Binv - 1) + pprimesq + 1)/2)/(A*(nprimedotpprimesq*(Binv - 1) + pprimesq + 1))
    Heff_prmp2 = np.sqrt(A*(nprimedotpprimesq*(Binv - 1) + pprimesq + 1))*(A*(Binv_prmp2*nprimedotpprimesq + nprimedotpprimesq_prmp2*(Binv - 1) + pprimesq_prmp2)/2 + A_prmp2*(nprimedotpprimesq*(Binv - 1) + pprimesq + 1)/2)/(A*(nprimedotpprimesq*(Binv - 1) + pprimesq + 1))
    Heff_prmp3 = np.sqrt(A*(nprimedotpprimesq*(Binv - 1) + pprimesq + 1))*(A*(Binv_prmp3*nprimedotpprimesq + nprimedotpprimesq_prmp3*(Binv - 1) + pprimesq_prmp3)/2 + A_prmp3*(nprimedotpprimesq*(Binv - 1) + pprimesq + 1)/2)/(A*(nprimedotpprimesq*(Binv - 1) + pprimesq + 1))
    Hreal_prmq1 = Heff_prmq1/np.sqrt(2*eta*(Heff - 1) + 1)
    Hreal_prmq2 = Heff_prmq2/np.sqrt(2*eta*(Heff - 1) + 1)
    Hreal_prmq3 = Heff_prmq3/np.sqrt(2*eta*(Heff - 1) + 1)
    Hreal_prmp1 = Heff_prmp1/np.sqrt(2*eta*(Heff - 1) + 1)
    Hreal_prmp2 = Heff_prmp2/np.sqrt(2*eta*(Heff - 1) + 1)
    Hreal_prmp3 = Heff_prmp3/np.sqrt(2*eta*(Heff - 1) + 1)
    return np.array([Hreal_prmq1, Hreal_prmq2, Hreal_prmq3, Hreal_prmp1, Hreal_prmp2, Hreal_prmp3])

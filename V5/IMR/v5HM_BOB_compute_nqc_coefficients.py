import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from IMR.v5HM_BOB_unoptimized_nqc_terms import v5HM_BOB_unoptimized_nqc_terms
def v5HM_BOB_unoptimized_compute_nqc_coefficients(t_peak_strain, t_attach, h22_fine, dynamics_fine, h22NRpeak, omega22NRpeak, omegaQNM, tau):
    time = dynamics_fine[:,0]

    r = dynamics_fine[:,1]

    omega_orb = dynamics_fine[:,6]

    pr = dynamics_fine[:,3]

    rOmega = r * omega_orb

    q1 = pr*pr / (rOmega*rOmega)

    q2 = q1 / r

    q3 = q2 / np.sqrt(r)

    p1 = -pr / rOmega

    p2 = -p1 * pr * pr

    amplitude = np.abs(h22_fine[:,1])

    phase = np.unwrap(np.angle(h22_fine[:,1]))

    idx = np.argmin(np.abs(time - t_attach))

    N = 5

    left = np.max((0, idx - N))

    right = np.min((idx + N, len(time)))

    Q1 = q1 * amplitude

    Q2 = q2 * amplitude

    Q3 = q3 * amplitude

    intrp_Q1 = InterpolatedUnivariateSpline(time[left:right], Q1[left:right])

    intrp_Q2 = InterpolatedUnivariateSpline(time[left:right], Q2[left:right])

    intrp_Q3 = InterpolatedUnivariateSpline(time[left:right], Q3[left:right])

    Q = np.zeros([3,3])

    Q[:, 0] = intrp_Q1.derivatives(t_attach)[:-1]

    Q[:, 1] = intrp_Q2.derivatives(t_attach)[:-1]

    Q[:, 2] = intrp_Q3.derivatives(t_attach)[:-1]

    intrp_P1 = InterpolatedUnivariateSpline(time[left:right], p1[left:right])

    intrp_P2 = InterpolatedUnivariateSpline(time[left:right], p2[left:right])

    P = np.zeros([2,2])

    P[:, 0] = -intrp_p1.derivatives(t_attach)[1:-1]

    P[:, 1] = -intrp_p2.derivatives(t_attach)[1:-1]

    intrp_amp = InterpolatedUnivariateSpline(time[left:right], amplitude[left:right])

    amps_inspiral = intrp_amp.derivatives(t_attach)[:-1]

    intrp_phase = InterpolatedUnivariateSpline(time[left:right], phase[left:right])

    omega, omegaDot = intrp_phase.derivatives(t_attach)[1:-1]

    if omega * omegaDot > 0.0:

        omega = np.abs(omega)

        omegaDot = np.abs(omegaDot)

    else:

        omega = np.abs(omega)

        omegaDot = -np.abs(omegaDot)

    omegas_inspiral = [omega,omegaDot]

    amps_BOB, omegas_BOB = v5HM_BOB_unoptimized_nqc_terms(t_attach,t_peak_strain,h22NRpeak,omega22NRpeak,omegaQNM,tau)

    rhs_ai = np.array(amps_BOB) - np.array(amps_inspiral)

    rhs_bi = np.array(omegas_BOB) - np.array(omegas_inspiral)

    res_ai = np.linalg.solve(Q,amps)

    res_bi = np.linalg.solve(P,omegas)

    return {"a": res_ai, "b" : res_bi}
import Derivatives.pyseobnr_derivatives as H
import Radiation.pyseobnr_waveforms as RR

def get_rhs(t,z,chi_1,chi_2,m_1,m_2,verbose = False):
    """
    Compute the RHS of the EOB evolution equations.
    In particular this function returns
    :math:`\\dot{r}`, :math:`\\dot{\phi}` , :math:`\\dot{p}_{r}` and :math:`\\dot{p}_{\\phi}` .

    See for example Eq(2) of [Buades2021]_ .

    The Hamiltonian is given by Eq(14) in [SEOBNRv5HM]_ doc
    and explicitly spelled out in Section I.C of [SEOBNRv5HM-theory]_
    and the RR force is described in Eq(43) of [SEOBNRv5HM]_ document, both
    contained in [DCC_T2300060]_ .

    """
    q = z[:2]
    p = z[2:]

    dynamics = H.dynamics(q, p, chi_1, chi_2, m_1, m_2)
    H_val = dynamics[4]
    omega = dynamics[3]
    pcirc = [0,p[1]]
    omega_circ = H.omega(q, pcirc, chi_1, chi_2, m_1, m_2)

    xi = dynamics[5]

    RR_f = RR.RR_force(m_1,m_2,q, p, omega, omega_circ, H_val,chi_1,chi_2)
    deriv = [xi * dynamics[2], dynamics[3], -dynamics[0] * xi + RR_f[0], -dynamics[1] + RR_f[1]]
    if not verbose:
        return deriv
    else:
        return deriv, RR_f[1], omega, omega_circ, xi, dynamics[0], m_1*m_2/(m_1 + m_2)/(m_1 + m_2), p[0],p[1]
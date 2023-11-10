from numpy import zeros,array,pi,sqrt,absolute,exp,log
from scipy.special import factorial2

euler_gamma=0.5772156649015329
I = 1j
ell_max = 8
PN_limit = 11

ylms = [
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [
        0.0,
        0.3454941494713355,
        0.0,
        0.3231801841141506,
        4.055554344078009e-17,
        0.32028164857621516,
        8.000317029649286e-17,
        0.31937046138540076,
        1.1920022960993937e-16,
    ],
    [
        0.0,
        0.0,
        0.3862742020231896,
        0.0,
        0.33452327177864466,
        4.447491653287784e-17,
        0.32569524293385776,
        8.531070568407408e-17,
        0.32254835519288305,
    ],
    [
        0.0,
        0.0,
        0.0,
        0.4172238236327842,
        0.0,
        0.34594371914684025,
        4.8749755603568385e-17,
        0.331899519333737,
        9.130149959227967e-17,
    ],
    [
        0.0,
        0.0,
        0.0,
        0.0,
        0.4425326924449826,
        0.0,
        0.3567812628539981,
        5.310595586255289e-17,
        0.3382915688890245,
    ],
    [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.46413220344085826,
        0.0,
        0.3669287245764378,
        5.745136532071183e-17,
    ],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.48308411358006625, 0.0, 0.3764161087284946],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5000395635705508, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5154289843972844],
]

LOOKUP_TABLE = [
        1,
        1,
        2,
        6,
        24,
        120,
        720,
        5040,
        40320,
        362880,
        3628800,
        39916800,
        479001600,
        6227020800,
        87178291200,
        1307674368000,
        20922789888000,
        355687428096000,
        6402373705728000,
        121645100408832000,
        2432902008176640000,
    ]

def calculate_multipole_prefix(m1, m2, l, m):
    """
    Calculates the Newtonian multipole prefactors, see Eq. 25-27 in https://dcc.ligo.org/LIGO-T2300060 (SEOBNRv5HM.pdf).
    See also Sec. 2A of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.79.064004
    """
    n = 0.0
    

    totalMass = m1 + m2

    epsilon = (l + m) % 2

    x1 = m1 / totalMass
    x2 = m2 / totalMass

    eta = m1 * m2 / (totalMass * totalMass)
    
    if absolute(m % 2) == 0:
        sign = 1
    else:
        sign = -1

    #
    # Eq. 7 of Damour, Iyer and Nagar 2008.
    # For odd m, c is proportional to dM = m1-m2. In the equal-mass case, c = dM = 0.
    # In the equal-mass unequal-spin case, however, when spins are different, the odd m term is generally not zero.
    # In this case, c can be written as c0 * dM, while spins terms in PN expansion may take the form chiA/dM.
    # Although the dM's cancel analytically, we can not implement c and chiA/dM with the possibility of dM -> 0.
    # Therefore, for this case, we give numerical values of c0 for relevant modes, and c0 is calculated as
    # c / dM in the limit of dM -> 0. Consistently, for this case, we implement chiA instead of chiA/dM
    # below.
    # Note that for equal masses and odd m modes we currently only have non-zero c up to l=5. If new spinning terms
    # are added to modes above l=5 this will need to be revisited.
    if m1 != m2 or sign == 1:
        c = pow(x2, l + epsilon - 1) + sign * pow(x1, l + epsilon - 1)

    else:
        if l == 2:
            c = -1.0
        elif l == 3:
            c = -1.0
        elif l == 4:
            c = -0.5
        elif l == 5:
            c = -0.5
        else:
            c = 0.0

    # Eqs 5 and 6. Dependent on the value of epsilon (parity), we get different n
    if epsilon == 0:
        n = 1j * m
        n = n ** l

        mult1 = 8.0 * pi / factorial2(2 * l + 1)
        mult2 = ((l + 1) * (l + 2)) / (l * (l - 1))
        mult2 = sqrt(mult2)

        n *= mult1
        n *= mult2

    elif epsilon == 1:

        n = 1j * m
        n = n ** l
        n = -n

        mult1 = 16.0 * pi / factorial2(2 * l + 1)

        mult2 = (2 * l + 1) * (l + 2) * (l * l - m * m)
        mult2 /= (2 * l - 1) * (l + 1) * l * (l - 1)
        mult2 = sqrt(mult2)

        n *= 1j * mult1
        n *= mult2

    prefix = n * eta * c
    return prefix

def compute_newtonian_prefixes(m1, m2):
    """
    Loop to set the Newtonian multipole prefactors, see Eq. 25-27 in https://dcc.ligo.org/LIGO-T2300060 (SEOBNRv5HM.pdf).
    """
    prefixes = zeros([9,9],dtype = complex)
    for l in range(2, ell_max + 1):
        for m in range(1, l + 1):
            prefixes[l][m] = calculate_multipole_prefix(m1, m2, l, m)
    return prefixes



def compute_tail(omega, H):
    """
    Calculate the resummed Tail effects, see Eq. 32 in https://dcc.ligo.org/LIGO-T2300060 (SEOBNRv5HM.pdf).
    See also Sec. 2B of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.79.064004
    and Eq. (42) of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.87.084035
    """
    Tlm = zeros([ell_max+1,ell_max+1])
    tlmprod_fac = 1.0
    for m in range(1,ell_max+1):
        k = m*omega
        hathatk = H*k
        hathatksq4 = 4.0 * hathatk * hathatk
        hathatk4pi = 4.0 * pi * hathatk
        tlmprod_fac = 1.0
        tlmprefac = sqrt(hathatk4pi / (1.0 - exp(-hathatk4pi)))

        for j in range(1,ell_max+1):
            z2 = LOOKUP_TABLE[j]
            tlmprod_fac*= hathatksq4 + j**2
            if m>j:
                continue
            Tlm[j][m] = tlmprefac*sqrt(tlmprod_fac)/z2
    return Tlm

def compute_rho_coeffs(nu,dm, a,chiS,chiA, extra_PN_terms):

    """
    Compute the amplitude residual coefficients.
    See Sec. 2C and 2D of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.79.064004
    See Eq. 59-63 of https://dcc.ligo.org/LIGO-T2300060 (SEOBNRv5_theory.pdf) for new terms, rest copied from SEOBNRv4HM LAL code.

    Coefficients can be found in:
    - https://dcc.ligo.org/LIGO-T2300060 (SEOBNRv5_theory.pdf and SEOBNRv5HM.pdf)
    - Appendix A of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.98.084028
    - Eqs. (2.4) to (2.6) of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.95.044028
    - https://journals.aps.org/prd/pdf/10.1103/PhysRevD.83.064003
    - https://journals.aps.org/prd/pdf/10.1103/PhysRevD.86.024011
    """
    rho_coeffs = zeros([ell_max+1,ell_max+1,PN_limit])
    rho_coeffs_log = zeros([ell_max+1,ell_max+1,PN_limit])
    f_coeffs = zeros([ell_max+1,ell_max+1,PN_limit])
    f_coeffs_vh = zeros([ell_max+1,ell_max+1,PN_limit])
    nu2 = nu*nu
    nu3 = nu*nu2

    a2 = a*a
    a3 = a*a2
    dm2 = dm*dm
    atemp = a

    chiA2 = chiA * chiA
    chiS2 = chiS * chiS
    chiA3 = chiA2 * chiA
    chiS3 = chiS2 * chiS

    m1Plus3nu = -1.0 + 3.0 * nu
    m1Plus3nu2 = m1Plus3nu * m1Plus3nu
    m1Plus3nu3 = m1Plus3nu * m1Plus3nu2

    # (2,2) mode begins
    rho_coeffs[2,2][2] = -43.0 / 42.0 + (55.0 * nu) / 84.0
    rho_coeffs[2,2][3]  = (-2.0 * (chiS + chiA * dm - chiS * nu)) / 3.0

    rho_coeffs[2,2][4]  = (
        -20555.0 / 10584.0
        + 0.5 * (chiS + chiA * dm) * (chiS + chiA * dm)
        - (33025.0 * nu) / 21168.0
        + (19583.0 * nu2) / 42336.0
    )

    rho_coeffs[2,2][5]  = (
        -34.0 / 21.0 + 49.0 * nu / 18.0 + 209.0 * nu2 / 126.0
    ) * chiS + (-34.0 / 21.0 - 19.0 * nu / 42.0) * dm * chiA

    rho_coeffs[2,2][6]  = (
        1556919113.0 / 122245200.0
        + (89.0 * a2) / 252.0
        - (48993925.0 * nu) / 9779616.0
        - (6292061.0 * nu2) / 3259872.0
        + (10620745.0 * nu3) / 39118464.0
        + (41.0 * nu * pi * pi) / 192.0
    )
    rho_coeffs_log[2,2][6] = -428.0 / 105.0

    # See https://dcc.ligo.org/T1600383
    rho_coeffs[2,2][7]  = (
        a3 / 3.0
        + chiA
        * dm
        * (18733.0 / 15876.0 + (50140.0 * nu) / 3969.0 + (97865.0 * nu2) / 63504.0)
        + chiS
        * (
            18733.0 / 15876.0
            + (74749.0 * nu) / 5292.0
            - (245717.0 * nu2) / 63504.0
            + (50803.0 * nu3) / 63504.0
        )
    )

    rho_coeffs[2,2][8]  = (
        -387216563023.0 / 160190110080.0 +
            (18353.0 * a2) / 21168.0 - a2 * a2 / 8.0
    )

    rho_coeffs_log[2,2][8] = 9202.0 / 2205.0
    rho_coeffs[2,2][10]  = -16094530514677.0 / 533967033600.0
    rho_coeffs_log[2,2][10] = 439877.0 / 55566.0
    # (2,2) mode ends

    # We set test-spin terms to 0 (as in SEOBNRv4HM)
    # as no major improvement was found when trying to include them
    a = 0.0
    a2 = 0.0
    a3 = 0.0
    # (2,1) mode begins
    if dm2:
        rho_coeffs[2,1][1] = 0.0
        rho_coeffs[2,1][2] = -59.0 / 56 + (23.0 * nu) / 84.0
        rho_coeffs[2,1][3] = 0.0

        rho_coeffs[2,1][4] = (
            -47009.0 / 56448.0
            - (865.0 * a2) / 1792.0
            - (405.0 * a2 * a2) / 2048.0
            - (10993.0 * nu) / 14112.0
            + (617.0 * nu2) / 4704.0
        )
        rho_coeffs[2,1][5] = (
            (-98635.0 * a) / 75264.0
            + (2031.0 * a * a2) / 7168.0
            - (1701.0 * a2 * a3) / 8192.0
        )
        rho_coeffs[2,1][6] = (
            7613184941.0 / 2607897600.0
            + (9032393.0 * a2) / 1806336.0
            + (3897.0 * a2 * a2) / 16384.0
            - (15309.0 * a3 * a3) / 65536.0
        )
        rho_coeffs_log[2,1][6] = -107.0 / 105.0
        rho_coeffs[2,1][7] = (
            (-3859374457.0 * a) / 1159065600.0
            - (55169.0 * a3) / 16384.0
            + (18603.0 * a2 * a3) / 65536.0
            - (72171.0 * a2 * a2 * a3) / 262144.0
        )
        rho_coeffs_log[2,1][7] = 107.0 * a / 140.0
        rho_coeffs[2,1][8] = -1168617463883.0 / 911303737344.0
        rho_coeffs_log[2,1][8] = 6313.0 / 5880.0
        rho_coeffs[2,1][10] = -63735873771463.0 / 16569158860800.0
        rho_coeffs_log[2,1][10] = 5029963.0 / 5927040.0

        f_coeffs[2,1][1] = (-3.0 * (chiS + chiA / dm)) / (2.0)

        f_coeffs[2,1][3] = (
            (
                chiS * dm * (427.0 + 79.0 * nu)
                + chiA * (147.0 + 280.0 * dm * dm + 1251.0 * nu)
            )
            / 84.0
            / dm
        )
        # RC: New terms for SEOBNRv4HM, they are put to zero if use_hm == 0


        # RC: This terms are in Eq.A11 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
        f_coeffs[2,1][4] = (
            (-3.0 - 2.0 * nu) * chiA2
            + (-3.0 + nu / 2.0) * chiS2
            + (-6.0 + 21.0 * nu / 2.0) * chiS * chiA / dm
        )
        f_coeffs[2,1][5] = (
            (3.0 / 4.0 - 3.0 * nu) * chiA3 / dm
            + (
                -81.0 / 16.0
                + 1709.0 * nu / 1008.0
                + 613.0 * nu2 / 1008.0
                + (9.0 / 4.0 - 3 * nu) * chiA2
            )
            * chiS
            + 3.0 / 4.0 * chiS3
            + (
                -81.0 / 16.0
                - 703.0 * nu2 / 112.0
                + 8797.0 * nu / 1008.0
                + (9.0 / 4.0 - 6.0 * nu) * chiS2
            )
            * chiA
            / dm
        )
        '''
        This was in SEOBNRv4HM
        f_coeffs[2,1][6] = (
            (4163.0 / 252.0 - 9287.0 * nu /
                1008.0 - 85.0 * nu2 / 112.0) * chiA2
            + (4163.0 / 252.0 - 2633.0 * nu / 1008.0 + 461.0 * nu2 / 1008.0)
            * chiS2
            + (4163.0 / 126.0 - 1636.0 * nu / 21.0 + 1088.0 * nu2 / 63.0)
            * chiS
            * chiA
            / dm
        )
        '''
        f_coeffs[2,1][6] = (((16652 - 9287*nu + 720*nu**2)*chiA**2)/1008 +
                            ((16652 - 39264*nu + 9487*nu**2)*chiA*chiS)/(504*dm) +
                            ((16652 - 2633*nu + 1946*nu**2)*chiS**2)/1008)
    else:
        f_coeffs[2,1][1] = -3.0 * chiA / 2.0
        f_coeffs[2,1][3] = (
            chiS * dm * (427.0 + 79.0 * nu)
            + chiA * (147.0 + 280.0 * dm * dm + 1251.0 * nu)
        ) / 84.0
        # New terms for SEOBNRv4HM, they are put to zero if use_hm == 0

        # RC: This terms are in Eq.A11 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
        f_coeffs[2,1][4] = (-6 + 21 * nu / 2.0) * chiS * chiA
        f_coeffs[2,1][5] = (3.0 / 4.0 - 3.0 * nu) * chiA3 + (
            -81.0 / 16.0
            - 703.0 * nu2 / 112.0
            + 8797.0 * nu / 1008.0
            + (9.0 / 4.0 - 6.0 * nu) * chiS2
        ) * chiA
        '''
        This was in SEOBNRv4HM
        f_coeffs[2,1][6] = (
            (4163.0 / 126.0 - 1636.0 * nu / 21.0 + 1088.0 * nu2 / 63.0)
            * chiS
            * chiA
        )
        '''
        f_coeffs[2,1][6] = ((16652 - 39264*nu + 9487*nu**2)*chiA*chiS)/504

    # (2,1) mode ends

    # (3,3) mode begins

    if dm2:
        rho_coeffs[3,3][2] = -7.0 / 6.0 + (2.0 * nu) / 3.0
        rho_coeffs[3,3][3] = 0.0
        rho_coeffs[3,3][4] = (
            -6719.0 / 3960.0
            + a2 / 2.0
            - (1861.0 * nu) / 990.0
            + (149.0 * nu2) / 330.0
        )
        rho_coeffs[3,3][5] = (-4.0 * a) / 3.0
        rho_coeffs[3,3][6] = (
                3203101567.0 / 227026800.0
                + (5.0 * a2) / 36.0
                + (-129509.0 / 25740.0 + 41.0 / 192.0 * pi * pi) * nu
                - 274621.0 / 154440.0 * nu2
                + 12011.0 / 46332.0 * nu3
            )
        rho_coeffs_log[3,3][6] = -26.0 / 7.0
        rho_coeffs[3,3][7] = (5297.0 * a) / 2970.0 + a * a2 / 3.0
        rho_coeffs[3,3][8] = -57566572157.0 / 8562153600.0
        rho_coeffs_log[3,3][8] = 13.0 / 3.0


        rho_coeffs[3,3][10] = -903823148417327.0 / 30566888352000.0
        rho_coeffs_log[3,3][10] = 87347.0 / 13860.0

        f_coeffs[3,3][3]= (
            chiS * dm * (-4.0 + 5.0 * nu) + chiA * (-4.0 + 19.0 * nu)
        ) / (2.0 * dm)

        f_coeffs[3,3][4]= (
            3.0 / 2.0 * chiS2 * dm
            + (3.0 - 12 * nu) * chiA * chiS
            + dm * (3.0 / 2.0 - 6.0 * nu) * chiA2
        ) / (dm)
        f_coeffs[3,3][5] = (
            dm * (241.0 / 30.0 * nu2 + 11.0 /
                    20.0 * nu + 2.0 / 3.0) * chiS
            + (407.0 / 30.0 * nu2 - 593.0 / 60.0 * nu + 2.0 / 3.0) * chiA
        ) / (dm)
        f_coeffs[3,3][6] = (
            dm * (6.0 * nu2 - 27.0 / 2.0 * nu - 7.0 / 4.0) * chiS2
            + (44.0 * nu2 - 1.0 * nu - 7.0 / 2.0) * chiA * chiS
            + dm * (-12 * nu2 + 11.0 / 2.0 * nu - 7.0 / 4.0) * chiA2
        ) / dm
        f_coeffs_vh[3,3][6] =  (
                dm * (593.0 / 108.0 * nu - 81.0 / 20.0) * chiS
                + (7339.0 / 540.0 * nu - 81.0 / 20.0) * chiA
            ) / (dm)

    else:
        f_coeffs[3,3][3] = chiA * 3.0 / 8.0

        # RC: This terms are in Eq.A10 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
        f_coeffs[3,3][4] = (3.0 - 12 * nu) * chiA * chiS
        f_coeffs[3,3][5] = (
            407.0 / 30.0 * nu2 - 593.0 / 60.0 * nu + 2.0 / 3.0
        ) * chiA
        f_coeffs[3,3][6] = (44.0 * nu2 - 1.0 * nu - 7.0 / 2.0) * chiA * chiS
        f_coeffs_vh[3,3][6] = (7339.0 / 540.0 * nu - 81.0 / 20.0) * chiA

    # (3,3) mode ends

    # (3,2) mode begins
    rho_coeffs[3,2][1] = (4.0 * chiS * nu) / (-3.0 * m1Plus3nu)
    rho_coeffs[3,2][2] = (328.0 - 1115.0 * nu + 320.0 *
                      nu2) / (270.0 * m1Plus3nu)

    rho_coeffs[3,2][3] = (
        2.0
        * (
            45.0 * a * m1Plus3nu3
            - a
            * nu
            * (328.0 - 2099.0 * nu + 5.0 * (733.0 + 20.0 * a2) * nu2 - 960.0 * nu3)
        )
    ) / (405.0 * m1Plus3nu3)

    rho_coeffs[3,2][4] = a2 / 3.0 + (
        -1444528.0
        + 8050045.0 * nu
        - 4725605.0 * nu2
        - 20338960.0 * nu3
        + 3085640.0 * nu2 * nu2
    ) / (1603800.0 * m1Plus3nu2)
    rho_coeffs[3,2][5] = (-2788.0 * a) / 1215.0
    rho_coeffs[3,2][6] = 5849948554.0 / 940355325.0 + (488.0 * a2) / 405.0
    rho_coeffs_log[3,2][6] = -104.0 / 63.0
    rho_coeffs[3,2][8] = -10607269449358.0 / 3072140846775.0
    rho_coeffs_log[3,2][8] = 17056.0 / 8505.0
    # (3,2) mode ends

    # (3,1) mode begins
    if dm2:

        rho_coeffs[3,1][2] = -13.0 / 18.0 - (2.0 * nu) / 9.0
        rho_coeffs[3,1][3] = 0.0
        rho_coeffs[3,1][4] = (
            101.0 / 7128.0
            - (5.0 * a2) / 6.0
            - (1685.0 * nu) / 1782.0
            - (829.0 * nu2) / 1782.0
        )
        rho_coeffs[3,1][5] = (4.0 * a) / 9.0
        rho_coeffs[3,1][6] = 11706720301.0 / 6129723600.0 - (49.0 * a2) / 108.0
        rho_coeffs_log[3,1][6] = -26.0 / 63.0
        rho_coeffs[3,1][7] = (-2579.0 * a) / 5346.0 + a * a2 / 9.0
        rho_coeffs[3,1][8] = 2606097992581.0 / 4854741091200.0
        rho_coeffs_log[3,1][8] = 169.0 / 567.0

        f_coeffs[3,1][3] = (
            chiA * (-4.0 + 11.0 * nu) + chiS * dm * (-4.0 + 13.0 * nu)
        ) / (2.0 * dm)

    else:
        f_coeffs[3,1][3]= -chiA * 5.0 / 8.0
    # (3,1) mode ends

    # (4,4) mode begins
    rho_coeffs[4,4][2] = (1614.0 - 5870.0 * nu + 2625.0 *
                      nu2) / (1320.0 * m1Plus3nu)
    rho_coeffs[4,4][3] = (
        chiA * (10.0 - 39.0 * nu) * dm + chiS *
                (10.0 - 41.0 * nu + 42.0 * nu2)
    ) / (15.0 * m1Plus3nu)

    # RC: This terms are in Eq.A8 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
    rho_coeffs[4,4][4] = (
        (
            -511573572.0
            + 2338945704.0 * nu
            - 313857376.0 * nu2
            - 6733146000.0 * nu3
            + 1252563795.0 * nu2 * nu2
        )
        / (317116800.0 * m1Plus3nu2)
        + chiS2 / 2.0
        + dm * chiS * chiA
        + dm2 * chiA2 / 2.0
    )
    rho_coeffs[4,4][5] = chiA * dm * (
        -8280.0 + 42716.0 * nu - 57990.0 * nu2 + 8955 * nu3
    ) / (6600.0 * m1Plus3nu2) + chiS * (
        -8280.0
        + 66284.0 * nu
        - 176418.0 * nu2
        + 128085.0 * nu3
        + 88650 * nu2 * nu2
    ) / (
        6600.0 * m1Plus3nu2
    )
    rho_coeffs[4,4][8] = -172066910136202271.0 / 19426955708160000.0
    rho_coeffs_log[4,4][8] = 845198.0 / 190575.0
    rho_coeffs[4,4][10] = -17154485653213713419357.0 / 568432724020761600000.0
    rho_coeffs_log[4,4][10] = 22324502267.0 / 3815311500.0

    rho_coeffs[4,4][6] = 16600939332793.0 / 1098809712000.0 + (217.0 * a2) / 3960.0
    rho_coeffs_log[4,4][6] = -12568.0 / 3465.0
    # (4,4) mode ends
    # (4,3) mode begins
    if dm2:
        rho_coeffs[4,3][1] = 0.0
        rho_coeffs[4,3][2] = (222.0 - 547.0 * nu + 160.0 * nu2) / (
            176.0 * (-1.0 + 2.0 * nu)
        )
        rho_coeffs[4,3][4] = -6894273.0 / 7047040.0 + (3.0 * a2) / 8.0
        rho_coeffs[4,3][5] = (-12113.0 * a) / 6160.0
        rho_coeffs[4,3][6] = 1664224207351.0 / 195343948800.0
        rho_coeffs_log[4,3][6] = -1571.0 / 770.0
        f_coeffs[4,3][1] = (5.0 * (chiA - chiS * dm) * nu) / (
            2.0 * dm * (-1.0 + 2.0 * nu)
        )

    else:
        f_coeffs[4,3][1] = -5.0 * chiA / 4.0
    # (4,3) mode ends

    # (4,2) mode begins
    rho_coeffs[4,2][2] = (1146.0 - 3530.0 * nu + 285.0 *
                      nu2) / (1320.0 * m1Plus3nu)
    rho_coeffs[4,2][3] = (
        chiA * (10.0 - 21.0 * nu) * dm + chiS *
                (10.0 - 59.0 * nu + 78.0 * nu2)
    ) / (15.0 * m1Plus3nu)
    rho_coeffs[4,2][4] = a2 / 2.0 + (
        -114859044.0
        + 295834536.0 * nu
        + 1204388696.0 * nu2
        - 3047981160.0 * nu3
        - 379526805.0 * nu2 * nu2
    ) / (317116800.0 * m1Plus3nu2)
    rho_coeffs[4,2][5] = (-7.0 * a) / 110.0
    rho_coeffs[4,2][6] = 848238724511.0 / 219761942400.0 + (2323.0 * a2) / 3960.0
    rho_coeffs_log[4,2][6] = -3142.0 / 3465.0
    # (4,2) mode ends

    # (4,1) mode begins
    if dm2:

        rho_coeffs[4,1][1] = 0.0
        rho_coeffs[4,1][2] = (602.0 - 1385.0 * nu + 288.0 * nu2) / (
            528.0 * (-1.0 + 2.0 * nu)
        )
        rho_coeffs[4,1][4] = -7775491.0 / 21141120.0 + (3.0 * a2) / 8.0
        rho_coeffs[4,1][5] = (-20033.0 * a) / 55440.0 - (5 * a * a2) / 6.0
        rho_coeffs[4,1][6] = 1227423222031.0 / 1758095539200.0
        rho_coeffs_log[4,1][6] = -1571.0 / 6930.0
        f_coeffs[4,1][1] = (5.0 * (chiA - chiS * dm) * nu) / (
            2.0 * dm * (-1.0 + 2.0 * nu)
        )

    else:
        f_coeffs[4,1][1] = -5.0 * chiA / 4.0

    # (4,1) mode ends

    # (5,5) mode begins
    if dm2:

        rho_coeffs[5,5][2] = (487.0 - 1298.0 * nu + 512.0 * nu2) / (
            390.0 * (-1.0 + 2.0 * nu)
        )
        rho_coeffs[5,5][3] = (-2.0 * a) / 3.0
        rho_coeffs[5,5][4] = -3353747.0 / 2129400.0 + a2 / 2.0
        rho_coeffs[5,5][5] = -241.0 * a / 195.0



        # RC: This terms are in Eq.A9 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
        rho_coeffs[5,5][6] = 190606537999247.0 / 11957879934000.0
        rho_coeffs_log[5,5][6] = -1546.0 / 429.0
        rho_coeffs[5,5][8] = -1213641959949291437.0 / 118143853747920000.0
        rho_coeffs_log[5,5][8] = 376451.0 / 83655.0
        rho_coeffs[5,5][10] = -150082616449726042201261.0 / 4837990810977324000000.0
        rho_coeffs_log[5,5][10] = 2592446431.0 / 456756300.0

        f_coeffs[5,5][3] = chiA / dm * (
            10.0 / (3.0 * (-1.0 + 2.0 * nu))
            - 70.0 * nu / (3.0 * (-1.0 + 2.0 * nu))
            + 110.0 * nu2 / (3.0 * (-1.0 + 2.0 * nu))
        ) + chiS * (
            10.0 / (3.0 * (-1.0 + 2.0 * nu))
            - 10.0 * nu / (-1.0 + 2.0 * nu)
            + 10 * nu2 / (-1.0 + 2.0 * nu)
        )
        f_coeffs[5,5][4] = (
            chiS2
            * (-5.0 / (2.0 * (-1.0 + 2.0 * nu)) + 5.0 * nu / (-1.0 + 2.0 * nu))
            + chiA
            * chiS
            / dm
            * (
                -5.0 / (-1.0 + 2.0 * nu)
                + 30.0 * nu / (-1.0 + 2.0 * nu)
                - 40.0 * nu2 / (-1.0 + 2.0 * nu)
            )
            + chiA2
            * (
                -5.0 / (2.0 * (-1.0 + 2.0 * nu))
                + 15.0 * nu / (-1.0 + 2.0 * nu)
                - 20.0 * nu2 / (-1.0 + 2.0 * nu)
            )
        )

    else:

        # RC: This terms are in Eq.A12 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
        f_coeffs[5,5][3] = chiA * (
            10.0 / (3.0 * (-1.0 + 2.0 * nu))
            - 70.0 * nu / (3.0 * (-1.0 + 2.0 * nu))
            + 110.0 * nu2 / (3.0 * (-1.0 + 2.0 * nu))
        )
        f_coeffs[5,5][4] = (
            chiA
            * chiS
            * (
                -5.0 / (-1.0 + 2.0 * nu)
                + 30.0 * nu / (-1.0 + 2.0 * nu)
                - 40.0 * nu2 / (-1.0 + 2.0 * nu)
            )
        )
    # (5,5) mode ends

    # (5,4) mode begins
    rho_coeffs[5,4][2] = (
        -17448.0 + 96019.0 * nu - 127610.0 * nu2 + 33320.0 * nu3
    ) / (13650.0 * (1.0 - 5.0 * nu + 5.0 * nu2))
    rho_coeffs[5,4][3] = (-2.0 * a) / 15.0
    rho_coeffs[5,4][4] = -16213384.0 / 15526875.0 + (2.0 * a2) / 5.0
    # (5,4) mode ends

    # (5,3)  mode begins
    if dm2:

        rho_coeffs[5,3][2] = (375.0 - 850.0 * nu + 176.0 * nu2) / (
            390.0 * (-1.0 + 2.0 * nu)
        )
        rho_coeffs[5,3][3] = (-2.0 * a) / 3.0
        rho_coeffs[5,3][4] = -410833.0 / 709800.0 + a2 / 2.0
        rho_coeffs[5,3][5] = -103.0 * a / 325.0
    # (5,3) mode ends

    # (5,2) mode begins
    rho_coeffs[5,2][2] = (
        -15828.0 + 84679.0 * nu - 104930.0 * nu2 + 21980.0 * nu3
    ) / (13650.0 * (1.0 - 5.0 * nu + 5.0 * nu2))
    rho_coeffs[5,2][3] = (-2.0 * a) / 15.0
    rho_coeffs[5,2][4] = -7187914.0 / 15526875.0 + (2.0 * a2) / 5.0
    # (5,2) mode ends

    # (5,1) mode begins
    if dm2:
        rho_coeffs[5,1][2] = (319.0 - 626.0 * nu + 8.0 * nu2) / (
            390.0 * (-1.0 + 2.0 * nu)
        )
        rho_coeffs[5,1][3] = (-2.0 * a) / 3.0
        rho_coeffs[5,1][4] = -31877.0 / 304200.0 + a2 / 2.0
        rho_coeffs[5,1][5] = 139.0 * a / 975.0
    # (5,1) mode ends

    # (6,6) mode begins
    rho_coeffs[6,6][2] = (-106.0 + 602.0 * nu - 861.0 * nu2 + 273.0 * nu3) / (
        84.0 * (1.0 - 5.0 * nu + 5.0 * nu2)
    )
    rho_coeffs[6,6][3] = (-2.0 * a) / 3.0
    rho_coeffs[6,6][4] = -1025435.0 / 659736.0 + a2 / 2.0
    # (6,6) mode ends

    # (6,5) mode begins
    if dm2:
        rho_coeffs[6,5][2] = (-185.0 + 838.0 * nu - 910.0 * nu2 + 220.0 * nu3) / (
            144.0 * (dm2 + 3.0 * nu2)
        )
        rho_coeffs[6,5][3] = -2.0 * a / 9.0
    # (6,5) mode ends

    # (6,4) mode begins
    rho_coeffs[6,4][2] = (-86.0 + 462.0 * nu - 581.0 * nu2 + 133.0 * nu3) / (
        84.0 * (1.0 - 5.0 * nu + 5.0 * nu2)
    )
    rho_coeffs[6,4][3] = (-2.0 * a) / 3.0
    rho_coeffs[6,4][4] = -476887.0 / 659736.0 + a2 / 2.0
    # (6,4) mode ends

    # (6,3) mode begins
    if dm2:
        rho_coeffs[6,3][2] = (-169.0 + 742.0 * nu - 750.0 * nu2 + 156.0 * nu3) / (
            144.0 * (dm2 + 3.0 * nu2)
        )
        rho_coeffs[6,3][3] = -2.0 * a / 9.0
    # (6,3) mode ends

    # (6,2) mode begins
    rho_coeffs[6,2][2] = (-74.0 + 378.0 * nu - 413.0 * nu2 + 49.0 * nu3) / (
        84.0 * (1.0 - 5.0 * nu + 5.0 * nu2)
    )
    rho_coeffs[6,2][3] = (-2.0 * a) / 3.0
    rho_coeffs[6,2][4] = -817991.0 / 3298680.0 + a2 / 2.0
    # (6,2) mode ends

    # (6,1) mode begins
    if dm2:
        rho_coeffs[6,1][2] = (-161.0 + 694.0 * nu - 670.0 * nu2 + 124.0 * nu3) / (
            144.0 * (dm2 + 3.0 * nu2)
        )
        rho_coeffs[6,1][3] = -2.0 * a / 9.0
    # (6,1) mode ends

    # l=7 modes begin
    if dm2:

        rho_coeffs[7,7][2] = (-906.0 + 4246.0 * nu - 4963.0 * nu2 + 1380.0 * nu3) / (
            714.0 * (dm2 + 3.0 * nu2)
        )
        rho_coeffs[7,7][3] = -2.0 * a / 3.0

    rho_coeffs[7,6][2] = (
        2144.0 - 16185.0 * nu + 37828.0 * nu2 - 29351.0 * nu3 + 6104.0 * nu2 * nu2
    ) / (1666.0 * (-1 + 7 * nu - 14 * nu2 + 7 * nu3))

    if dm2:
        rho_coeffs[7,5][2] = (-762.0 + 3382.0 * nu - 3523.0 * nu2 + 804.0 * nu3) / (
            714.0 * (dm2 + 3.0 * nu2)
        )
        rho_coeffs[7,5][3] = -2.0 * a / 3.0

    rho_coeffs[7,4][2] = (
        17756.0
        - 131805.0 * nu
        + 298872.0 * nu2
        - 217959.0 * nu3
        + 41076.0 * nu2 * nu2
    ) / (14994.0 * (-1.0 + 7.0 * nu - 14.0 * nu2 + 7.0 * nu3))

    if dm2:
        rho_coeffs[7,3][2] = (-666.0 + 2806.0 * nu - 2563.0 * nu2 + 420.0 * nu3) / (
            714.0 * (dm2 + 3.0 * nu2)
        )
        rho_coeffs[7,3][3] = -2.0 * a / 3.0

    rho_coeffs[7,2][2] = (
        16832.0
        - 123489.0 * nu
        + 273924.0 * nu2
        - 190239.0 * nu3
        + 32760.0 * nu2 * nu2
    ) / (14994.0 * (-1.0 + 7.0 * nu - 14.0 * nu2 + 7.0 * nu3))

    if dm2:
        rho_coeffs[7,1][2] = (-618.0 + 2518.0 * nu - 2083.0 * nu2 + 228.0 * nu3) / (
            714.0 * (dm2 + 3.0 * nu2)
        )
        rho_coeffs[7,1][3] = -2.0 * a / 3.0

    # l =7 modes end

    # l=8 modes begin
    rho_coeffs[8,8][2] = (
        3482.0 - 26778.0 * nu + 64659.0 * nu2 -
            53445.0 * nu3 + 12243.0 * nu2 * nu2
    ) / (2736.0 * (-1.0 + 7.0 * nu - 14.0 * nu2 + 7.0 * nu3))

    if dm2:
        rho_coeffs[8,7][2] = (
            23478.0
            - 154099.0 * nu
            + 309498.0 * nu2
            - 207550.0 * nu3
            + 38920 * nu2 * nu2
        ) / (18240.0 * (-1 + 6 * nu - 10 * nu2 + 4 * nu3))

    rho_coeffs[8,6][2] = (
        1002.0 - 7498.0 * nu + 17269.0 * nu2 - 13055.0 * nu3 + 2653.0 * nu2 * nu2
    ) / (912.0 * (-1.0 + 7.0 * nu - 14.0 * nu2 + 7.0 * nu3))

    if dm2:
        rho_coeffs[8,5][2] = (
            4350.0
            - 28055.0 * nu
            + 54642.0 * nu2
            - 34598.0 * nu3
            + 6056.0 * nu2 * nu2
        ) / (3648.0 * (-1.0 + 6.0 * nu - 10.0 * nu2 + 4.0 * nu3))

    rho_coeffs[8,4][2] = (
        2666.0 - 19434.0 * nu + 42627.0 * nu2 - 28965.0 * nu3 + 4899.0 * nu2 * nu2
    ) / (2736.0 * (-1.0 + 7.0 * nu - 14.0 * nu2 + 7.0 * nu3))

    if dm2:
        rho_coeffs[8,3][2] = (
            20598.0
            - 131059.0 * nu
            + 249018.0 * nu2
            - 149950.0 * nu3
            + 24520.0 * nu2 * nu2
        ) / (18240.0 * (-1.0 + 6.0 * nu - 10.0 * nu2 + 4.0 * nu3))

    rho_coeffs[8,2][2] = (
        2462.0 - 17598.0 * nu + 37119.0 * nu2 - 22845.0 * nu3 + 3063.0 * nu2 * nu2
    ) / (2736.0 * (-1.0 + 7.0 * nu - 14.0 * nu2 + 7.0 * nu3))

    if dm2:
        rho_coeffs[8,1][2] = (
            20022.0
            - 126451.0 * nu
            + 236922.0 * nu2
            - 138430.0 * nu3
            + 21640.0 * nu2 * nu2
        ) / (18240.0 * (-1.0 + 6.0 * nu - 10.0 * nu2 + 4.0 * nu3))

    if extra_PN_terms:
        # NB: One needs to be careful when adding new terms here.
        # If terms at the given order already exist then it may be
        # appropriate to _add_ new terms to values defined above,
        # i.e with "+="". However, for brand new terms they should
        # just be defined with "=" since the coefficients are not
        # reset to 0 at every step.
        # These terms were not in SEOBNRv4HM
        # Add the NLO  spin-squared term at 3PN
        # We need to get rid of the test-spin term
        rho_coeffs[2,2][6] -=(89.0 * atemp*atemp) / 252.0
        rho_coeffs[2,2][6] += (
            ((178. - 457.*nu- 972.*nu**2)*chiA**2)/504. +
            (dm*(178. - 781.*nu)*chiA*chiS)/252. +
            ((178. - 1817.*nu+ 560.*nu**2)*chiS**2)/504.)
        # Add LO spin-cubed term at 3.5 PN
        # We need to get rid of the test-spin term
        # We need atemp because we have set a to zero afer the (2,2) mode
        rho_coeffs[2,2][7] -= atemp**3/3.0
        rho_coeffs[2,2][7] += (((dm - 4.*dm*nu)*chiA**3)/3. +
                (1. - 3.*nu - 4.*nu**2)*chiA**2*chiS +
                ( dm + 2.*dm*nu)*chiA*chiS**2 +
                (1./3. + nu)*chiS**3)

        # NB: from now on we don't need to worry about subtracting out the
        # test-spin terms since they were set to 0 for modes above (2,2)

        # All the known spin terms in (3,2)
        rho_coeffs[3,2][2] +=  -(16.*nu**2*chiS**2)/(9.*(1. - 3.*nu)**2)
        rho_coeffs[3,2][3] += ((dm*(2. + 13.*nu)*chiA)/(9. - 27.*nu) +
                                ((90. - 1478.*nu + 2515.*nu**2 + 3035.*nu**3)*chiS)/
                                (405.*(1. - 3.*nu)**2) - (320.*nu**3*chiS**3)/
                                (81.*(-1. + 3.*nu)**3))
        rho_coeffs[3,2][4] +=  (((1. - 9.*nu + 12.*nu**2)*chiA**2)/(3. - 9.*nu) -
                                (2.*dm*(-9. + 44.*nu + 25.*nu**2)*chiA*chiS)/
                                (27.*(1. - 3.*nu)**2) + ((-81. + 387.*nu - 1435.*nu**2 +
                                1997.*nu**3 + 2452.*nu**4)*chiS**2)/(243.*(-1 + 3.*nu)**3))
        rho_coeffs[3,2][5] += ((dm*(-245344. + 1128531.*nu - 1514740.*nu**2 +889673.*nu**3)*chiA)/(106920.*(1 - 3*nu)**2) +
                            ((2208096. - 20471053.*nu + 70519165.*nu**2 - 101706029.*nu**3 +40204523.*nu**4 + 11842250.*nu**5)*chiS)/
                            (962280.*(-1. + 3.*nu)**3) - (8.*nu*(1. - 9.*nu + 12.*nu**2)*chiA**2*chiS)/(9.*(1. - 3.*nu)**2) -
                            (16.*dm*nu*(-9. + 46.*nu + 38.*nu**2)*chiA*chiS**2)/(81.*(-1. + 3.*nu)**3) + (8.*nu*(-243. + 1269.*nu - 5029.*nu**2 +
                            5441.*nu**3 + 12022.*nu**4)*chiS**3)/(2187.*(1. - 3.*nu)**4))

        # All the known spin terms in (4,3)
        # NB: odd m so these are f-coeffs and we must treat the dm->0 issue as always

        if dm2:

            f_coeffs[4,3][3] = ((nu*(-2661. + 3143.*nu)*chiA)/(132.*dm*(-1. + 2.*nu)) +
                            (23.*nu*(87. + 23.*nu)*chiS)/(132.*(-1. + 2.*nu)))


            f_coeffs[4,3][4] = (((9. - 74.*nu + 72.*nu**2)*chiA**2)/(6. - 12.*nu) +
                            ((18. - 108.*nu + 137.*nu**2)*chiA*chiS)/(6.*dm - 12.*dm*nu) + ((9. + 2.*nu + 35.*nu**2)*
                            chiS**2)/(6. - 12.*nu))


        else:
            f_coeffs[4,3][3] = ((nu*(-2661. + 3143.*nu)*chiA)/(132.*(-1. + 2.*nu)))
            f_coeffs[4,3][4] = ((18. - 108.*nu + 137.*nu**2)*chiA*chiS)/(6. - 12.*nu)

    return rho_coeffs, rho_coeffs_log, f_coeffs, f_coeffs_vh

def compute_delta_coeffs( nu, dm, a, chiS, chiA):

    """
    Compute the phase residual coefficients.
    See Sec. 2B of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.79.064004
    See Eq. 59-63 of https://dcc.ligo.org/LIGO-T2300060 (SEOBNRv5_theory.pdf) for new terms, rest copied from SEOBNRv4HM LAL code

    Coefficients can be found in:
    - https://dcc.ligo.org/LIGO-T2300060 (SEOBNRv5_theory.pdf and SEOBNRv5HM.pdf)
    - Appendix A of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.98.084028
    - Eqs. (2.4) to (2.6) of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.95.044028
    - https://journals.aps.org/prd/pdf/10.1103/PhysRevD.83.064003
    - https://journals.aps.org/prd/pdf/10.1103/PhysRevD.86.024011
    """
    
    delta_coeffs = zeros([ell_max+1,ell_max+1,PN_limit])
    delta_coeffs_vh = zeros([ell_max+1,ell_max+1,PN_limit])

    nu2 = nu*nu
    nu3 = nu*nu2

    a2 = a*a
    a3 = a*a2
    dm2 = dm*dm

    chiA2 = chiA * chiA
    chiS2 = chiS * chiS
    chiA3 = chiA2 * chiA
    chiS3 = chiS2 * chiS

    m1Plus3nu = -1.0 + 3.0 * nu
    m1Plus3nu2 = m1Plus3nu * m1Plus3nu
    m1Plus3nu3 = m1Plus3nu * m1Plus3nu2

    aDelta = 0.0

    #(2,2) mode begins
    delta_coeffs_vh[2,2][3] = 7.0 / 3.0
    # See https://dcc.ligo.org/T1600383

    delta_coeffs_vh[2,2][6] = (
        -4.0 / 3.0 * (dm * chiA + chiS * (1 - 2 * nu)) + (428.0 * pi) / 105.0
    )

    delta_coeffs[2,2][8] = (20.0 * aDelta) / 63.0
    delta_coeffs_vh[2,2][9] = -2203.0 / 81.0 + (1712.0 * pi * pi) / 315.0
    delta_coeffs[2,2][5] = -24.0 * nu
    delta_coeffs[2,2][6] = 0.0

    # (2,2) mode ends

    # (2,1) mode begins

    delta_coeffs_vh[2,1][3] = 2.0 / 3.0
    delta_coeffs_vh[2,1][6] = (-17.0 * aDelta) / 35.0 + (107.0 * pi) / 105.0
    delta_coeffs_vh[2,1][7] = (3.0 * aDelta * aDelta) / 140.0
    delta_coeffs_vh[2,1][9] = -272.0 / 81.0 + (214.0 * pi * pi) / 315.0
    """
    This was in SEOBNRv4HM
    delta_coeffs[2,1][5] = -493.0 * nu / 42.0
    """
    delta_coeffs[2,1][5] = -25.0 * nu /2

    # (2,1) mode ends

    # (3,3) mode begins
    # l = 3, Eqs. A9a - A9c for rho, Eqs. A15b and A15c for f,
    # Eqs. 22 - 24 of DIN and Eqs. 27c - 27e of PBFRT for delta
    delta_coeffs_vh[3,3][3] = 13.0 / 10.0
    delta_coeffs_vh[3,3][6] = (-81.0 * aDelta) / 20.0 + (39.0 * pi) / 7.0
    delta_coeffs_vh[3,3][9] = -227827.0 / 3000.0 + (78.0 * pi * pi) / 7.0
    delta_coeffs[3,3][5] = -80897.0 * nu / 2430.0
    # (3,3) mode ends

    # (3,2) mode begins
    delta_coeffs_vh[3,2][3] = (10.0 + 33.0 * nu) / (-15.0 * m1Plus3nu)
    delta_coeffs_vh[3,2][4] = 4.0 * aDelta
    delta_coeffs_vh[3,2][6] = (-136.0 * aDelta) / 45.0 + (52.0 * pi) / 21.0
    delta_coeffs_vh[3,2][9] = -9112.0 / 405.0 + (208.0 * pi * pi) / 63.0
    # (3,2) mode ends

    # (3,1) mode begins
    delta_coeffs_vh[3,1][3] = 13.0 / 30.0
    delta_coeffs_vh[3,1][6] = (61.0 * aDelta) / 20.0 + (13.0 * pi) / 21.0
    delta_coeffs_vh[3,1][7] = (-24.0 * aDelta * aDelta) / 5.0
    delta_coeffs_vh[3,1][9] = -227827.0 / 81000.0 + (26.0 * pi * pi) / 63.0
    delta_coeffs[3,1][5] = -17.0 * nu / 10.0
    # (3,1) mode ends

    # (4,4) mode begins
    delta_coeffs_vh[4,4][3] = (112.0 + 219.0 * nu) / (-120.0 * m1Plus3nu)
    delta_coeffs_vh[4,4][6] = (-464.0 * aDelta) / 75.0 + (25136.0 * pi) / 3465.0

    # RC: This terms are in Eq.A15 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
    delta_coeffs_vh[4,4][9] = -55144.0 / 375.0 + 201088.0 * pi * pi / 10395.0

    # SO: what is going on with delta44v5: it's declared but never set in the LAL code!

    # (4,4) mode ends

    # (4,3) mode begins
    delta_coeffs_vh[4,3][3] = (486.0 + 4961.0 * nu) / (810.0 * (1.0 - 2.0 * nu))
    delta_coeffs_vh[4,3][4] = (11.0 * aDelta) / 4.0
    delta_coeffs_vh[4,3][6] = 1571.0 * pi / 385.0
    # (4,3) mode ends

    # (4,2) mode begins
    delta_coeffs_vh[4,2][3] = (7.0 * (1.0 + 6.0 * nu)) / (-15.0 * m1Plus3nu)
    delta_coeffs_vh[4,2][6] = (212.0 * aDelta) / 75.0 + (6284.0 * pi) / 3465.0
    # (4,2) mode ends

    # (4,1) mode begins
    delta_coeffs_vh[4,1][3] = (2.0 + 507.0 * nu) / (10.0 * (1.0 - 2.0 * nu))
    delta_coeffs_vh[4,1][4] = (11.0 * aDelta) / 12.0
    delta_coeffs_vh[4,1][6] = 1571.0 * pi / 3465.0
    # (4,1) mode ends

    # l=5 modes begin
    delta_coeffs_vh[5,5][3] = (96875.0 + 857528.0 * nu) / (131250.0 * (1 - 2 * nu))

    # RC: This terms are in Eq.A16 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
    delta_coeffs_vh[5,5][6] = 3865.0 / 429.0 * pi
    delta_coeffs_vh[5,5][9] = (
        -7686949127.0 + 954500400.0 * pi * pi
    ) / 31783752.0

    delta_coeffs_vh[5,4][3] = 8.0 / 15.0
    delta_coeffs_vh[5,4][4] = 12.0 * aDelta / 5.0

    delta_coeffs_vh[5,3][3] = 31.0 / 70.0

    delta_coeffs_vh[5,2][3] = 4.0 / 15.0
    delta_coeffs_vh[5,2][4] = 6.0 * aDelta / 5.0

    delta_coeffs_vh[5,1][3] = 31.0 / 210.0
    # ell = 5 modes end

    # l=6 modes begin
    delta_coeffs_vh[6,6][3] = 43.0 / 70.0
    delta_coeffs_vh[6,5][3] = 10.0 / 21.0
    delta_coeffs_vh[6,4][3] = 43.0 / 105.0
    delta_coeffs_vh[6,3][3] = 2.0 / 7.0
    delta_coeffs_vh[6,2][3] = 43.0 / 210.0
    delta_coeffs_vh[6,1][3] = 2.0 / 21.0
    # l=6 modes end

    # l=7 modes begin
    delta_coeffs_vh[7,7][3] = 19.0 / 36.0
    delta_coeffs_vh[7,5][3] = 95.0 / 252.0
    delta_coeffs_vh[7,3][3] = 19.0 / 84.0
    delta_coeffs_vh[7,1][3] = 19.0 / 252.0
    # l=7 modes end

    return delta_coeffs, delta_coeffs_vh

def compute_deltalm_single(vs, vhs, l, m,delta_coeffs,delta_coeffs_vh):
    """
    Compute the  full \delta_{\ell m} contribution for a given mode
    """
    delta = 0
    for j in range(PN_limit):
        delta += delta_coeffs[l,m,j]*vs[j] + delta_coeffs_vh[l,m,j]*vhs[j]
    return delta

def compute_delta(v, vh, nu, delta_coeffs,delta_coeffs_vh):
    """
    Compute the  full \delta_{\ell m} contribution for all modes
    """
    vs = zeros(PN_limit)
    vhs = zeros(PN_limit)
    deltalm = zeros([ell_max+1,ell_max+1])
    for i in range(PN_limit):
        vs[i] = v**i
        vhs[i] = vh**i
    for l in range(2,ell_max+1):
        for m in range(1,l+1):
            deltalm[l,m] = compute_deltalm_single(vs,vhs,l,m,delta_coeffs,delta_coeffs_vh)
    return deltalm

def compute_extra_flm_terms(l, m, vh,f_coeffs_vh):
    """
    Compute the complex term in f_{33}. See last term in Eq(A10) of https://arxiv.org/pdf/1803.10701.pdf
    """
    vh3 = vh**3
    extra_term = 0.0
    if l==3 and m==3:
        extra_term = I*vh3 * vh3 * f_coeffs_vh[3,3][6]
    return extra_term

def compute_rholm_single(nu, vs, vh, l, m, rho_coeffs, rho_coeffs_log, extra_coeffs, extra_coeffs_log, f_coeffs, f_coeffs_vh):
    """
    Compute the full \rho_{\ell m}​ contribution for a given mode
    """
    v = vs[1]
    rho_final  = 0.0
    rho = 1.0

    rho_log = 0.0
    f = 0.0
    f_final = 0.0
    eulerlogxabs =  euler_gamma + log(2.0 * m * v)
    # This is just computing rho_coeffs \dot v

    if m%2:
        # For odd m modes we need to compute f
        for j in range(1,PN_limit):
            rho+= (rho_coeffs[l,m,j]+extra_coeffs[l,m,j]+ (rho_coeffs_log[l,m,j]+extra_coeffs_log[l,m,j])*eulerlogxabs)*vs[j]
            f += f_coeffs[l,m,j]*vs[j]
    else:
        # For even m modes we only need rho
        for j in range(1,PN_limit):
            rho+= (rho_coeffs[l,m,j]+extra_coeffs[l,m,j]+(rho_coeffs_log[l,m,j]+extra_coeffs_log[l,m,j])*eulerlogxabs)*vs[j]

    # Deal with the complex amplitude term
    f_final = f
    #if l==3 and m==3:
    #    extra_flm_term = compute_extra_flm_terms(l,m,vh,f_coeffs_vh)
    #    f_final += extra_flm_term


    rho_final = rho**l
    if absolute(nu-0.25)<1e-14 and m%2:
        rho_final = f_final
    else:
        rho_final += f_final
    return rho_final

def compute_rholm(v, vh, nu,  rho_coeffs, rho_coeffs_log, extra_coeffs, extra_coeffs_log, f_coeffs, f_coeffs_vh):
    vs = zeros(PN_limit)
    vs[0] = 0
    vs[1] = v
    for i in range(2,PN_limit):
        vs[i] = v*vs[i-1]
    rholm = zeros([ell_max+1,ell_max+1],dtype = complex)
    for l in range(2,ell_max+1):
        for m in range(1,l+1):
            rholm[l,m] = compute_rholm_single(nu, vs,vh,l,m, rho_coeffs, rho_coeffs_log, extra_coeffs, extra_coeffs_log, f_coeffs, f_coeffs_vh)
    return rholm

def  EOBFluxCalculateNewtonianMultipoleAbs(x, phi, l, m, params):
    """
    Compute the Newtonian multipole (optimised for the flux calculation).
    """
    param = params[l,m]
    epsilon = (l + m) % 2


    y = ylms[m][l-epsilon]
    multipole = param * x ** ((l + epsilon) / 2.0)
    multipole *= y
    return multipole

def update_rho_coeffs(rho_coeffs, extra_coeffs):
    rho_coeffs_new = zeros([ell_max+1,ell_max+1,PN_limit])
    temp = 0.0
    for l in range(2,ell_max+1):
        for m in range(1,l+1):
            for i in range(PN_limit):
                temp = extra_coeffs[l,m,i]
                if absolute(temp)>1e-15:
                    rho_coeffs_new[l,m,i] = rho_coeffs[l,m,i] + temp
    return rho_coeffs_new

def compute_flux(m1, m2, r, phi, pr, pphi, omega, omega_circ, H,chi1,chi2):
    """
    Compute the full flux. See Eq(43) in the .
    """
    # Note the "abs" in the prefixes.
    prefixes = compute_newtonian_prefixes(m1,m2)
    prefixes_abs = array([absolute(prefixes[i]) for i in range(len(prefixes))])
    v = omega**(1./3)
    vh3 = H*omega
    vh = vh3**(1./3)
    omega2 = omega*omega
    nu = m1*m2/((m1+m2)**2)
    flux = 0.0
    # Precompute the tail
    Tlm = compute_tail(omega, H)

    # Precompute the source term
    Slm = 0.0
    source1 = (H * H - 1.0) / (2.0 * nu) + 1.0 # H_eff
    source2 = v*pphi
    v_phi = omega/omega_circ**(2./3)
    v_phi2 = v_phi*v_phi

    rholmPwrl  = 1.0
    chiA = 0.5*(chi1 - chi2)
    chiS = 0.5*(chi1 + chi2)
    dm = (m1 - m2)/(m1 + m2)
    a = .5
    # Assume that the spin params have already been updated appropriately
    rho_coeffs, rho_coeffs_log, f_coeffs, f_coeffs_vh = compute_rho_coeffs(nu,dm,a,chiS,chiA,1)
    
    extra_coeffs = zeros([ell_max+1,ell_max+1,PN_limit])
    extra_coeffs_log = zeros([ell_max+1,ell_max+1,PN_limit])
    coeffs_arrays = [
        21.2,
        -411.0,
        12.0,
        -215.0,
        1.65,
        26.5,
        80.0,
        -3.56,
        15.6,
        -216.0,
        -2.61,
        1.25,
        -35.7,
        0.333,
        -6.5,
        (98 - (1312549797426453052 / 176264081083715625) / nu),
        (18778864 / 12629925) / nu,
        -0.654,
        -3.69,
        18.5 - (2465107182496333 / 460490801971200) / nu,
        (174381 / 67760) / nu,
    ]

    (
        h22_v8,
        h22_v10,
        h33_v8,
        h33_v10,
        h21_v6,
        h21_v8,
        h21_v10,
        h44_v6,
        h44_v8,
        h44_v10,
        h55_v4,
        h55_v6,
        h55_v8,
        h32_v6,
        h32_v8,
        h32_v10,
        h32_vlog10,
        h43_v4,
        h43_v6,
        h43_v8,
        h43_vlog8,
    ) = nu * array(coeffs_arrays)
    extra_coeffs[2,2,8] = h22_v8
    extra_coeffs[2,2,10] = h22_v10
    extra_coeffs[3,3,8] = h33_v8
    extra_coeffs[3,3,10] = h33_v10
    extra_coeffs[2,1,6] = h21_v6
    extra_coeffs[2,1,8] = h21_v8
    extra_coeffs[2,1,10] = h21_v10
    extra_coeffs[4,4,6] = h44_v6
    extra_coeffs[4,4,8] = h44_v8
    extra_coeffs[4,4,10] = h44_v10
    extra_coeffs[5,5,4] = h55_v4
    extra_coeffs[5,5,6] = h55_v6
    extra_coeffs[5,5,8] = h55_v8
    extra_coeffs[3,2,6] = h32_v6
    extra_coeffs[3,2,8] = h32_v8
    extra_coeffs[3,2,10] = h32_v10
    extra_coeffs_log[3,2,10] = h32_vlog10
    extra_coeffs[4,3,4] = h43_v4
    extra_coeffs[4,3,6] = h43_v6
    extra_coeffs[4,3,8] = h43_v8
    extra_coeffs_log[4,3,8] = h43_vlog8

    rholm = compute_rholm(v, vh, nu,  rho_coeffs, rho_coeffs_log, extra_coeffs, extra_coeffs_log, f_coeffs, f_coeffs_vh)
    
    for l in range(2,ell_max+1):
        for m in range(1,l+1):
            # Assemble the waveform

            hNewton = EOBFluxCalculateNewtonianMultipoleAbs(
                v_phi2, pi/2, l, m, prefixes_abs
            )
            if ((l + m) % 2) == 0:
                Slm = source1

            else:
                Slm = source2

            tail = Tlm[l,m]
            rholmPwrl = rholm[l,m]
            hlm = tail*Slm*hNewton*rholmPwrl

            flux += m*m*omega2*absolute(hlm)**2
    return rholm#-flux/(8*pi)

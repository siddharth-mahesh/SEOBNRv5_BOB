from numpy import log,sqrt,array
def evaluate_H(q,p, chi_1, chi_2, m_1, m_2,verbose = False):
        """
        Evaluate the Hamiltonian and xi

        Args:
          q ([:]): Canonical positions (r,phi).
          p ([:]): Canonical momenta  (prstar,pphi).
          chi1 (): Dimensionless z-spin of the primary.
          chi2 (): Dimensionless z-spin of the secondary.
          m_1 (): Primary mass component.
          m_2 (): Secondary mass component.
          
        Returns:
           (tuple)  H,xi

        """
        # Coordinate definitions
        M = m_1 + m_2
        X_1 = m_1/M
        X_2 = m_2/M
        nu = X_1*X_2
        
        r = q[0]
        phi = q[1]
        u = 1/r
        prst = p[0]
        L = p[1]

        pphi = L


        r2 = r*r
        r3 = r2*r
        r4 = r2*r2
        r5 = r*r4

        L2 = L*L
        L4 = L2*L2
        lr = log(r)

        nu2 = nu*nu
        nu3 = nu2*nu
        nu4 = nu3*nu

        prst2 = prst*prst
        prst4 = prst2*prst2
        prst6 = prst4*prst2
        prst8 = prst6*prst2

        # Actual Hamiltonian expressions
        d5 = 0
        par = array([4.17877875e01, -3.02193382e03, 3.34144394e04, -1.69019140e05, 3.29523262e05])
        a6 = par[0] + par[1] * nu + par[2] * nu ** 2 + par[3] * nu ** 3 + par[4] * nu ** 4
        Dbpm = r*(6730497718123.02*nu3 + 22295347200.0*nu2*d5 + 133772083200.0*nu2*r2 + 1822680546449.21*nu2*r + 80059249540278.2*nu2 + 22295347200.0*nu*d5*r - 193226342400.0*nu*d5 + 2589101062873.81*nu*r2 + 10611661054566.2*nu*r - 12049908701745.2*nu + 5107745331375.71*r2 - 326837426.241486*r*(14700.0*nu + 42911.0) - 39476764256925.6*r - (-5041721180160.0*nu2 - 25392914995744.3*nu - 879923036160.0*r2 - 283115520.0*r*(14700.0*nu + 42911.0) + 104186110149937.0)*lr + 5787938193408.0*lr**2 + 275059053208689.0)/(55296.0*nu*(14515200.0*nu3 - 42636451.6032331*nu2 - 7680.0*nu*(315.0*d5 + 890888.810272497) + 4331361844.61149*nu + 1002013764.01019) - 967680.0*r3*(-138240.0*nu2 - 2675575.66847905*nu - 5278341.3229329) - 9216.0*r2*(-197773496.793534*nu2 - 7680.0*nu*(315.0*d5 + 405152.309729121) + 2481453539.84635*nu + 5805304367.87913) + r*(5927865218923.02*nu3 + 70778880.0*nu2*(315.0*d5 + 2561145.80918574) - 138141470005001.0*nu2 - 4718592.0*nu*(40950.0*d5 + 86207832.4415642) + 450172889755120.0*nu + 86618264430493.3*(1 - 0.496948781616935*nu)**2 + 188440788778196.0) + 5787938193408.0*r*lr**2 + (-1698693120.0*nu*(11592.0*nu + 69847.0) + 879923036160.0*r3 + 283115520.0*r2*(14700.0*nu + 42911.0) + 49152.0*r*(102574080.0*nu2 + 409207698.136075*nu - 2119671837.36038))*lr)

        Apm = 7680.0*r4*(-5416406.59541186*nu2 + 28.0*nu*(1920.0*a6 + 733955.307463037) + 2048.0*nu*(756.0*nu + 336.0*r + 407.0)*lr - 7.0*r*(-185763.092693281*nu2 + 938918.400156317*nu - 245760.0) - 3440640.0)/(241555486248.807*nu4 + 1120.0*nu3*(-17833256.898555*r2 - 163683964.822551*r - 1188987459.03162) + 7.0*nu2*(-39321600.0*a6*(3.0*r + 59.0) + 745857848.115604*a6 + 1426660551.8844*r5 - 3089250703.76879*r4 - 6178501407.53758*r3 + 2064783811.32587*r2 + 122635399361.987*r + 276057889687.011) + 67645734912.0*nu2*lr**2 + 53760.0*nu*(7680.0*a6*(r4 + 2.0*r3 + 4.0*r2 + 8.0*r + 16.0) + 128.0*r*(-6852.34813868015*r4 + 4264.6962773603*r3 + 8529.39255472061*r2 + 13218.7851094412*r - 33722.4297811176) + 113485.217444961*r*(-r4 + 2.0*r3 + 4.0*r2 + 8.0*r + 16.0) + 148.04406601634*r*(349.0*r4 + 1926.0*r3 + 3852.0*r2 + 7704.0*r + 36400.0)) + 32768.0*nu*(-1882456.23663972*nu2 - 38842241.4769507*nu + 161280.0*r5 + 480.0*r4*(756.0*nu + 1079.0) + 960.0*r3*(756.0*nu + 1079.0) + 1920.0*r2*(588.0*nu + 1079.0) + 240.0*r*(-3024.0*nu2 - 7466.27061066206*nu + 17264.0) + 13447680.0)*lr + 13212057600.0*r5)

        ap = chi_1*X_1 + chi_2*X_2

        ap2 = ap*ap

        xi = sqrt(Dbpm)*r2*(Apm + ap2/r2)/(ap2 + r2)

        pr = prst/xi

        flagNLOSS2 = 1.00000000000000

        delta = X_1 - X_2

        am = chi_1*X_1 - chi_2*X_2

        apam = am*ap

        am2 = am*am

        apamd = apam*delta

        dSO = (-7.71251231383957 * am ** 3- 17.2294679794015 * am ** 2 * ap- 238.430383378296 * am ** 2 * nu+ 69.5461667822545 * am ** 2- 10.5225438990315 * am * ap ** 2+ 362.767393298729 * am * ap * nu- 85.8036338010274 * am * ap- 1254.66845939312 * am * nu ** 2+ 472.431937787377 * am * nu- 39.742317057316 * am- 7.58458103577458 * ap ** 3- 42.7601129678844 * ap ** 2 * nu+ 18.1783435552183 * ap ** 2- 201.905934468847 * ap * nu ** 2- 90.5790079104259 * ap * nu+ 49.6299175121658 * ap+ 478.546231305475 * nu ** 3+ 679.521769948995 * nu ** 2- 177.334831768076 * nu- 37.6897780220529)

        QSalign2 = flagNLOSS2*pr**4*(-0.46875*am2*(4.0*nu2 - 5.0*nu + 1.0) - 0.15625*ap2*(32.0*nu2 - 33.0*nu - 5.0) + 0.3125*apamd*(18.0*nu - 1.0))/r3

        flagQPN55 = 1.00000000000000

        flagQPN5 = 1.00000000000000

        flagQPN4 = 1.00000000000000

        Qpm = flagQPN4*(0.121954868780449*nu*prst8/r + prst6*(6.0*nu3 - 5.4*nu2 - 2.78300763695006*nu)/r2 + prst4*(10.0*nu3 - 131.0*nu2 + 92.7110442849544*nu)/r3) + flagQPN5*(prst8*(-6.0*nu4 + 3.42857142857143*nu3 + 3.33842023648322*nu2 + 1.38977750996128*nu)/r2 + prst6*(-14.0*nu4 + 188.0*nu3 - 89.5298327361234*nu2 - 33.9782122170436*nu)/r3 + prst4*(602.318540416564*nu3 + nu2*(118.4*lr - 1796.13660498019) + nu*(452.542166996693 - 51.6952380952381*lr))/r4) + flagQPN55*(1.48275342024365*nu*prst8/r**2.5 - 11.3175085791863*nu*prst6/r**3.5 + 147.443752990146*nu*prst4/r**4.5) + prst4*(-6.0*nu2 + 8.0*nu)/r2

        Qq = QSalign2 + Qpm

        Bnpa = -r*(r + 2.0)/(ap2*r*(r + 2.0) + r4)

        flagNLOSS = 1.00000000000000

        BnpSalign2 = flagNLOSS*(0.1875*am2*(4.0*nu - 1.0) + ap2*(3.0*nu + 2.8125) - 2.625*apamd)/r3 + flagNLOSS2*(0.015625*am2*(4.0*nu2 + 115.0*nu - 37.0) + 0.015625*ap2*(-1171.0*nu - 861.0) + 0.03125*apamd*(26.0*nu + 449.0))/r4

        Bnp = Apm*Dbpm + BnpSalign2 + ap2/r2 - 1.0

        ASalignCal2 = 0.0

        ASalign2 = flagNLOSS*(0.125*am2*(4.0*nu + 1.0) + 1.125*ap2 - 1.25*apamd)/r4 + flagNLOSS2*(0.046875*am2*(28.0*nu2 - 27.0*nu - 3.0) - 0.390625*ap2*(7.0*nu + 9.0) - 1.21875*apamd*(2.0*nu - 3.0))/r**5

        A = (ASalign2 + ASalignCal2 + Apm + ap2/r2)/(ap2*(1.0 + 2.0/r)/r2 + 1.0)

        lap = ap

        Heven = sqrt(A*(Bnpa*L2*lap**2/r2 + L2/r2 + Qq + prst2*(Bnp + 1.0)/xi**2 + 1.0))

        lam = am

        Ga3 = 0.0416666666666667*L*ap2*delta*lam/r2 + L*lap*(-0.25*ap2 + 0.208333333333333*apamd)/r2

        SOcalib = L*nu*dSO*lap/r3

        flagNLOSO2 = 1.00000000000000

        flagNLOSO = 1.00000000000000

        gam = flagNLOSO*(L2*(0.46875 - 0.28125*nu)/r2 + (0.34375*nu + 0.09375)/r) + flagNLOSO2*(L4*(0.29296875*nu2 - 0.3515625*nu - 0.41015625)/r4 + L2*(-0.798177083333333*nu2 - 0.2734375*nu - 0.23046875)/r3 + (0.536458333333333*nu2 - 0.03125*nu + 0.078125)/r2) + 0.25

        gap = flagNLOSO*(L2*(-1.40625*nu - 0.46875)/r2 + (0.71875*nu - 0.09375)/r) + flagNLOSO2*(L4*(1.34765625*nu2 + 0.5859375*nu + 0.41015625)/r4 + L2*(-2.07161458333333*nu2 - 2.0859375*nu + 0.23046875)/r3 + (0.567708333333333*nu2 - 5.53125*nu - 0.078125)/r2) + 1.75

        Hodd = (Ga3 + L*delta*gam*lam + L*gap*lap + SOcalib)/(2.0*ap2 + 2.0*r2 + r*(ap2 + r2 - 2.0*r))

        Heff = Heven + Hodd

        # Evaluate H_real/nu
        H = M * sqrt(1+2*nu*(Heff-1)) / nu
        if not verbose:
            return H*nu,xi
        else:
            return H*nu,xi,A,Bnp,Bnpa,Qq,Heven,Hodd,QSalign2,Qpm,Ga3,gam,gap,SOcalib,u,nu,ap,am,r,phi,prst,L,chi_1,chi_2,m_1,m_2



def omega(q, p, chi_1, chi_2, m_1, m_2):

        """
        Compute the orbital frequency from the Hamiltonian.

        Args:
          q (double[:]): Canonical positions (r,phi).
          p (double[:]): Canonical momenta  (prstar,pphi).
          chi1 (double): Dimensionless z-spin of the primary.
          chi2 (double): Dimensionless z-spin of the secondary.
          m_1 (double): Primary mass component.
          m_2 (double): Secondary mass component.

        Returns:
           (double) dHdpphi
        """

        # Coordinate definitions

        r = q[0]
        phi = q[1]

        prst = p[0]
        L = p[1]

        pphi = L
        M = m_1 + m_2
        X_1 = m_1/M
        X_2 = m_2/M
        nu = X_1*X_2
        ap = chi_1*X_1 + chi_2*X_2

        am = chi_1*X_1 - chi_2*X_2

        d5 = 0
        par = array([4.17877875e01, -3.02193382e03, 3.34144394e04, -1.69019140e05, 3.29523262e05])
        a6 = par[0] + par[1] * nu + par[2] * nu ** 2 + par[3] * nu ** 3 + par[4] * nu ** 4
        dSO = (-7.71251231383957 * am ** 3- 17.2294679794015 * am ** 2 * ap- 238.430383378296 * am ** 2 * nu+ 69.5461667822545 * am ** 2- 10.5225438990315 * am * ap ** 2+ 362.767393298729 * am * ap * nu- 85.8036338010274 * am * ap- 1254.66845939312 * am * nu ** 2+ 472.431937787377 * am * nu- 39.742317057316 * am- 7.58458103577458 * ap ** 3- 42.7601129678844 * ap ** 2 * nu+ 18.1783435552183 * ap ** 2- 201.905934468847 * ap * nu ** 2- 90.5790079104259 * ap * nu+ 49.6299175121658 * ap+ 478.546231305475 * nu ** 3+ 679.521769948995 * nu ** 2- 177.334831768076 * nu- 37.6897780220529)

        # Extra quantities used in the Jacobian
        z0 = r**2
        z1 = chi_1*X_1
        z2 = chi_2*X_2
        z3 = z1 + z2
        z4 = z3**2
        z5 = 2.0*z4
        z6 = z0 + z4
        z7 = r**3
        z8 = 1/z7
        z9 = 1/z0
        z10 = z4*z9
        z11 = z1 - z2
        z12 = z11*(X_1 - X_2)
        z13 = 5.0*X_1 - 5.0*X_2
        z14 = z11*z3
        z15 = -1.40625*nu - 0.46875
        z16 = 2.0*pphi
        z17 = z16*z9
        z18 = nu**2
        z19 = -2.0859375*nu - 2.07161458333333*z18 + 0.23046875
        z20 = z16*z8
        z21 = 0.5859375*nu + 1.34765625*z18 + 0.41015625
        z22 = r**4
        z23 = 1/z22
        z24 = 4.0*pphi**3*z23
        z25 = 0.46875 - 0.28125*nu
        z26 = -0.2734375*nu - 0.798177083333333*z18 - 0.23046875
        z27 = -0.3515625*nu + 0.29296875*z18 - 0.41015625
        z28 = r**(-1)
        z29 = pphi**2
        z30 = z29*z9
        z31 = z29*z8
        z32 = pphi**4*z23
        z33 = r + 2.0
        z34 = z33*z4
        z35 = z28/(r*z34 + z22)
        z36 = z11**2
        z37 = r**5
        z38 = 1/z37
        z39 = 0.015625*z4
        z40 = nu**4
        z41 = log(r)
        z42 = z41**2
        z43 = nu**3
        z44 = 756.0*nu
        z45 = z44 + 1079.0
        z46 = 8.0*r + 4.0*z0 + 2.0*z7 + 16.0
        z47 = (2048.0*nu*z41*(336.0*r + z44 + 407.0) + 28.0*nu*(1920.0*a6 + 733955.307463037) - 7.0*r*(938918.400156317*nu - 185763.092693281*z18 - 245760.0) - 5416406.59541186*z18 - 3440640.0)/(32768.0*nu*z41*(-38842241.4769507*nu + 240.0*r*(-7466.27061066206*nu - 3024.0*z18 + 17264.0) + 1920.0*z0*(588.0*nu + 1079.0) - 1882456.23663972*z18 + 480.0*z22*z45 + 161280.0*z37 + 960.0*z45*z7 + 13447680.0) + 53760.0*nu*(7680.0*a6*(z22 + z46) + 113485.217444961*r*(-z22 + z46) + 148.04406601634*r*(7704.0*r + 3852.0*z0 + 349.0*z22 + 1926.0*z7 + 36400.0) + 128.0*r*(13218.7851094412*r + 8529.39255472061*z0 - 6852.34813868015*z22 + 4264.6962773603*z7 - 33722.4297811176)) + 67645734912.0*z18*z42 + 7.0*z18*(-39321600.0*a6*(3.0*r + 59.0) + 745857848.115604*a6 + 122635399361.987*r + 2064783811.32587*z0 - 3089250703.76879*z22 + 1426660551.8844*z37 - 6178501407.53758*z7 + 276057889687.011) + 13212057600.0*z37 + 241555486248.807*z40 + 1120.0*z43*(-163683964.822551*r - 17833256.898555*z0 - 1188987459.03162))
        z48 = z22*z47
        z49 = z10 + z23*(-0.25*z13*z14 + 0.125*z36*(4.0*nu + 1.0) + 1.125*z4) + z38*(0.03125*z14*(2.0*nu - 3.0)*(-39.0*X_1 + 39.0*X_2) + 0.046875*z36*(-27.0*nu + 28.0*z18 - 3.0) - z39*(175.0*nu + 225.0)) + 7680.0*z48
        z50 = prst**4
        z51 = prst**6
        z52 = prst**8
        z53 = nu*z52
        z54 = 4.0*z18
        z55 = nu*r
        z56 = nu*z0
        z57 = r*z18
        z58 = 14700.0*nu + 42911.0
        z59 = r*z58
        z60 = z0*z18
        z61 = 283115520.0*z58
        z62 = z41*(-25392914995744.3*nu - r*z61 - 879923036160.0*z0 - 5041721180160.0*z18 + 104186110149937.0)
        z63 = r*z42
        z64 = z0*(-630116198.873299*nu - 197773496.793534*z18 + 5805304367.87913)
        z65 = z7*(-2675575.66847905*nu - 138240.0*z18 - 5278341.3229329)
        z66 = nu*(-2510664218.28128*nu - 42636451.6032331*z18 + 14515200.0*z43 + 1002013764.01019)
        z67 = r*(43393301259014.8*nu + 43133561885859.3*z18 + 5927865218923.02*z43 + 86618264430493.3*(1 - 0.496948781616935*nu)**2 + 188440788778196.0)
        z68 = z41*(-1698693120.0*nu*(11592.0*nu + 69847.0) + 49152.0*r*(409207698.136075*nu + 102574080.0*z18 - 2119671837.36038) + z0*z61 + 879923036160.0*z7)
        z69 = 0.000130208333333333*z10 + z48
        z70 = -12049908701745.2*nu - 39476764256925.6*r + 5107745331375.71*z0 + 80059249540278.2*z18 + 5787938193408.0*z42 + 6730497718123.02*z43 + 10611661054566.2*z55 + 2589101062873.81*z56 + 1822680546449.21*z57 - 326837426.241486*z59 + 133772083200.0*z60 - z62 + 275059053208689.0
        z71 = 5787938193408.0*z63 - 9216.0*z64 - 967680.0*z65 + 55296.0*z66 + z67 + z68


        # Evaluate Hamiltonian
        H,xi =  evaluate_H(q,p, chi_1, chi_2, m_1, m_2)

        # Heff Jacobian expressions

        dHeffdpphi = z49*(z49*(147.443752990146*nu*r**(-4.5)*z50 - 11.3175085791863*nu*r**(-3.5)*z51 + 1.69542100694444e-8*prst**2*z38*z6**2*z71*(z10 + z23*(0.03125*z12*z3*(26.0*nu + 449.0) + 0.015625*z36*(115.0*nu + z54 - 37.0) + z39*(-1171.0*nu - 861.0)) + 7680.0*z37*z47*z70/z71 + z8*(0.125*z14*(-21.0*X_1 + 21.0*X_2) + 0.0625*z36*(12.0*nu - 3.0) + z4*(3.0*nu + 2.8125)))/(z69**2*z70) + 1.48275342024365*r**(-2.5)*z53 + z23*z50*(nu*(452.542166996693 - 51.6952380952381*z41) + z18*(118.4*z41 - 1796.13660498019) + 602.318540416564*z43) + 0.121954868780449*z28*z53 - z29*z34*z35 + z30 + z50*z8*(92.7110442849544*nu - 131.0*z18 + 10.0*z43) + z50*z9*(8.0*nu - 6.0*z18) + z51*z8*(-33.9782122170436*nu - 89.5298327361234*z18 - 14.0*z40 + 188.0*z43) + z51*z9*(-2.78300763695006*nu - 5.4*z18 + 6.0*z43) + z52*z9*(1.38977750996128*nu + 3.33842023648322*z18 - 6.0*z40 + 3.42857142857143*z43) + 1.0 + 1.27277314139085e-19*z50*z6**4*(0.0625*z11*z13*z3*(18.0*nu - 1.0) - 0.46875*z36*(-5.0*nu + z54 + 1.0) - 0.15625*z4*(-33.0*nu + 32.0*z18 - 5.0))*(z63 - 1.59227685093395e-9*z64 - 1.67189069348064e-7*z65 + 9.55366110560367e-9*z66 + 1.72773095804465e-13*z67 + 1.72773095804465e-13*z68)**2/(r**13*z69**4*(-0.0438084424460039*nu - 0.143521050466841*r + 0.0185696317637669*z0 + 0.291062041428379*z18 + 0.0210425293255724*z42 + 0.0244692826489756*z43 + 0.0385795738434214*z55 + 0.00941289164152486*z56 + 0.00662650629087394*z57 - 1.18824456940711e-6*z59 + 0.000486339502879429*z60 - 3.63558293513537e-15*z62 + 1)**2))/(z10*(2.0*z28 + 1.0) + 1.0))**(-0.5)*(-pphi*z33*z35*z5 + 2.0*pphi*z9)/(z10*(4.0*z28 + 2.0) + 2.0) + (nu*dSO*z3*z8 + pphi*z12*(z17*z25 + z20*z26 + z24*z27) + pphi*z3*(z15*z17 + z19*z20 + z21*z24) + 0.0416666666666667*z10*z12 + z12*(z25*z30 + z26*z31 + z27*z32 + z28*(0.34375*nu + 0.09375) + z9*(-0.03125*nu + 0.536458333333333*z18 + 0.078125) + 0.25) + z3*z9*(0.0416666666666667*z13*z14 - 0.25*z4) + z3*(z15*z30 + z19*z31 + z21*z32 + z28*(0.71875*nu - 0.09375) + z9*(-5.53125*nu + 0.567708333333333*z18 - 0.078125) + 1.75))/(r*(-2.0*r + z6) + 2.0*z0 + z5)

        # Compute H Jacobian

        omega = M * M * dHeffdpphi / (nu*H)

        return omega

def grad(q, p,chi_1,chi_2,m_1,m_2):

        """
        Compute the gradient of the Hamiltonian in polar coordinates.

        Args:
          q (double[:]): Canonical positions (r,phi).
          p (double[:]): Canonical momenta  (prstar,pphi).
          chi1 (double): Dimensionless z-spin of the primary.
          chi2 (double): Dimensionless z-spin of the secondary.
          m_1 (double): Primary mass component.
          m_2 (double): Secondary mass component.

        Returns:
           (tuple) dHdr, dHdphi, dHdpr, dHdpphi

        """

        # Coordinate definitions

        r = q[0]
        phi = q[1]

        prst = p[0]
        L = p[1]

        pphi = L
        
        M = m_1 + m_2
        X_1 = m_1/M
        X_2 = m_2/M
        nu = X_1*X_2
        
        ap = chi_1*X_1 + chi_2*X_2

        am = chi_1*X_1 - chi_2*X_2

        d5 = 0
        par = array([4.17877875e01, -3.02193382e03, 3.34144394e04, -1.69019140e05, 3.29523262e05])
        a6 = par[0] + par[1] * nu + par[2] * nu ** 2 + par[3] * nu ** 3 + par[4] * nu ** 4
        dSO = (-7.71251231383957 * am ** 3- 17.2294679794015 * am ** 2 * ap- 238.430383378296 * am ** 2 * nu+ 69.5461667822545 * am ** 2- 10.5225438990315 * am * ap ** 2+ 362.767393298729 * am * ap * nu- 85.8036338010274 * am * ap- 1254.66845939312 * am * nu ** 2+ 472.431937787377 * am * nu- 39.742317057316 * am- 7.58458103577458 * ap ** 3- 42.7601129678844 * ap ** 2 * nu+ 18.1783435552183 * ap ** 2- 201.905934468847 * ap * nu ** 2- 90.5790079104259 * ap * nu+ 49.6299175121658 * ap+ 478.546231305475 * nu ** 3+ 679.521769948995 * nu ** 2- 177.334831768076 * nu- 37.6897780220529)

        # Extra quantities used in the Jacobian
        x0 = r**2
        x1 = chi_1*X_1
        x2 = chi_2*X_2
        x3 = x1 + x2
        x4 = x3**2
        x5 = 2.0*x4
        x6 = 2.0*r
        x7 = x0 + x4
        x8 = r*(-x6 + x7)
        x9 = (2.0*x0 + x5 + x8)**(-1)
        x10 = r**4
        x11 = x10**(-1)
        x12 = 3.0*nu
        x13 = pphi*x3
        x14 = r**3
        x15 = x14**(-1)
        x16 = x15*x4
        x17 = x1 - x2
        x18 = X_1 - X_2
        x19 = x17*x18
        x20 = pphi*x19
        x21 = 5.0*X_1 - 5.0*X_2
        x22 = x17*x3
        x23 = 0.0416666666666667*x21*x22 - 0.25*x4
        x24 = 2.0*x15
        x25 = x0**(-1)
        x26 = 0.71875*nu - 0.09375
        x27 = -1.40625*nu - 0.46875
        x28 = pphi**2
        x29 = x15*x28
        x30 = 2.0*x29
        x31 = nu**2
        x32 = -2.0859375*nu - 2.07161458333333*x31 + 0.23046875
        x33 = 3.0*x11
        x34 = x28*x33
        x35 = 0.5859375*nu + 1.34765625*x31 + 0.41015625
        x36 = pphi**4
        x37 = r**5
        x38 = x37**(-1)
        x39 = 4.0*x38
        x40 = x36*x39
        x41 = 0.34375*nu + 0.09375
        x42 = 0.46875 - 0.28125*nu
        x43 = -0.2734375*nu - 0.798177083333333*x31 - 0.23046875
        x44 = -0.3515625*nu + 0.29296875*x31 - 0.41015625
        x45 = nu*dSO*x15
        x46 = x25*x4
        x47 = 0.0416666666666667*x46
        x48 = x23*x25
        x49 = r**(-1)
        x50 = x25*x28
        x51 = x11*x36
        x52 = x3*(x25*(-5.53125*nu + 0.567708333333333*x31 - 0.078125) + x26*x49 + x27*x50 + x29*x32 + x35*x51 + 1.75)
        x53 = x19*(x25*(-0.03125*nu + 0.536458333333333*x31 + 0.078125) + x29*x43 + x41*x49 + x42*x50 + x44*x51 + 0.25)
        x54 = prst**4
        x55 = r**(-4.5)
        x56 = nu*x55
        x57 = prst**6
        x58 = r**(-3.5)
        x59 = nu*x58
        x60 = r**(-2.5)
        x61 = prst**8
        x62 = nu*x61
        x63 = nu*x49
        x64 = 0.121954868780449*x61
        x65 = 8.0*nu - 6.0*x31
        x66 = x25*x65
        x67 = nu**3
        x68 = 92.7110442849544*nu - 131.0*x31 + 10.0*x67
        x69 = x15*x54
        x70 = -2.78300763695006*nu - 5.4*x31 + 6.0*x67
        x71 = x25*x70
        x72 = nu**4
        x73 = -33.9782122170436*nu - 89.5298327361234*x31 + 188.0*x67 - 14.0*x72
        x74 = x15*x57
        x75 = 1.38977750996128*nu + 3.33842023648322*x31 + 3.42857142857143*x67 - 6.0*x72
        x76 = x61*x75
        x77 = log(r)
        x78 = nu*(452.542166996693 - 51.6952380952381*x77) + x31*(118.4*x77 - 1796.13660498019) + 602.318540416564*x67
        x79 = x11*x54
        x80 = r + 2.0
        x81 = x4*x80
        x82 = r*x4
        x83 = x10 + x80*x82
        x84 = x83**(-1)
        x85 = x28*x49*x84
        x86 = r**(-13)
        x87 = x7**4
        x88 = 756.0*nu
        x89 = 336.0*r + x88 + 407.0
        x90 = 2048.0*nu*x77*x89 + 28.0*nu*(1920.0*a6 + 733955.307463037) - 7.0*r*(938918.400156317*nu - 185763.092693281*x31 - 245760.0) - 5416406.59541186*x31 - 3440640.0
        x91 = x77**2
        x92 = x31*x91
        x93 = x67*(-163683964.822551*r - 17833256.898555*x0 - 1188987459.03162)
        x94 = x31*(-39321600.0*a6*(3.0*r + 59.0) + 745857848.115604*a6 + 122635399361.987*r + 2064783811.32587*x0 - 3089250703.76879*x10 - 6178501407.53758*x14 + 1426660551.8844*x37 + 276057889687.011)
        x95 = 588.0*nu + 1079.0
        x96 = x88 + 1079.0
        x97 = x14*x96
        x98 = -38842241.4769507*nu + 240.0*r*(-7466.27061066206*nu - 3024.0*x31 + 17264.0) + 1920.0*x0*x95 + 480.0*x10*x96 - 1882456.23663972*x31 + 161280.0*x37 + 960.0*x97 + 13447680.0
        x99 = nu*x77
        x100 = x98*x99
        x101 = 8.0*r
        x102 = 4.0*x0 + x101 + 2.0*x14 + 16.0
        x103 = 7680.0*a6
        x104 = 128.0*r
        x105 = 7704.0*r
        x106 = 148.04406601634*r
        x107 = 113485.217444961*r
        x108 = nu*(x103*(x10 + x102) + x104*(13218.7851094412*r + 8529.39255472061*x0 - 6852.34813868015*x10 + 4264.6962773603*x14 - 33722.4297811176) + x106*(3852.0*x0 + 349.0*x10 + x105 + 1926.0*x14 + 36400.0) + x107*(-x10 + x102))
        x109 = (32768.0*x100 + 53760.0*x108 + 13212057600.0*x37 + 241555486248.807*x72 + 67645734912.0*x92 + 1120.0*x93 + 7.0*x94)**(-1)
        x110 = x10*x109*x90
        x111 = x110 + 0.000130208333333333*x46
        x112 = x111**(-4)
        x113 = r*x91
        x114 = -630116198.873299*nu - 197773496.793534*x31 + 5805304367.87913
        x115 = x0*x114
        x116 = -2675575.66847905*nu - 138240.0*x31 - 5278341.3229329
        x117 = x116*x14
        x118 = nu*(-2510664218.28128*nu - 42636451.6032331*x31 + 14515200.0*x67 + 1002013764.01019)
        x119 = 43393301259014.8*nu + 43133561885859.3*x31 + 5927865218923.02*x67 + 86618264430493.3*(1 - 0.496948781616935*nu)**2 + 188440788778196.0
        x120 = r*x119
        x121 = 14700.0*nu + 42911.0
        x122 = 283115520.0*x121
        x123 = -1698693120.0*nu*(11592.0*nu + 69847.0) + 49152.0*r*(409207698.136075*nu + 102574080.0*x31 - 2119671837.36038) + x0*x122 + 879923036160.0*x14
        x124 = x123*x77
        x125 = (x113 - 1.59227685093395e-9*x115 - 1.67189069348064e-7*x117 + 9.55366110560367e-9*x118 + 1.72773095804465e-13*x120 + 1.72773095804465e-13*x124)**2
        x126 = nu*r
        x127 = nu*x0
        x128 = r*x31
        x129 = r*x121
        x130 = x0*x31
        x131 = 5041721180160.0*x31 - 104186110149937.0
        x132 = -25392914995744.3*nu - r*x122 - 879923036160.0*x0 - x131
        x133 = x132*x77
        x134 = -0.0438084424460039*nu - 0.143521050466841*r + 0.0185696317637669*x0 + 0.0385795738434214*x126 + 0.00941289164152486*x127 + 0.00662650629087394*x128 - 1.18824456940711e-6*x129 + 0.000486339502879429*x130 - 3.63558293513537e-15*x133 + 0.291062041428379*x31 + 0.0244692826489756*x67 + 0.0210425293255724*x91 + 1
        x135 = x134**(-2)
        x136 = 4.0*x31
        x137 = x17**2
        x138 = -0.46875*x137*(-5.0*nu + x136 + 1.0) + 0.0625*x17*x21*x3*(18.0*nu - 1.0) - 0.15625*x4*(-33.0*nu + 32.0*x31 - 5.0)
        x139 = x112*x125*x135*x138*x54*x87
        x140 = x7**2
        x141 = prst**2
        x142 = x111**(-2)
        x143 = 1822680546449.21*x31
        x144 = 5787938193408.0*x91
        x145 = -12049908701745.2*nu + r*x143 - 39476764256925.6*r + 5107745331375.71*x0 + 10611661054566.2*x126 + 2589101062873.81*x127 - 326837426.241486*x129 + 133772083200.0*x130 - x133 + x144 + 80059249540278.2*x31 + 6730497718123.02*x67 + 275059053208689.0
        x146 = x145**(-1)
        x147 = 0.0625*x137
        x148 = 0.125*x22
        x149 = -1171.0*nu - 861.0
        x150 = 0.015625*x4
        x151 = x137*(115.0*nu + x136 - 37.0)
        x152 = 0.03125*x22
        x153 = x18*(26.0*nu + 449.0)
        x154 = 5787938193408.0*x113 - 9216.0*x115 - 967680.0*x117 + 55296.0*x118 + x120 + x124
        x155 = x154**(-1)
        x156 = x145*x155*x37
        x157 = x109*x90
        x158 = x11*(x149*x150 + 0.015625*x151 + x152*x153) + x15*(x147*(12.0*nu - 3.0) + x148*(-21.0*X_1 + 21.0*X_2) + x4*(x12 + 2.8125)) + 7680.0*x156*x157 + x46
        x159 = x140*x141*x142*x146*x154*x158
        x160 = 1.27277314139085e-19*x139*x86 + 1.69542100694444e-8*x159*x38 + x25*x76 + x50 + 147.443752990146*x54*x56 + x54*x66 - 11.3175085791863*x57*x59 + x57*x71 + 1.48275342024365*x60*x62 + x63*x64 + x68*x69 + x73*x74 + x78*x79 - x81*x85 + 1.0
        x161 = x46*(2.0*x49 + 1.0) + 1.0
        x162 = x161**(-1)
        x163 = x137*(4.0*nu + 1.0)
        x164 = -x21*x22
        x165 = x137*(-27.0*nu + 28.0*x31 - 3.0)
        x166 = x152*(-39.0*X_1 + 39.0*X_2)
        x167 = x11*(0.125*x163 + 0.25*x164 + 1.125*x4) + 7680.0*x110 + x38*(-x150*(175.0*nu + 225.0) + 0.046875*x165 + x166*(2.0*nu - 3.0)) + x46
        x168 = x162*x167
        x169 = (x160*x168)**(-0.5)
        x170 = 4.0*x49 + 2.0
        x171 = r**(-6)
        x172 = x15*x5
        x173 = -6572428.80109422*nu + 1300341.64885296*x31 + 2048.0*x63*x89 + 688128.0*x99 + 1720320.0
        x174 = x31*x49
        x175 = 4.0*x14
        x176 = 6.0*x0 + x101 + 8.0
        x177 = 1.31621673590926e-19*x90*(53760.0*nu*(3740417.71815805*r + 2115968.85907902*x0 - 938918.400156317*x10 + x103*(x175 + x176) + x104*(17058.7851094412*r + 12794.0888320809*x0 - 27409.3925547206*x14 + 13218.7851094412) + x106*(5778.0*x0 + x105 + 1396.0*x14 + 7704.0) + x107*(-x175 + x176) + 1057984.42953951*x14 + 2888096.47013111) + 66060288000.0*x10 + 135291469824.0*x174*x77 + 7.0*x31*(-117964800.0*a6 + 4129567622.65173*r - 18535504222.6128*x0 + 7133302759.42198*x10 - 12357002815.0752*x14 + 122635399361.987) + 32768.0*x63*x98 + 7.0*x67*(-5706642207.53758*r - 26189434371.6082) + 32768.0*x99*(-1791904.9465589*nu + 3840.0*r*x95 + 2880.0*x0*x96 + 806400.0*x10 - 725760.0*x31 + 1920.0*x97 + 4143360.0))/(1.35654132757922e-7*x100 + 2.22557561555966e-7*x108 + 0.0546957463279941*x37 + x72 + 0.28004222119933*x92 + 4.63661586574928e-9*x93 + 2.8978849160933e-11*x94)**2
        x178 = -7680.0*x10*x109*x173 + x10*x177 - 30720.0*x109*x14*x90 + x172
        x179 = 11575876386816.0*x77
        x180 = 5807150888816.34*nu + 10215490662751.4*r + 5178202125747.62*x126 + 267544166400.0*x128 - x132*x49 + x143 + x179*x49 + x77*(4161798144000.0*nu + 1759846072320.0*r + 12148770078720.0) - 53501685054374.1
        x181 = x125*x138*x54*x86*x87
        x182 = -18432.0*r*x114 - 2903040.0*x0*x116 + x119 + x123*x49 + x144 + x179 + x77*(20113376778784.3*nu + 2639769108480.0*x0 + 566231040.0*x129 + x131)
        x183 = x141*x158
        x184 = -x178
        x185 = x154*x38
        x186 = prst**3
        x187 = prst**5
        x188 = prst**7
        x189 = 4.0*x186
        x190 = 6.0*x187
        x191 = x167*x169/(x170*x46 + 2.0)
        x192 = 2.0*pphi
        x193 = x192*x25
        x194 = x15*x192
        x195 = 4.0*pphi**3*x11


        # Evaluate Hamiltonian
        H,xi =  evaluate_H(q,p,chi_1,chi_2,m_1,m_2)

        # Heff Jacobian expressions
        dHeffdr = 0.5*x169*(x160*x162*(-x171*(-x150*(875.0*nu + 1125.0) + 0.234375*x165 + x166*(10.0*nu - 15.0)) - x178 - x38*(0.5*x163 + x164 + 4.5*x4)) - x160*x167*(-x11*x5 - x16*x170)/x161**2 + x168*(-663.496888455656*nu*r**(-5.5)*x54 - nu*x25*x64 + 39.6112800271521*nu*x55*x57 + 6.78168402777778e-8*x11*x141*x142*x146*x154*x158*x7 + x11*x54*(118.4*x174 - 51.6952380952381*x63) + 7.59859378406358e-45*x112*x135*x138*x154*x182*x54*x86*x87 - 9.25454462627843e-34*x112*x180*x181/x134**3 - 2.24091649004576e-37*x135*x140*x142*x154*x180*x183*x38 + 1.69542100694444e-8*x140*x141*x142*x146*x154*x38*(38400.0*x10*x109*x145*x155*x90 + 7680.0*x109*x145*x155*x173*x37 + 7680.0*x109*x155*x180*x37*x90 - x11*(x147*(36.0*nu - 9.0) + x148*(-63.0*X_1 + 63.0*X_2) + x4*(9.0*nu + 8.4375)) - x156*x177 - x172 - x38*(x148*x153 + 0.0625*x149*x4 + 0.0625*x151) - 2.29252167428035e-22*x145*x157*x182*x37/x125) + 1.69542100694444e-8*x140*x141*x142*x146*x158*x182*x38 - 8.47710503472222e-8*x159*x171 - x24*x76 + x25*x28*x4*x80*x84 + x28*x4*x49*x80*(x175 + x81 + x82)/x83**2 - x30 - x33*x57*x73 - x39*x54*x78 - x4*x85 - 3.70688355060912*x58*x62 - 2.0*x65*x69 - 3.0*x68*x79 - 2.0*x70*x74 - 4.41515887225116e-12*x140*x146*x183*x184*x185/x111**3 - 6.62902677807736e-23*x135*x181*x184/x111**5 + 1.01821851311268e-18*x112*x125*x135*x138*x54*x7**3/r**12 - 1.65460508380811e-18*x139/r**14)) + x9*(-dSO*x11*x12*x13 + pphi*x17*x18*(-x15*(-0.0625*nu + 1.07291666666667*x31 + 0.15625) - x25*x41 - x30*x42 - x34*x43 - x40*x44) + pphi*x3*(-x15*(-11.0625*nu + 1.13541666666667*x31 - 0.15625) - x25*x26 - x27*x30 - x32*x34 - x35*x40) - x13*x23*x24 - 0.0833333333333333*x16*x20) - 0.25*(r*(x6 - 2.0) + x6 + x7)*(pphi*x52 + pphi*x53 + x13*x45 + x13*x48 + x20*x47)/(x7 + 0.5*x8)**2

        dHeffdpr = x191*(11.8620273619492*nu*x188*x60 + 3.39084201388889e-8*prst*x140*x142*x146*x158*x185 + x11*x189*x78 + 5.09109256556341e-19*x112*x125*x135*x138*x186*x86*x87 + x15*x189*x68 + x15*x190*x73 + 589.775011960583*x186*x56 - 67.9050514751178*x187*x59 + 8.0*x188*x25*x75 + 0.975638950243592*x188*x63 + x189*x66 + x190*x71)

        dHeffdpphi = x191*(2.0*pphi*x25 - pphi*x49*x5*x80*x84) + x9*(x13*(x193*x27 + x194*x32 + x195*x35) + x19*x47 + x20*(x193*x42 + x194*x43 + x195*x44) + x3*x45 + x3*x48 + x52 + x53)

        # Compute H Jacobian
        dHdr = M * M * dHeffdr / (nu*H)
        dHdpr = M * M * dHeffdpr / (nu*H)
        dHdpphi = M * M * dHeffdpphi / (nu*H)

        return dHdr, dHdpr, dHdpphi

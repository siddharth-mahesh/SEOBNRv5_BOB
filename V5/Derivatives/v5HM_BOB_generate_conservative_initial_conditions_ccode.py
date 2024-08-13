import sympy as sp
from nrpy.c_codegen import c_codegen
def v5HM_BOB_generate_conservative_initial_conditions_ccode():
    m1, m2, r, pphi, chi1, chi2, a6, dSO = sp.symbols('m1 m2 r pphi chi1 chi2 a6 dSO',real = True)
    prstar = 0.
    u = 1/r

    M = m1+m2

    delta = (m1-m2)/M

    nu = m1*m2/(M**2)

    gam = (sp.Rational(1,4)+(pphi**2/r**2)*(sp.Rational(15,32)-sp.Rational(9,32)*nu)+(1/r)*(sp.Rational(11,32)*nu+sp.Rational(3,32))+(pphi**4/r**4)*(sp.Rational(75,256)*nu**2-sp.Rational(45,128)*nu-sp.Rational(105,256))+(pphi**2/r**3)*(-sp.Rational(613,768)*nu**2-sp.Rational(35,128)*nu-sp.Rational(59,256))+(1/r**2)*(sp.Rational(103,192)*nu**2-sp.Rational(1,32)*nu+sp.Rational(5,64)))

    gap = (sp.Rational(7,4)+(pphi**2/r**2)*(-sp.Rational(45,32)*nu-sp.Rational(15,32))+(1/r)*(sp.Rational(23,32)*nu-sp.Rational(3,32))+(pphi**4/r**4)*(sp.Rational(345,256)*nu**2+sp.Rational(75,128)*nu+sp.Rational(105,256))+(pphi**2/r**3)*(-sp.Rational(1591,768)*nu**2-sp.Rational(267,128)*nu+sp.Rational(59,256))+(1/r**2)*(sp.Rational(109,192)*nu**2-sp.Rational(177,32)*nu-sp.Rational(5,64)))

    am = (m1*chi1-m2*chi2)/M

    ap = (m1*chi1+m2*chi2)/M

    Qnos = (0.121954868780449*nu*prstar**8/r+prstar**6*(6.0*nu**3-5.4*nu**2-2.78300763695006*nu)/r**2+prstar**4*(10.0*nu**3-131.0*nu**2+92.7110442849544*nu)/r**3)+(prstar**8*(-6.0*nu**4+3.42857142857143*nu**3+3.33842023648322*nu**2+1.38977750996128*nu)/r**2+prstar**6*(-14.0*nu**4+188.0*nu**3-89.5298327361234*nu**2-33.9782122170436*nu)/r**3+prstar**4*(602.318540416564*nu**3+nu**2*(118.4*sp.log(r)-1796.13660498019)+nu*(452.542166996693-51.6952380952381*sp.log(r)))/r**4)+(1.48275342024365*nu*prstar**8/r**2.5-11.3175085791863*nu*prstar**6/r**3.5+147.443752990146*nu*prstar**4/r**4.5)+prstar**4*(-6.0*nu**2+8.0*nu)/r**2

    d5 = 0

    Dnons = r*(6730497718123.02*nu**3+22295347200.0*nu**2*d5+133772083200.0*nu**2*r**2+1822680546449.21*nu**2*r+80059249540278.2*nu**2+22295347200.0*nu*d5*r-193226342400.0*nu*d5+2589101062873.81*nu*r**2+10611661054566.2*nu*r-12049908701745.2*nu+5107745331375.71*r**2-326837426.241486*r*(14700.0*nu+42911.0)-39476764256925.6*r-(-5041721180160.0*nu**2-25392914995744.3*nu-879923036160.0*r**2-283115520.0*r*(14700.0*nu+42911.0)+104186110149937.0)*sp.log(r)+5787938193408.0*sp.log(r)**2+275059053208689.0)/(55296.0*nu*(14515200.0*nu**3-42636451.6032331*nu**2-7680.0*nu*(315.0*d5+890888.810272497)+4331361844.61149*nu+1002013764.01019)-967680.0*r**3*(-138240.0*nu**2-2675575.66847905*nu-5278341.3229329)-9216.0*r**2*(-197773496.793534*nu**2-7680.0*nu*(315.0*d5+405152.309729121)+2481453539.84635*nu+5805304367.87913)+r*(5927865218923.02*nu**3+70778880.0*nu**2*(315.0*d5+2561145.80918574)-138141470005001.0*nu**2-4718592.0*nu*(40950.0*d5+86207832.4415642)+450172889755120.0*nu+86618264430493.3*(1-0.496948781616935*nu)**2+188440788778196.0)+5787938193408.0*r*sp.log(r)**2+(-1698693120.0*nu*(11592.0*nu+69847.0)+879923036160.0*r**3+283115520.0*r**2*(14700.0*nu+42911.0)+49152.0*r*(102574080.0*nu**2+409207698.136075*nu-2119671837.36038))*sp.log(r))

    Anons = 7680.0*r**4*(-5416406.59541186*nu**2+28.0*nu*(1920.0*a6+733955.307463037)+2048.0*nu*(756.0*nu+336.0*r+407.0)*sp.log(r)-7.0*r*(-185763.092693281*nu**2+938918.400156317*nu-245760.0)-3440640.0)/(241555486248.807*nu**4+1120.0*nu**3*(-17833256.898555*r**2-163683964.822551*r-1188987459.03162)+7.0*nu**2*(-39321600.0*a6*(3.0*r+59.0)+745857848.115604*a6+1426660551.8844*r**5-3089250703.76879*r**4-6178501407.53758*r**3+2064783811.32587*r**2+122635399361.987*r+276057889687.011)+67645734912.0*nu**2*sp.log(r)**2+53760.0*nu*(7680.0*a6*(r**4+2.0*r**3+4.0*r**2+8.0*r+16.0)+128.0*r*(-6852.34813868015*r**4+4264.6962773603*r**3+8529.39255472061*r**2+13218.7851094412*r-33722.4297811176)+113485.217444961*r*(-r**4+2.0*r**3+4.0*r**2+8.0*r+16.0)+148.04406601634*r*(349.0*r**4+1926.0*r**3+3852.0*r**2+7704.0*r+36400.0))+32768.0*nu*(-1882456.23663972*nu**2-38842241.4769507*nu+161280.0*r**5+480.0*r**4*(756.0*nu+1079.0)+960.0*r**3*(756.0*nu+1079.0)+1920.0*r**2*(588.0*nu+1079.0)+240.0*r*(-3024.0*nu**2-7466.27061066206*nu+17264.0)+13447680.0)*sp.log(r)+13212057600.0*r**5)

    xi = sp.sqrt(Dnons)*(Anons+ap**2*u*u)/(1+ap**2*u*u)

    pr = prstar/xi

    QalignSS = ((pr**4)/(r**3))*(ap**2*(-5*nu*nu+nu*sp.Rational(165,32)+sp.Rational(25,32))+delta*ap*am*(nu*sp.Rational(45,8)-sp.Rational(5,16))+am**2*(-nu*nu*sp.Rational(15,8)+nu*sp.Rational(75,32)-sp.Rational(15,32)))

    BnpalignSS = ((1/r**3)*(ap**2*(3*nu+sp.Rational(45,16))-delta*ap*am*sp.Rational(21,8)+am**2*(nu*sp.Rational(3,4)-sp.Rational(3,16)))+(1/r**4)*(ap**2*(-nu*sp.Rational(1171,64)-sp.Rational(861,64))+delta*ap*am*(nu*sp.Rational(13,16)+sp.Rational(449,32))+am**2*(nu*nu*sp.Rational(1,16)+nu*sp.Rational(115,64)-sp.Rational(37,64))))

    AalignSS = ((1/r**4)*(ap**2*sp.Rational(9,8)-delta*ap*am*sp.Rational(5,4)+am**2*(nu*sp.Rational(1,2)+sp.Rational(1,8)))+(1/r**5)*(ap**2*(-nu*sp.Rational(175,64)-sp.Rational(225,64))+delta*ap*am*(-nu*sp.Rational(39,16)+sp.Rational(117,32))+am**2*(nu*nu*sp.Rational(21,16)-nu*sp.Rational(81,64)-sp.Rational(9,64))))

    Qalign = Qnos+QalignSS

    Balignnp = -1+ap**2*u*u+Anons*Dnons+BnpalignSS

    Bkerreqnp = -(1+2/r)/(r**2+ap**2*(1+2/r))

    Aalign = (ap**2*u*u+Anons+AalignSS)/(1+ap**2*(1+2/r)/(r**2))

    Galigna3 = pphi*(delta*am*ap**2-ap**3)/(4*r**2)

    SOcalib = nu*dSO*ap*pphi*(u**3)

    Heven = sp.sqrt(Aalign*(1+pphi*pphi/r**2+(1+Balignnp)*pr*pr+Bkerreqnp*pphi*pphi*ap**2/r**2+Qalign))

    Hodd = (pphi*(gap*ap+delta*gam*am)+SOcalib+Galigna3)/(r**3+ap**2*(r+2))

    Heff = Hodd+Heven

    Hreal = sp.sqrt(1+2*nu*(Heff-1))

    gam_prmpphi = pphi**3*(75*nu**2/64 - 45*nu/32 - 105/64)/r**4 + pphi*(15/16 - 9*nu/16)/r**2 + pphi*(-613*nu**2/384 - 35*nu/64 - 59/128)/r**3
    gap_prmpphi = pphi**3*(345*nu**2/64 + 75*nu/32 + 105/64)/r**4 + pphi*(-45*nu/16 - 15/16)/r**2 + pphi*(-1591*nu**2/384 - 267*nu/64 + 59/128)/r**3
    Galigna3_prmpphi = (am*ap**2*delta - ap**3)/(4*r**2)
    SOcalib_prmpphi = ap*dSO*nu*u**3
    Heven_prmpphi = sp.sqrt(Aalign*(Bkerreqnp*ap**2*pphi**2/r**2 + Qalign + pphi**2/r**2 + pr**2*(Balignnp + 1) + 1))*(2*Bkerreqnp*ap**2*pphi/r**2 + 2*pphi/r**2)/(2*(Bkerreqnp*ap**2*pphi**2/r**2 + Qalign + pphi**2/r**2 + pr**2*(Balignnp + 1) + 1))
    Hodd_prmpphi = (Galigna3_prmpphi + SOcalib_prmpphi + am*delta*gam + ap*gap + pphi*(am*delta*gam_prmpphi + ap*gap_prmpphi))/(ap**2*(r + 2) + r**3)
    Heff_prmpphi = Heven_prmpphi + Hodd_prmpphi
    Hreal_prmpphi = Heff_prmpphi*nu/sp.sqrt(nu*(2*Heff - 2) + 1)
    u_prmr = -1/r**2
    gam_prmr = -pphi**4*(75*nu**2/64 - 45*nu/32 - 105/64)/r**5 + pphi**2*(9*nu/16 - 15/16)/r**3 - pphi**2*(-613*nu**2/256 - 105*nu/128 - 177/256)/r**4 - (11*nu/32 + 3/32)/r**2 - (103*nu**2/96 - nu/16 + 5/32)/r**3
    gap_prmr = -pphi**4*(345*nu**2/64 + 75*nu/32 + 105/64)/r**5 - pphi**2*(-45*nu/16 - 15/16)/r**3 - pphi**2*(-1591*nu**2/256 - 801*nu/128 + 177/256)/r**4 - (23*nu/32 - 3/32)/r**2 - (109*nu**2/96 - 177*nu/16 - 5/32)/r**3
    Qnos_prmr = -3.70688355060913*nu*prstar**8/r**3.5 - 0.121954868780449*nu*prstar**8/r**2 + 39.6112800271521*nu*prstar**6/r**4.5 - 663.496888455657*nu*prstar**4/r**5.5 - prstar**8*(-12.0*nu**4 + 6.85714285714286*nu**3 + 6.67684047296644*nu**2 + 2.77955501992256*nu)/r**3 - prstar**6*(12.0*nu**3 - 10.8*nu**2 - 5.56601527390012*nu)/r**3 - prstar**6*(-42.0*nu**4 + 564.0*nu**3 - 268.58949820837*nu**2 - 101.934636651131*nu)/r**4 - prstar**4*(-12.0*nu**2 + 16.0*nu)/r**3 + prstar**4*(118.4*nu**2/r - 51.6952380952381*nu/r)/r**4 - prstar**4*(30.0*nu**3 - 393.0*nu**2 + 278.133132854863*nu)/r**4 - prstar**4*(2409.27416166626*nu**3 + 4*nu**2*(118.4*sp.log(r) - 1796.13660498019) + 4*nu*(452.542166996693 - 51.6952380952381*sp.log(r)))/r**5
    Dnons_prmr = r*(22295347200.0*d5*nu + 267544166400.0*nu**2*r + 1822680546449.21*nu**2 + 5178202125747.62*nu*r + 5807150888816.35*nu + 10215490662751.4*r - (-4161798144000.0*nu - 1759846072320.0*r - 12148770078720.0)*sp.log(r) - 53501685054374.0 - (-5041721180160.0*nu**2 - 25392914995744.3*nu - 879923036160.0*r**2 + r*(-4161798144000.0*nu - 12148770078720.0) + 104186110149937.0)/r + 11575876386816.0*sp.log(r)/r)/(nu*(802632499200.0*nu**3 - 2357625227852.38*nu**2 + 55296.0*nu*(-2419200.0*d5 - 6842026062.89278) + 239506984559637.0*nu + 55407353094707.5) + r**3*(133772083200.0*nu**2 + 2589101062873.81*nu + 5107745331375.71) - r**2*(-1822680546449.21*nu**2 + 9216.0*nu*(-2419200.0*d5 - 3111569738.71965) + 22869075823224.0*nu + 53501685054374.1) + r*(5927865218923.02*nu**3 + nu**2*(22295347200.0*d5 + 181275031890860.0) - 138141470005001.0*nu**2 - nu*(193226342400.0*d5 + 406779588496105.0) + 450172889755120.0*nu + 86618264430493.3*(1 - 0.496948781616935*nu)**2 + 188440788778196.0) + 5787938193408.0*r*sp.log(r)**2 + (nu*(-19691250647040.0*nu - 118648618352640.0) + 879923036160.0*r**3 + r**2*(4161798144000.0*nu + 12148770078720.0) + r*(5041721180160.0*nu**2 + 20113376778784.4*nu - 104186110149937.0))*sp.log(r)) + r*(-1.7694999344523e-13*nu**3 - 2.98505426338587e-26*nu**2*(22295347200.0*d5 + 181275031890860.0) + 4.123597839888195e-12*nu**2 + 2.98505426338587e-26*nu*(193226342400.0*d5 + 406779588496105.0) - 1.343790503824258e-11*nu + 2.98505426338587e-26*r**2*(-401316249600.0*nu**2 - 7767303188621.42*nu - 15323235994127.1) + 2.98505426338587e-26*r*(-3645361092898.42*nu**2 + 18432.0*nu*(-2419200.0*d5 - 3111569738.71965) + 45738151646447.9*nu + 107003370108748.0) - 2.58560219525329e-12*(1 - 0.496948781616935*nu)**2 - 2.98505426338587e-26*(5041721180160.0*nu**2 + 20113376778784.4*nu + 2639769108480.0*r**2 + r*(8323596288000.0*nu + 24297540157440.0) - 104186110149937.0)*sp.log(r) - 1.72773095804465e-13*sp.log(r)**2 - 3.45546191608929e-13*sp.log(r) - 5.625059799381501e-12 - 2.98505426338587e-26*(nu*(-19691250647040.0*nu - 118648618352640.0) + 879923036160.0*r**3 + r**2*(4161798144000.0*nu + 12148770078720.0) + r*(5041721180160.0*nu**2 + 20113376778784.4*nu - 104186110149937.0))/r)*(22295347200.0*d5*nu**2 + 22295347200.0*d5*nu*r - 193226342400.0*d5*nu + 6730497718123.02*nu**3 + 133772083200.0*nu**2*r**2 + 1822680546449.21*nu**2*r + 80059249540278.2*nu**2 + 2589101062873.81*nu*r**2 + 10611661054566.2*nu*r - 12049908701745.2*nu + 5107745331375.71*r**2 + r*(-4804510165749.84*nu - 14024920797448.4) - 39476764256925.6*r - (-5041721180160.0*nu**2 - 25392914995744.3*nu - 879923036160.0*r**2 + r*(-4161798144000.0*nu - 12148770078720.0) + 104186110149937.0)*sp.log(r) + 5787938193408.0*sp.log(r)**2 + 275059053208689.0)/(nu*(0.138673301680058*nu**3 - 0.407334209362761*nu**2 + 9.55366110560367e-9*nu*(-2419200.0*d5 - 6842026062.89278) + 41.3803631891606*nu + 9.57289992450369) + r**3*(0.0231122169466764*nu**2 + 0.447327005983337*nu + 0.882480973482581) - r**2*(-0.314910160672604*nu**2 + 1.59227685093395e-9*nu*(-2419200.0*d5 - 3111569738.71965) + 3.95116102816545*nu + 9.24365175759969) + r*(1.02417562538494*nu**3 + 1.72773095804465e-13*nu**2*(22295347200.0*d5 + 181275031890860.0) - 23.86712943174367*nu**2 - 1.72773095804465e-13*nu*(193226342400.0*d5 + 406779588496105.0) + 77.77776381023421*nu + 14.9653056988661*(1 - 0.496948781616935*nu)**2 + 32.55749845304421) + r*sp.log(r)**2 + (1.72773095804465e-13*nu*(-19691250647040.0*nu - 118648618352640.0) + 0.152027027027027*r**3 + 1.72773095804465e-13*r**2*(4161798144000.0*nu + 12148770078720.0) + 1.72773095804465e-13*r*(5041721180160.0*nu**2 + 20113376778784.4*nu - 104186110149937.0))*sp.log(r))**2 + (22295347200.0*d5*nu**2 + 22295347200.0*d5*nu*r - 193226342400.0*d5*nu + 6730497718123.02*nu**3 + 133772083200.0*nu**2*r**2 + 1822680546449.21*nu**2*r + 80059249540278.2*nu**2 + 2589101062873.81*nu*r**2 + 10611661054566.2*nu*r - 12049908701745.2*nu + 5107745331375.71*r**2 + r*(-4804510165749.84*nu - 14024920797448.4) - 39476764256925.6*r - (-5041721180160.0*nu**2 - 25392914995744.3*nu - 879923036160.0*r**2 + r*(-4161798144000.0*nu - 12148770078720.0) + 104186110149937.0)*sp.log(r) + 5787938193408.0*sp.log(r)**2 + 275059053208689.0)/(nu*(802632499200.0*nu**3 - 2357625227852.38*nu**2 + 55296.0*nu*(-2419200.0*d5 - 6842026062.89278) + 239506984559637.0*nu + 55407353094707.5) + r**3*(133772083200.0*nu**2 + 2589101062873.81*nu + 5107745331375.71) - r**2*(-1822680546449.21*nu**2 + 9216.0*nu*(-2419200.0*d5 - 3111569738.71965) + 22869075823224.0*nu + 53501685054374.1) + r*(5927865218923.02*nu**3 + nu**2*(22295347200.0*d5 + 181275031890860.0) - 138141470005001.0*nu**2 - nu*(193226342400.0*d5 + 406779588496105.0) + 450172889755120.0*nu + 86618264430493.3*(1 - 0.496948781616935*nu)**2 + 188440788778196.0) + 5787938193408.0*r*sp.log(r)**2 + (nu*(-19691250647040.0*nu - 118648618352640.0) + 879923036160.0*r**3 + r**2*(4161798144000.0*nu + 12148770078720.0) + r*(5041721180160.0*nu**2 + 20113376778784.4*nu - 104186110149937.0))*sp.log(r))
    Anons_prmr = r**4*(-7.12916500937039e-13*nu**2 + 1.31621673590926e-19*nu*(53760.0*a6 + 20550748.608965) + 1.31621673590926e-19*nu*(1548288.0*nu + 688128.0*r + 833536.0)*sp.log(r) - 1.31621673590926e-19*r*(-1300341.64885297*nu**2 + 6572428.80109422*nu - 1720320.0) - 4.52862795023884e-13)*(nu**3*(39946495452.7632*r + 183326040601.257) - nu**2*(-825753600.0*a6 + 49933119315.954*r**4 - 86499019705.5264*r**3 - 129748529558.289*r**2 + 28906973358.5622*r + 858447795533.909) - 135291469824.0*nu**2*sp.log(r)/r - nu*(-23781703680.0*nu**2 - 58717141288.8419*nu + 26424115200.0*r**4 + 32768.0*r**3*(1451520.0*nu + 2071680.0) + 32768.0*r**2*(2177280.0*nu + 3107520.0) + 32768.0*r*(2257920.0*nu + 4143360.0) + 135769620480.0)*sp.log(r) - nu*(412876800.0*a6*(4*r**3 + 6.0*r**2 + 8.0*r + 8.0) - 50476253192.4036*r**4 + 56877242932.0441*r**3 + 113754485864.088*r**2 + 53760.0*r*(-3508402.24700424*r**3 + 1637643.37050636*r**2 + 2183524.49400847*r + 1692004.49400847) + 53760.0*r*(-453940.869779844*r**3 + 680911.304669766*r**2 + 907881.739559688*r + 907881.739559688) + 53760.0*r*(206669.516158811*r**3 + 855398.613442412*r**2 + 1140531.48458988*r + 1140531.48458988) + 201084856528.177*r + 155264066234.248) - nu*(-61684325962.2103*nu**2 - 1272782568716.72*nu + 5284823040.0*r**5 + 32768.0*r**4*(362880.0*nu + 517920.0) + 32768.0*r**3*(725760.0*nu + 1035840.0) + 32768.0*r**2*(1128960.0*nu + 2071680.0) + 32768.0*r*(-725760.0*nu**2 - 1791904.94655889*nu + 4143360.0) + 440653578240.0)/r - 66060288000.0*r**4)/(nu**4 + nu**3*(-0.0826859618738229*r**2 - 0.758939668264987*r - 5.51287811672293) + nu**2*(-0.00113949471516654*a6*(3.0*r + 59.0) + 0.0216141020760402*a6 + 0.0413429809369115*r**5 - 0.0895229301648219*r**4 - 0.179045860329644*r**3 + 0.0598350586183487*r**2 + 3.5538327399018*r + 7.99983994492537) + 0.280042221199329*nu**2*sp.log(r)**2 + nu*(0.00170924207274981*a6*(r**4 + 2.0*r**3 + 4.0*r**2 + 8.0*r + 16.0) + 2.22557561555965e-7*r*(-877100.561751059*r**4 + 545881.123502118*r**3 + 1091762.24700424*r**2 + 1692004.49400847*r - 4316471.01198305) + 2.22557561555965e-7*r*(-113485.217444961*r**4 + 226970.434889922*r**3 + 453940.869779844*r**2 + 907881.739559688*r + 1815763.47911938) + 2.22557561555965e-7*r*(51667.3790397027*r**4 + 285132.871147471*r**3 + 570265.742294942*r**2 + 1140531.48458988*r + 5388804.00299478)) + nu*(-0.255362968236103*nu**2 - 5.26911058192953*nu + 0.0218782985311977*r**5 + 1.35654132757922e-7*r**4*(362880.0*nu + 517920.0) + 1.35654132757922e-7*r**3*(725760.0*nu + 1035840.0) + 1.35654132757922e-7*r**2*(1128960.0*nu + 2071680.0) + 1.35654132757922e-7*r*(-725760.0*nu**2 - 1791904.94655889*nu + 4143360.0) + 1.82423336800605)*sp.log(r) + 0.054695746327994*r**5)**2 + r**4*(9986623863.19079*nu**2 + 5284823040.0*nu*sp.log(r) - 50476253192.4036*nu + 7680.0*nu*(1548288.0*nu + 688128.0*r + 833536.0)/r + 13212057600.0)/(241555486248.807*nu**4 + nu**3*(-19973247726.3816*r**2 - 183326040601.257*r - 1331665954115.41) + nu**2*(-275251200.0*a6*(3.0*r + 59.0) + 5221004936.80923*a6 + 9986623863.1908*r**5 - 21624754926.3815*r**4 - 43249509852.7631*r**3 + 14453486679.2811*r**2 + 858447795533.909*r + 1932405227809.08) + 67645734912.0*nu**2*sp.log(r)**2 + nu*(412876800.0*a6*(r**4 + 2.0*r**3 + 4.0*r**2 + 8.0*r + 16.0) + 53760.0*r*(-877100.561751059*r**4 + 545881.123502118*r**3 + 1091762.24700424*r**2 + 1692004.49400847*r - 4316471.01198305) + 53760.0*r*(-113485.217444961*r**4 + 226970.434889922*r**3 + 453940.869779844*r**2 + 907881.739559688*r + 1815763.47911938) + 53760.0*r*(51667.3790397027*r**4 + 285132.871147471*r**3 + 570265.742294942*r**2 + 1140531.48458988*r + 5388804.00299478)) + nu*(-61684325962.2103*nu**2 - 1272782568716.72*nu + 5284823040.0*r**5 + 32768.0*r**4*(362880.0*nu + 517920.0) + 32768.0*r**3*(725760.0*nu + 1035840.0) + 32768.0*r**2*(1128960.0*nu + 2071680.0) + 32768.0*r*(-725760.0*nu**2 - 1791904.94655889*nu + 4143360.0) + 440653578240.0)*sp.log(r) + 13212057600.0*r**5) + r**3*(-166392010611.052*nu**2 + 30720.0*nu*(53760.0*a6 + 20550748.608965) + 30720.0*nu*(1548288.0*nu + 688128.0*r + 833536.0)*sp.log(r) - 30720.0*r*(-1300341.64885297*nu**2 + 6572428.80109422*nu - 1720320.0) - 105696460800.0)/(241555486248.807*nu**4 + nu**3*(-19973247726.3816*r**2 - 183326040601.257*r - 1331665954115.41) + nu**2*(-275251200.0*a6*(3.0*r + 59.0) + 5221004936.80923*a6 + 9986623863.1908*r**5 - 21624754926.3815*r**4 - 43249509852.7631*r**3 + 14453486679.2811*r**2 + 858447795533.909*r + 1932405227809.08) + 67645734912.0*nu**2*sp.log(r)**2 + nu*(412876800.0*a6*(r**4 + 2.0*r**3 + 4.0*r**2 + 8.0*r + 16.0) + 53760.0*r*(-877100.561751059*r**4 + 545881.123502118*r**3 + 1091762.24700424*r**2 + 1692004.49400847*r - 4316471.01198305) + 53760.0*r*(-113485.217444961*r**4 + 226970.434889922*r**3 + 453940.869779844*r**2 + 907881.739559688*r + 1815763.47911938) + 53760.0*r*(51667.3790397027*r**4 + 285132.871147471*r**3 + 570265.742294942*r**2 + 1140531.48458988*r + 5388804.00299478)) + nu*(-61684325962.2103*nu**2 - 1272782568716.72*nu + 5284823040.0*r**5 + 32768.0*r**4*(362880.0*nu + 517920.0) + 32768.0*r**3*(725760.0*nu + 1035840.0) + 32768.0*r**2*(1128960.0*nu + 2071680.0) + 32768.0*r*(-725760.0*nu**2 - 1791904.94655889*nu + 4143360.0) + 440653578240.0)*sp.log(r) + 13212057600.0*r**5)
    xi_prmr = -2*sp.sqrt(Dnons)*ap**2*u*u_prmr*(Anons + ap**2*u**2)/(ap**2*u**2 + 1)**2 + sp.sqrt(Dnons)*(Anons_prmr + 2*ap**2*u*u_prmr)/(ap**2*u**2 + 1) + Dnons_prmr*(Anons + ap**2*u**2)/(sp.sqrt(Dnons)*(2*ap**2*u**2 + 2))
    pr_prmr = -prstar*xi_prmr/xi**2
    QalignSS_prmr = pr**4*(-3*am**2*(-15*nu**2/8 + 75*nu/32 - 15/32) - 3*am*ap*delta*(45*nu/8 - 5/16) - 3*ap**2*(-5*nu**2 + 165*nu/32 + 25/32))/r**4 + pr**3*pr_prmr*(4*am**2*(-15*nu**2/8 + 75*nu/32 - 15/32) + 4*am*ap*delta*(45*nu/8 - 5/16) + 4*ap**2*(-5*nu**2 + 165*nu/32 + 25/32))/r**3
    BnpalignSS_prmr = -(3*am**2*(3*nu/4 - 3/16) - 63*am*ap*delta/8 + 3*ap**2*(3*nu + 45/16))/r**4 + (-4*am**2*(nu**2/16 + 115*nu/64 - 37/64) - 4*am*ap*delta*(13*nu/16 + 449/32) - 4*ap**2*(-1171*nu/64 - 861/64))/r**5
    AalignSS_prmr = (-4*am**2*(nu/2 + 1/8) + 5*am*ap*delta - 9*ap**2/2)/r**5 - (5*am**2*(21*nu**2/16 - 81*nu/64 - 9/64) + 5*am*ap*delta*(117/32 - 39*nu/16) + 5*ap**2*(-175*nu/64 - 225/64))/r**6
    Qalign_prmr = QalignSS_prmr + Qnos_prmr
    Balignnp_prmr = Anons*Dnons_prmr + Anons_prmr*Dnons + BnpalignSS_prmr + 2*ap**2*u*u_prmr
    Bkerreqnp_prmr = (-1 - 2/r)*(2*ap**2/r**2 - 2*r)/(ap**2*(1 + 2/r) + r**2)**2 + 2/(r**2*(ap**2*(1 + 2/r) + r**2))
    Aalign_prmr = (ap**2*(2 + 4/r)/r**3 + 2*ap**2/r**4)*(AalignSS + Anons + ap**2*u**2)/(ap**2*(1 + 2/r)/r**2 + 1)**2 + (AalignSS_prmr + Anons_prmr + 2*ap**2*u*u_prmr)/(ap**2*(1 + 2/r)/r**2 + 1)
    Galigna3_prmr = pphi*(-am*ap**2*delta + ap**3)/(2*r**3)
    SOcalib_prmr = 3*ap*dSO*nu*pphi*u**2*u_prmr
    Heven_prmr = sp.sqrt(Aalign*(Bkerreqnp*ap**2*pphi**2/r**2 + Qalign + pphi**2/r**2 + pr**2*(Balignnp + 1) + 1))*(Aalign*(Balignnp_prmr*pr**2 - 2*Bkerreqnp*ap**2*pphi**2/r**3 + Bkerreqnp_prmr*ap**2*pphi**2/r**2 + Qalign_prmr - 2*pphi**2/r**3 + pr*pr_prmr*(2*Balignnp + 2))/2 + Aalign_prmr*(Bkerreqnp*ap**2*pphi**2/r**2 + Qalign + pphi**2/r**2 + pr**2*(Balignnp + 1) + 1)/2)/(Aalign*(Bkerreqnp*ap**2*pphi**2/r**2 + Qalign + pphi**2/r**2 + pr**2*(Balignnp + 1) + 1))
    Hodd_prmr = (-ap**2 - 3*r**2)*(Galigna3 + SOcalib + pphi*(am*delta*gam + ap*gap))/(ap**2*(r + 2) + r**3)**2 + (Galigna3_prmr + SOcalib_prmr + pphi*(am*delta*gam_prmr + ap*gap_prmr))/(ap**2*(r + 2) + r**3)
    Heff_prmr = Heven_prmr + Hodd_prmr
    Hreal_prmr = Heff_prmr*nu/sp.sqrt(nu*(2*Heff - 2) + 1)
    return c_codegen([Hreal_prmr/nu,Hreal_prmpphi/nu],['const double Hreal_prmr','const double Hreal_prmpphi'],verbose = False, include_braces = False)
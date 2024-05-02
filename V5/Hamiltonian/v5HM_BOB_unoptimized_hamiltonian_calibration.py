import numpy as np
def v5HM_BOB_unoptimized_hamiltonian_calibration(m1, m2, r, prstar, pphi, chi1, chi2,a6,dSO,verbose = False):
    u = 1/r
    M = m1+m2
    delta = (m1-m2)/M
    nu = m1*m2/(M**2)
    gam=(np.divide(1,4)+(pphi**2/r**2)*(np.divide(15,32)-np.divide(9,32)*nu)+(1/r)*(np.divide(11,32)*nu+np.divide(3,32))+(pphi**4/r**4)*(np.divide(75,256)*nu**2-np.divide(45,128)*nu-np.divide(105,256))+(pphi**2/r**3)*(-np.divide(613,768)*nu**2-np.divide(35,128)*nu-np.divide(59,256))+(1/r**2)*(np.divide(103,192)*nu**2-np.divide(1,32)*nu+np.divide(5,64)))
    gap=(np.divide(7,4)+(pphi**2/r**2)*(-np.divide(45,32)*nu-np.divide(15,32))+(1/r)*(np.divide(23,32)*nu-np.divide(3,32))+(pphi**4/r**4)*(np.divide(345,256)*nu**2+np.divide(75,128)*nu+np.divide(105,256))+(pphi**2/r**3)*(-np.divide(1591,768)*nu**2-np.divide(267,128)*nu+np.divide(59,256))+(1/r**2)*(np.divide(109,192)*nu**2-np.divide(177,32)*nu-np.divide(5,64)))
    am = (m1*chi1 - m2*chi2)/M
    ap = (m1*chi1 + m2*chi2)/M
    Qnos = (0.121954868780449*nu*prstar**8/r + prstar**6*(6.0*nu**3 - 5.4*nu**2 - 2.78300763695006*nu)/r**2 + prstar**4*(10.0*nu**3 - 131.0*nu**2 + 92.7110442849544*nu)/r**3) + (prstar**8*(-6.0*nu**4 + 3.42857142857143*nu**3 + 3.33842023648322*nu**2 + 1.38977750996128*nu)/r**2 + prstar**6*(-14.0*nu**4 + 188.0*nu**3 - 89.5298327361234*nu**2 - 33.9782122170436*nu)/r**3 + prstar**4*(602.318540416564*nu**3 + nu**2*(118.4*np.log(r) - 1796.13660498019) + nu*(452.542166996693 - 51.6952380952381*np.log(r) ))/r**4) + (1.48275342024365*nu*prstar**8/r**2.5 - 11.3175085791863*nu*prstar**6/r**3.5 + 147.443752990146*nu*prstar**4/r**4.5) + prstar**4*(-6.0*nu**2 + 8.0*nu)/r**2
    d5 = 0
    Dnons = r*(6730497718123.02*nu**3 + 22295347200.0*nu**2*d5 + 133772083200.0*nu**2*r**2 + 1822680546449.21*nu**2*r + 80059249540278.2*nu**2 + 22295347200.0*nu*d5*r - 193226342400.0*nu*d5 + 2589101062873.81*nu*r**2 + 10611661054566.2*nu*r - 12049908701745.2*nu + 5107745331375.71*r**2 - 326837426.241486*r*(14700.0*nu + 42911.0) - 39476764256925.6*r - (-5041721180160.0*nu**2 - 25392914995744.3*nu - 879923036160.0*r**2 - 283115520.0*r*(14700.0*nu + 42911.0) + 104186110149937.0)*np.log(r) + 5787938193408.0*np.log(r)**2 + 275059053208689.0)/(55296.0*nu*(14515200.0*nu**3 - 42636451.6032331*nu**2 - 7680.0*nu*(315.0*d5 + 890888.810272497) + 4331361844.61149*nu + 1002013764.01019) - 967680.0*r**3*(-138240.0*nu**2 - 2675575.66847905*nu - 5278341.3229329) - 9216.0*r**2*(-197773496.793534*nu**2 - 7680.0*nu*(315.0*d5 + 405152.309729121) + 2481453539.84635*nu + 5805304367.87913) + r*(5927865218923.02*nu**3 + 70778880.0*nu**2*(315.0*d5 + 2561145.80918574) - 138141470005001.0*nu**2 - 4718592.0*nu*(40950.0*d5 + 86207832.4415642) + 450172889755120.0*nu + 86618264430493.3*(1 - 0.496948781616935*nu)**2 + 188440788778196.0) + 5787938193408.0*r*np.log(r)**2 + (-1698693120.0*nu*(11592.0*nu + 69847.0) + 879923036160.0*r**3 + 283115520.0*r**2*(14700.0*nu + 42911.0) + 49152.0*r*(102574080.0*nu**2 + 409207698.136075*nu - 2119671837.36038))*np.log(r))
    Anons = 7680.0*r**4*(-5416406.59541186*nu**2 + 28.0*nu*(1920.0*a6 + 733955.307463037) + 2048.0*nu*(756.0*nu + 336.0*r + 407.0)*np.log(r) - 7.0*r*(-185763.092693281*nu**2 + 938918.400156317*nu - 245760.0) - 3440640.0)/(241555486248.807*nu**4 + 1120.0*nu**3*(-17833256.898555*r**2 - 163683964.822551*r - 1188987459.03162) + 7.0*nu**2*(-39321600.0*a6*(3.0*r + 59.0) + 745857848.115604*a6 + 1426660551.8844*r**5 - 3089250703.76879*r**4 - 6178501407.53758*r**3 + 2064783811.32587*r**2 + 122635399361.987*r + 276057889687.011) + 67645734912.0*nu**2*np.log(r)**2 + 53760.0*nu*(7680.0*a6*(r**4 + 2.0*r**3 + 4.0*r**2 + 8.0*r + 16.0) + 128.0*r*(-6852.34813868015*r**4 + 4264.6962773603*r**3 + 8529.39255472061*r**2 + 13218.7851094412*r - 33722.4297811176) + 113485.217444961*r*(-r**4 + 2.0*r**3 + 4.0*r**2 + 8.0*r + 16.0) + 148.04406601634*r*(349.0*r**4 + 1926.0*r**3 + 3852.0*r**2 + 7704.0*r + 36400.0)) + 32768.0*nu*(-1882456.23663972*nu**2 - 38842241.4769507*nu + 161280.0*r**5 + 480.0*r**4*(756.0*nu + 1079.0) + 960.0*r**3*(756.0*nu + 1079.0) + 1920.0*r**2*(588.0*nu + 1079.0) + 240.0*r*(-3024.0*nu**2 - 7466.27061066206*nu + 17264.0) + 13447680.0)*np.log(r) + 13212057600.0*r**5)
    xi = np.sqrt(Dnons) * ( Anons + ap**2*u*u ) / (1 + ap**2*u*u)
    pr = prstar/xi
    QalignSS=((pr**4)/(r**3))*(ap**2*(-5*nu*nu+nu*np.divide(165,32)+np.divide(25,32))+delta*ap*am*(nu*np.divide(45,8)-np.divide(5,16))+am**2*(-nu*nu*np.divide(15,8)+nu*np.divide(75,32)-np.divide(15,32)))
    BnpalignSS=((1/r**3)*(ap**2*(3*nu+np.divide(45,16))-delta*ap*am*np.divide(21,8)+am**2*(nu*np.divide(3,4)-np.divide(3,16)))+(1/r**4)*(ap**2*(-nu*np.divide(1171,64)-np.divide(861,64))+delta*ap*am*(nu*np.divide(13,16)+np.divide(449,32))+am**2*(nu*nu*np.divide(1,16)+nu*np.divide(115,64)-np.divide(37,64))))
    AalignSS=((1/r**4)*(ap**2*np.divide(9,8)-delta*ap*am*np.divide(5,4)+am**2*(nu*np.divide(1,2)+np.divide(1,8)))+(1/r**5)*(ap**2*(-nu*np.divide(175,64)-np.divide(225,64))+delta*ap*am*(-nu*np.divide(39,16)+np.divide(117,32))+am**2*(nu*nu*np.divide(21,16)-nu*np.divide(81,64)-np.divide(9,64))))
    Qalign = Qnos + QalignSS
    Balignnp = -1 + ap**2*u*u + Anons*Dnons + BnpalignSS
    Bkerreqnp = - (1 + 2/r)/(r**2 + ap**2*(1 + 2/r))
    Aalign = (ap**2*u*u + Anons + AalignSS )/(1 + ap**2*(1 + 2/r)/(r**2) )
    Galigna3 = pphi*(delta*am*ap**2 - ap**3)/(4*r**2)
    SOcalib = nu*dSO*ap*pphi*(u**3)
    Heven=np.sqrt(Aalign*(1+pphi*pphi/r**2+(1+Balignnp)*pr*pr+Bkerreqnp*pphi*pphi*ap**2/r**2+Qalign))
    Hodd = (pphi*(gap*ap + delta*gam*am) + SOcalib + Galigna3)/(r**3 + ap**2*(r + 2))
    Heff = Hodd + Heven
    Hreal = np.sqrt(1 + 2*nu*(Heff - 1))
    if not verbose:
        return Hreal,xi
    else:
        return Hreal,xi,Aalign,Balignnp,Bkerreqnp,Qalign,Heven,Hodd,QalignSS,Qnos,Galigna3,gam,gap,SOcalib,u,nu,ap,am,r,prstar,pphi,chi1,chi2,m1,m2
import numpy as np
def v5HM_unoptimized_IC_cons(m1,m2,r,pphi,chi1,chi2):
    prstar = 0
    u = 1/r
    
    M = m1+m2
    
    delta = (m1-m2)/M
    
    eta = m1*m2/(M**2)
    
    a6 = 41.7877875-3021.93382*eta+33414.4394*eta**2-169019.14*eta**3+329523.262*eta**4
    
    gam = (np.divide(1,4)+(pphi**2/r**2)*(np.divide(15,32)-np.divide(9,32)*eta)+(1/r)*(np.divide(11,32)*eta+np.divide(3,32))+(pphi**4/r**4)*(np.divide(75,256)*eta**2-np.divide(45,128)*eta-np.divide(105,256))+(pphi**2/r**3)*(-np.divide(613,768)*eta**2-np.divide(35,128)*eta-np.divide(59,256))+(1/r**2)*(np.divide(103,192)*eta**2-np.divide(1,32)*eta+np.divide(5,64)))
    
    gap = (np.divide(7,4)+(pphi**2/r**2)*(-np.divide(45,32)*eta-np.divide(15,32))+(1/r)*(np.divide(23,32)*eta-np.divide(3,32))+(pphi**4/r**4)*(np.divide(345,256)*eta**2+np.divide(75,128)*eta+np.divide(105,256))+(pphi**2/r**3)*(-np.divide(1591,768)*eta**2-np.divide(267,128)*eta+np.divide(59,256))+(1/r**2)*(np.divide(109,192)*eta**2-np.divide(177,32)*eta-np.divide(5,64)))
    
    am = (m1*chi1-m2*chi2)/M
    
    ap = (m1*chi1+m2*chi2)/M
    
    Qnosterm3 = (prstar**8)*(u*eta*(-np.divide(35772,175)+np.divide(21668992,45)*np.log(2)+np.divide(6591861,350)*np.log(3)-np.divide(27734375,126)*np.log(5))+u**2*(-6*eta**4+np.divide(24,7)*eta**3+eta**2*(np.divide(870976,525)+np.divide(703189497728,33075)*np.log(2)+np.divide(332067403089,39200)*np.log(3)-np.divide(13841287201,4320)*np.log(7)-np.divide(468490234375,42336)*np.log(5))+eta*(-np.divide(2222547,2450)-np.divide(1347019456,525)*np.log(2)+np.divide(278690984375,169344)*np.log(5)+np.divide(13841287201,17280)*np.log(7)-np.divide(346536085761,156800)*np.log(3)))+np.divide(5994461,12700800)*np.pi*eta*u**(np.divide(5,2)))
    
    Qnosterm2 = (prstar**6)*(u**2*(6*eta**3-np.divide(27,5)*eta**2+eta*(-np.divide(827,3)-np.divide(2358912,25)*np.log(2)+np.divide(1399437,50)*np.log(3)+np.divide(390625,18)*np.log(5)))+u**3*(-14*eta**4+188*eta**3+eta**2*(np.divide(154229,75)-np.divide(4998308864,1575)*np.log(2)+np.divide(26171875,18)*np.log(5)-np.divide(45409167,350)*np.log(3))+eta*(-np.divide(860317,1050)+np.divide(305146624,945)*np.log(2)+np.divide(35643726,175)*np.log(3)-np.divide(52468750,189)*np.log(5)))-np.divide(2723471,756000)*np.pi*eta*u**(np.divide(7,2)))
    
    Qnosterm1 = (prstar**4)*(2*(4-3*eta)*eta*u**2+u**3*(10*eta**3-131*eta**2+eta*(-np.divide(4348,15)+np.divide(496256,45)*np.log(2)-np.divide(33048,5)*np.log(3)))+u**4*((792-np.divide(615,32)*np.pi**2)*eta**3+eta**2*(-np.divide(592,5)*np.log(u)+np.divide(31633,512)*np.pi**2-np.divide(1184,5)*np.euler_gamma+np.divide(45683,105)+np.divide(33693536,105)*np.log(2)-np.divide(6396489,70)*np.log(3)-np.divide(9765625,126)*np.log(5))+eta*(np.divide(5428,105)*np.log(u)+np.divide(1249177,1050)-np.divide(93031,1536)*np.pi**2+np.divide(10856,105)*np.euler_gamma-np.divide(4396376,105)*np.log(2)+np.divide(9765625,504)*np.log(5)-np.divide(601911,280)*np.log(3)))+np.divide(88703,1890)*np.pi*eta*u**(np.divide(9,2)))
    
    Qnos = Qnosterm1+Qnosterm2+Qnosterm3
    
    d5 = 0
    
    Dnons = r*(6730497718123.02*eta**3+22295347200.0*eta**2*d5+133772083200.0*eta**2*r**2+1822680546449.21*eta**2*r+80059249540278.2*eta**2+22295347200.0*eta*d5*r-193226342400.0*eta*d5+2589101062873.81*eta*r**2+10611661054566.2*eta*r-12049908701745.2*eta+5107745331375.71*r**2-326837426.241486*r*(14700.0*eta+42911.0)-39476764256925.6*r-(-5041721180160.0*eta**2-25392914995744.3*eta-879923036160.0*r**2-283115520.0*r*(14700.0*eta+42911.0)+104186110149937.0)*np.log(r)+5787938193408.0*np.log(r)**2+275059053208689.0)/(55296.0*eta*(14515200.0*eta**3-42636451.6032331*eta**2-7680.0*eta*(315.0*d5+890888.810272497)+4331361844.61149*eta+1002013764.01019)-967680.0*r**3*(-138240.0*eta**2-2675575.66847905*eta-5278341.3229329)-9216.0*r**2*(-197773496.793534*eta**2-7680.0*eta*(315.0*d5+405152.309729121)+2481453539.84635*eta+5805304367.87913)+r*(5927865218923.02*eta**3+70778880.0*eta**2*(315.0*d5+2561145.80918574)-138141470005001.0*eta**2-4718592.0*eta*(40950.0*d5+86207832.4415642)+450172889755120.0*eta+86618264430493.3*(1-0.496948781616935*eta)**2+188440788778196.0)+5787938193408.0*r*np.log(r)**2+(-1698693120.0*eta*(11592.0*eta+69847.0)+879923036160.0*r**3+283115520.0*r**2*(14700.0*eta+42911.0)+49152.0*r*(102574080.0*eta**2+409207698.136075*eta-2119671837.36038))*np.log(r))
    
    Anons = 7680.0*r**4*(-5416406.59541186*eta**2+28.0*eta*(1920.0*a6+733955.307463037)+2048.0*eta*(756.0*eta+336.0*r+407.0)*np.log(r)-7.0*r*(-185763.092693281*eta**2+938918.400156317*eta-245760.0)-3440640.0)/(241555486248.807*eta**4+1120.0*eta**3*(-17833256.898555*r**2-163683964.822551*r-1188987459.03162)+7.0*eta**2*(-39321600.0*a6*(3.0*r+59.0)+745857848.115604*a6+1426660551.8844*r**5-3089250703.76879*r**4-6178501407.53758*r**3+2064783811.32587*r**2+122635399361.987*r+276057889687.011)+67645734912.0*eta**2*np.log(r)**2+53760.0*eta*(7680.0*a6*(r**4+2.0*r**3+4.0*r**2+8.0*r+16.0)+128.0*r*(-6852.34813868015*r**4+4264.6962773603*r**3+8529.39255472061*r**2+13218.7851094412*r-33722.4297811176)+113485.217444961*r*(-r**4+2.0*r**3+4.0*r**2+8.0*r+16.0)+148.04406601634*r*(349.0*r**4+1926.0*r**3+3852.0*r**2+7704.0*r+36400.0))+32768.0*eta*(-1882456.23663972*eta**2-38842241.4769507*eta+161280.0*r**5+480.0*r**4*(756.0*eta+1079.0)+960.0*r**3*(756.0*eta+1079.0)+1920.0*r**2*(588.0*eta+1079.0)+240.0*r*(-3024.0*eta**2-7466.27061066206*eta+17264.0)+13447680.0)*np.log(r)+13212057600.0*r**5)
    
    xi = np.sqrt(Dnons)*(Anons+ap**2*u*u)/(1+ap**2*u*u)
    
    pr = prstar/xi
    
    QalignSS = ((pr**4)/(r**3))*(ap**2*(-5*eta*eta+eta*np.divide(165,32)+np.divide(25,32))+delta*ap*am*(eta*np.divide(45,8)-np.divide(5,16))+am**2*(-eta*eta*np.divide(15,8)+eta*np.divide(75,32)-np.divide(15,32)))
    
    BnpalignSS = ((1/r**3)*(ap**2*(3*eta+np.divide(45,16))-delta*ap*am*np.divide(21,8)+am**2*(eta*np.divide(3,4)-np.divide(3,16)))+(1/r**4)*(ap**2*(-eta*np.divide(1171,64)-np.divide(861,64))+delta*ap*am*(eta*np.divide(13,16)+np.divide(449,32))+am**2*(eta*eta*np.divide(1,16)+eta*np.divide(115,64)-np.divide(37,64))))
    
    AalignSS = ((1/r**4)*(ap**2*np.divide(9,8)-delta*ap*am*np.divide(5,4)+am**2*(eta*np.divide(1,2)+np.divide(1,8)))+(1/r**5)*(ap**2*(-eta*np.divide(175,64)-np.divide(225,64))+delta*ap*am*(-eta*np.divide(39,16)+np.divide(117,32))+am**2*(eta*eta*np.divide(21,16)-eta*np.divide(81,64)-np.divide(9,64))))
    
    Qalign = Qnos+QalignSS
    
    Balignnp = -1+ap**2*u*u+Anons*Dnons+BnpalignSS
    
    Bkerreqnp = -(1+2/r)/(r**2+ap**2*(1+2/r))
    
    Aalign = (ap**2*u*u+Anons+AalignSS)/(1+ap**2*(1+2/r)/(r**2))
    
    Galigna3 = pphi*(delta*am*ap**2-ap**3)/(4*r**2)
    
    dSO = -7.71251231383957*am**3-17.2294679794015*am**2*ap-238.430383378296*am**2*eta+69.5461667822545*am**2-10.5225438990315*am*ap**2+362.767393298729*am*ap*eta-85.8036338010274*am*ap-1254.66845939312*am*eta**2+472.431937787377*am*eta-39.742317057316*am-7.58458103577458*ap**3-42.7601129678844*ap**2*eta+18.1783435552183*ap**2-201.905934468847*ap*eta**2-90.5790079104259*ap*eta+49.6299175121658*ap+478.546231305475*eta**3+679.521769948995*eta**2-177.334831768076*eta-37.6897780220529
    
    SOcalib = eta*dSO*ap*pphi*(u**3)
    
    Heven = np.sqrt(Aalign*(1+pphi*pphi/r**2+(1+Balignnp)*pr*pr+Bkerreqnp*pphi*pphi*ap**2/r**2+Qalign))
    
    Hodd = (pphi*(gap*ap+delta*gam*am)+SOcalib+Galigna3)/(r**3+ap**2*(r+2))
    
    Heff = Hodd+Heven
    
    Hreal = np.sqrt(1+2*eta*(Heff-1))
    
    
    gam_prmpphi_preq0 = pphi**3*(75*eta**2/64 - 45*eta/32 - 105/64)/r**4 + pphi*(15/16 - 9*eta/16)/r**2 + pphi*(-613*eta**2/384 - 35*eta/64 - 59/128)/r**3
    gap_prmpphi_preq0 = pphi**3*(345*eta**2/64 + 75*eta/32 + 105/64)/r**4 + pphi*(-45*eta/16 - 15/16)/r**2 + pphi*(-1591*eta**2/384 - 267*eta/64 + 59/128)/r**3
    Galigna3_prmpphi_preq0 = (am*ap**2*delta - ap**3)/(4*r**2)
    SOcalib_prmpphi_preq0 = ap*dSO*eta*u**3
    Heven_prmpphi_preq0 = np.sqrt(Aalign*(Bkerreqnp*ap**2*pphi**2/r**2 + Qalign + pphi**2/r**2 + pr**2*(Balignnp + 1) + 1))*(2*Bkerreqnp*ap**2*pphi/r**2 + 2*pphi/r**2)/(2*(Bkerreqnp*ap**2*pphi**2/r**2 + Qalign + pphi**2/r**2 + pr**2*(Balignnp + 1) + 1))
    Hodd_prmpphi_preq0 = (Galigna3_prmpphi_preq0 + SOcalib_prmpphi_preq0 + am*delta*gam + ap*gap + pphi*(am*delta*gam_prmpphi_preq0 + ap*gap_prmpphi_preq0))/(ap**2*(r + 2) + r**3)
    Heff_prmpphi_preq0 = Heven_prmpphi_preq0 + Hodd_prmpphi_preq0
    Hreal_prmpphi_preq0 = Heff_prmpphi_preq0*eta/np.sqrt(eta*(2*Heff - 2) + 1)
    gam_prmr_preq0 = pphi**3*(75*eta**2/64 - 45*eta/32 - 105/64)/r**4 + pphi*(15/16 - 9*eta/16)/r**2 + pphi*(-613*eta**2/384 - 35*eta/64 - 59/128)/r**3
    gap_prmr_preq0 = pphi**3*(345*eta**2/64 + 75*eta/32 + 105/64)/r**4 + pphi*(-45*eta/16 - 15/16)/r**2 + pphi*(-1591*eta**2/384 - 267*eta/64 + 59/128)/r**3
    Galigna3_prmr_preq0 = (am*ap**2*delta - ap**3)/(4*r**2)
    SOcalib_prmr_preq0 = ap*dSO*eta*u**3
    Heven_prmr_preq0 = np.sqrt(Aalign*(Bkerreqnp*ap**2*pphi**2/r**2 + Qalign + pphi**2/r**2 + pr**2*(Balignnp + 1) + 1))*(2*Bkerreqnp*ap**2*pphi/r**2 + 2*pphi/r**2)/(2*(Bkerreqnp*ap**2*pphi**2/r**2 + Qalign + pphi**2/r**2 + pr**2*(Balignnp + 1) + 1))
    Hodd_prmr_preq0 = (Galigna3_prmr_preq0 + SOcalib_prmr_preq0 + am*delta*gam + ap*gap + pphi*(am*delta*gam_prmr_preq0 + ap*gap_prmr_preq0))/(ap**2*(r + 2) + r**3)
    Heff_prmr_preq0 = Heven_prmr_preq0 + Hodd_prmr_preq0
    Hreal_prmr_preq0 = Heff_prmr_preq0*eta/np.sqrt(eta*(2*Heff - 2) + 1)
    return np.array([Hreal_prmr_preq0, Hreal_prmpphi_preq0])
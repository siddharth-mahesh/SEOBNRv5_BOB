import numpy as np
from scipy.special import gamma, factorial, factorial2
def v5HM_unoptimized_flux(m1, m2, r, phi, prstar, pphi, chi1, chi2,Hreal,Omega,Omega_circ,verbose = False):
    M = m1+m2
    eta = m1*m2/(M**2)
    TINYDOUBLE = 1e-100
    delta = (m1 - m2)/M
    min_delta_1em14 = np.divide(1,2)*(delta - 1e-14 + 0 - np.absolute(delta - 1e-14 - 0))
    max_delta_1em14 = np.divide(1,2)*(delta - 1e-14 + TINYDOUBLE + 0 + np.absolute(delta - 1e-14 + TINYDOUBLE - 0))
    noneqcond = max_delta_1em14/(delta - 1e-14+TINYDOUBLE)
    eqcond = min_delta_1em14/(delta - 1e-14-TINYDOUBLE)
    deltainvertible = delta*noneqcond + 1*eqcond
    deltainv = 1/deltainvertible
    chi_S = np.divide(1,2)*(chi1 + chi2)
    chi_A = np.divide(1,2)*(chi1 - chi2)
    vomega = Omega**(np.divide(1,3))
    vphi = Omega*(Omega_circ**(-np.divide(2,3)))
    Heff = (Hreal**2 - 1)/(2*eta) + 1
    eulerlog5 = np.euler_gamma + np.log(2*5*vomega)
    eulerlog4 = np.euler_gamma + np.log(2*4*vomega)
    eulerlog3 = np.euler_gamma + np.log(2*3*vomega)
    eulerlog2v2 = np.euler_gamma + np.log(2*2*vomega**2)
    eulerlog2 = np.euler_gamma + np.log(2*2*vomega)
    eulerlog1v2 = np.euler_gamma + np.log(2*1*vomega**2)
    eulerlog1 = np.euler_gamma + np.log(2*1*vomega)
    khat1 = 1*Omega*Hreal
    khat2 = 2*Omega*Hreal
    khat3 = 3*Omega*Hreal
    khat4 = 4*Omega*Hreal
    khat5 = 5*Omega*Hreal
    khat6 = 6*Omega*Hreal
    khat7 = 7*Omega*Hreal
    khat8 = 8*Omega*Hreal
    delta55=((np.divide(96875,131250)+np.divide(857528,131250)*eta)*(Omega*Hreal)/(1-2*eta)+np.divide(3865,429)*np.pi*(Omega*Hreal)**2+((np.divide(-7686949127,31783752)+np.divide(954500400,31783752)*np.pi**2)*(Omega*Hreal)**3))
    delta43=(((np.divide(4961,810)*eta+np.divide(3,5))*(Omega*Hreal))/(1-2*eta)+np.divide(1571,385)*np.pi*(Omega*Hreal)**2)
    delta44=((np.divide(112,120)+np.divide(219,120)*eta)*(Omega*Hreal)+np.divide(25136,3465)*(Omega*Hreal)**2+(np.divide(201088,10395)*np.pi**2-np.divide(55144,375))*(Omega*Hreal)**3)
    delta32=(((np.divide(11,5)*eta+np.divide(2,3))*(Hreal*Omega))/(1-3*eta)+np.divide(52,21)*np.pi*(Hreal*Omega)**2+((np.divide(208,63)*np.pi**2)-np.divide(9112,405))*(Hreal*Omega)**3)
    delta33=(np.divide(13,10)*(Hreal*Omega)+np.divide(39,7)*(Hreal*Omega)**2+(-np.divide(227827,3000)+np.divide(78,7)*np.pi**2)*(Hreal*Omega)**3-np.divide(80897,2430)*eta*vomega**5)
    delta21=(np.divide(2,3)*Omega*Hreal+(np.divide(107,105)*np.pi)*(Omega*Hreal)**2+((np.divide(214,315)*np.pi**2)-(np.divide(272,81)))*(Omega*Hreal)**3-(np.divide(25,2)*eta*vomega**5))
    delta22=(np.divide(7,3)*Omega*Hreal+(Omega*Hreal)**2*(((np.divide(8,3)*eta-np.divide(4,3))*chi_S)-(np.divide(4,3)*delta*chi_A)+(np.divide(428,105)*np.pi))+(Omega*Hreal)**3*((np.divide(1712,315)*np.pi**2)-(np.divide(2203,81)))-24*eta*vomega**5)
    fspin55_limit=(((np.divide(-70,3)*eta/((-1+2*eta))+np.divide(110,3)*eta**2/((-1+2*eta))+np.divide(10,3)/((-1+2*eta)))*(chi_A))*vomega**3+((-5.0/(-1.0+2.0*eta)+30.0*eta/(-1.0+2.0*eta)-40.0*eta**2/(-1.0+2.0*eta))*chi_A*chi_S)*vomega**4)
    fspin55=(((-70*eta/(3*(-1+2*eta))+110*eta**2/(3*(-1+2*eta))+10/(3*(-1+2*eta)))*(chi_A*deltainv)+(10/(3*(-1+2*eta))-10*eta/(-1+2*eta)+10*eta**2/(-1+2*eta))*chi_S)*vomega**3+(np.divide(5,2)*chi_S**2+(-5.0/(-1.0+2.0*eta)+30.0*eta/(-1.0+2.0*eta)-40.0*eta**2/(-1.0+2.0*eta))*chi_A*chi_S*deltainv+(-5.0/(2.0*(-1.0+2.0*eta))+15.0*eta/(-1.0+2.0*eta)-20.0*eta**2/(-1.0+2.0*eta))*chi_A**2)*vomega**4)
    fspin41_limit=(np.divide(5,2)*eta*vomega/(1-2*eta)*(-chi_A))
    fspin41=(np.divide(5,2)*eta*vomega/(1-2*eta)*(chi_S-chi_A*deltainv))
    fspin43_limit=(vomega/(1-2*eta)*(-np.divide(5,2)*eta*chi_A)+vomega**3/(1-2*eta)*((np.divide(887,44)*eta-np.divide(3143,132)*eta**2)*chi_A)+vomega**4/(1-2*eta)*((np.divide(137,6)*eta**2-18*eta+3)*chi_A*chi_S))
    fspin43=(vomega/(1-2*eta)*(np.divide(5,2)*eta*chi_S-np.divide(5,2)*eta*chi_A*deltainv)+vomega**3/(1-2*eta)*((np.divide(887,44)*eta-np.divide(3143,132)*eta**2)*chi_A*deltainv+(-np.divide(529,132)*eta**2-np.divide(667,44)*eta)*chi_S)+vomega**4/(1-2*eta)*((12*eta**2-np.divide(37,3)*eta+np.divide(3,2))*chi_A**2+(np.divide(137,6)*eta**2-18*eta+3)*chi_A*chi_S*deltainv+(np.divide(35,6)*eta**2+np.divide(1,3)*eta+np.divide(3,2))*chi_S**2))
    fspin31_limit = vomega**3 * (-chi_A * 5.0 / 8.0)
    fspin31=vomega**3*(chi_A*deltainv*(-4.0+11.0*eta)+chi_S*(-4.0+13.0*eta))/(2.0)
    fspin33_limit=vomega**3*((np.divide(19,2)*eta-2)*chi_A)+vomega**4*((3-12*eta)*(chi_A*chi_S))+vomega**5*((np.divide(407,30)*eta**2-np.divide(593,60)*eta+np.divide(2,3))*chi_A)+vomega**6*((44*eta**2-eta-np.divide(7,2))*(chi_A*chi_S))+1j*(Omega*Hreal)**2*(np.divide(7339,540)*eta-np.divide(81,20))*chi_A
    fspin33amp_limit=vomega**3*((np.divide(19,2)*eta-2)*chi_A)+vomega**4*((3-12*eta)*(chi_A*chi_S))+vomega**5*((np.divide(407,30)*eta**2-np.divide(593,60)*eta+np.divide(2,3))*chi_A)+vomega**6*((44*eta**2-eta-np.divide(7,2))*(chi_A*chi_S))
    fspin33=(vomega**3*(((np.divide(19,2)*eta-2)*chi_A*deltainv)+((np.divide(5,2)*eta-2)*chi_S))+vomega**4*((np.divide(3,2)-6*eta)*chi_A**2+(3-12*eta)*(chi_A*chi_S*deltainv)+np.divide(3,2)*chi_S**2)+vomega**5*(((np.divide(407,30)*eta**2-np.divide(593,60)*eta+np.divide(2,3))*chi_A*deltainv)+((np.divide(241,30)*eta**2+np.divide(11,20)*eta+np.divide(2,3))*chi_S))+vomega**6*((-12*eta**2+np.divide(11,2)*eta-np.divide(7,4))*chi_A**2+(44*eta**2-eta-np.divide(7,2))*(chi_A*chi_S*deltainv)+(6*eta**2-np.divide(27,2)*eta-np.divide(7,4))*chi_S**2)+1j*((Omega*Hreal)**2*(np.divide(7339,540)*eta-np.divide(81,20))*chi_A*deltainv+(np.divide(593,108)*eta-np.divide(81,20))*chi_S))
    fspin33amp=(vomega**3*(((np.divide(19,2)*eta-2)*chi_A*deltainv)+((np.divide(5,2)*eta-2)*chi_S))+vomega**4*((np.divide(3,2)-6*eta)*chi_A**2+(3-12*eta)*(chi_A*chi_S*deltainv)+np.divide(3,2)*chi_S**2)+vomega**5*(((np.divide(407,30)*eta**2-np.divide(593,60)*eta+np.divide(2,3))*chi_A*deltainv)+((np.divide(241,30)*eta**2+np.divide(11,20)*eta+np.divide(2,3))*chi_S))+vomega**6*((-12*eta**2+np.divide(11,2)*eta-np.divide(7,4))*chi_A**2+(44*eta**2-eta-np.divide(7,2))*(chi_A*chi_S*deltainv)+(6*eta**2-np.divide(27,2)*eta-np.divide(7,4))*chi_S**2))
    fspin21_limit=(-np.divide(3,2)*vomega*(chi_A)+vomega**3*((np.divide(131,84)*eta+np.divide(61,12))*(chi_A))+vomega**4*(+(np.divide(21,2)*eta-6)*(chi_A*chi_S))+vomega**5*((np.divide(-703,112)*eta**2+np.divide(8797,1008)*eta-np.divide(81,16))*(chi_A)+(np.divide(3,4)-3*eta)*(chi_A**3)+(np.divide(9,4)-6*eta)*(chi_A*chi_S**2))+vomega**6*((np.divide(9487,504)*eta**2-np.divide(1636,21)*eta+np.divide(4163,126))*(chi_A*chi_S)))
    fspin21=(-np.divide(3,2)*vomega*(chi_A*deltainv+chi_S)+vomega**3*((np.divide(131,84)*eta+np.divide(61,12))*(chi_A*deltainv)+(np.divide(79,84)*eta+np.divide(61,12))*chi_S)+vomega**4*((-2*eta-3)*chi_A**2+(np.divide(21,2)*eta-6)*(chi_A*chi_S*deltainv)+(np.divide(1,2)*eta-3)*chi_S**2)+vomega**5*((np.divide(-703,112)*eta**2+np.divide(8797,1008)*eta-np.divide(81,16))*(chi_A*deltainv)+(np.divide(613,1008)*eta**2+np.divide(1709,1008)*eta-np.divide(81,16))*chi_S+(np.divide(3,4)-3*eta)*(chi_A**3*deltainv)+(np.divide(9,4)-6*eta)*(chi_A*chi_S**2*deltainv)+(np.divide(9,4)-3*eta)*(chi_A**2*chi_S)+np.divide(3,4)*chi_S**3)+vomega**6*((np.divide(5,7)*eta**2-np.divide(9287,1008)*eta+np.divide(4163,252))*chi_A**2+(np.divide(139,72)*eta**2-np.divide(2633,1008)*eta+np.divide(4163,252))*chi_S**2+(np.divide(9487,504)*eta**2-np.divide(1636,21)*eta+np.divide(4163,126))*(chi_A*chi_S*deltainv)))
    rhoNS81=1+((20022-126451*eta+236922*eta**2-138430*eta**3+21640*eta**4)/(18240*(-1+6*eta-10*eta**2+4*eta**3)))*vomega**2
    rho82=1+((2462-17598*eta+37119*eta**2-22845*eta**3+3063*eta**4)/(2736*(-1+7*eta-14*eta**2+7*eta**3)))*vomega**2
    rhoNS83=1+((20598-131059*eta+249018*eta**2-149950*eta**3+24520*eta**4)/(18240*(-1+6*eta-10*eta**2+4*eta**3)))*vomega**2
    rho84=1+((2666-19434*eta+42627*eta**2-28965*eta**3+4899*eta**4)/(2736*(-1+7*eta-14*eta**2+7*eta**3)))*vomega**2
    rhoNS85=1+((4350-28055*eta+54642*eta**2-34598*eta**3+6056*eta**4)/(3648*(-1+6*eta-10*eta**2+4*eta**3)))*vomega**2
    rho86=1+((1002-7498*eta+17269*eta**2-13055*eta**3+2653*eta**4)/(912*(-1+7*eta-14*eta**2+7*eta**3)))*vomega**2
    rhoNS87=1+((23478-154099*eta+309498*eta**2-207550*eta**3+38920*eta**4)/(18240*(-1+6*eta-10*eta**2+4*eta**3)))*vomega**2
    rho88=1+((3482-26778*eta+64659*eta**2-53445*eta**3+12243*eta**4)/(2736*(-1+7*eta-14*eta**2+7*eta**3)))*vomega**2
    rhoNS71=1+((228*eta**3-2083*eta**2+2518*eta-618)/(714*(3*eta**2-4*eta+1)))*vomega**2
    rho72=1+((32760*eta**4-190239*eta**3+273924*eta**2-123489*eta+16832)/(14994*(7*eta**3-14*eta**2+7*eta-1)))*vomega**2
    rhoNS73=1+((420*eta**3-2563*eta**2+2806*eta-666)/(714*(3*eta**2-4*eta+1)))*vomega**2
    rho74=1+((41076*eta**4-217959*eta**3+298872*eta**2-131805*eta+17756)/(14994*(7*eta**3-14*eta**2+7*eta-1)))*vomega**2
    rhoNS75=1+((804*eta**3-3523*eta**2+3382*eta-762)/(714*(3*eta**2-4*eta+1)))*vomega**2
    rho76=1+((6104*eta**4-29351*eta**3+37828*eta**2-16185*eta+2144)/(1666*(7*eta**3-14*eta**2+7*eta-1)))*vomega**2
    rhoNS77=1+((1380*eta**3-4963*eta**2+4246*eta-906)/(714*(3*eta**2-4*eta+1)))*vomega**2
    rhoNS61=1+((124*eta**3-670*eta**2+694*eta-161)/(144*(3*eta**2-4*eta+1)))*vomega**2
    rho62=1+((49*eta**3-413*eta**2+378*eta-74)/(84*(5*eta**2-5*eta+1)))*vomega**2-np.divide(817991,3298680)*vomega**4
    rhoNS63=1+((156*eta**3-750*eta**2+742*eta-169)/(144*(3*eta**2-4*eta+1)))*vomega**2
    rho64=1+vomega**2*((133*eta**3-581*eta**2+462*eta-86)/(84*(5*eta**2-5*eta+1)))-np.divide(476887,659736)*vomega**4
    rhoNS65=1+((220*eta**3-910*eta**2+838*eta-185)/(144*(3*eta**2-4*eta+1)))*vomega**2
    rho66=1+((273*eta**3-861*eta**2+602*eta-106)/(84*(5*eta**2-5*eta+1)))*vomega**2-np.divide(1025435,659736)*vomega**4
    rhoNS51=1+((8*eta**2-626*eta+319)/(390*(2*eta-1)))*vomega**2-np.divide(31877,304200)*vomega**4
    rho52=1+((21980*eta**3-104930*eta**2+84679*eta-15828)/(13650*(5*eta**2-5*eta+1)))*vomega**2-np.divide(7187914,15526875)*vomega**4
    rhoNS53=1+((176*eta**2-850*eta+375)/(390*(2*eta-1)))*vomega**2-np.divide(410833,709800)*vomega**4
    rho54=1+((33320*eta**3-127610*eta**2+96019*eta-17448)/(13650*(5*eta**2-5*eta+1)))*vomega**2-np.divide(16213384,15526875)*vomega**4
    rhoNS55=1+vomega**2*(np.divide(487,390)/(-1+2*eta)-np.divide(649,195)*eta/(-1+2*eta)+np.divide(256,195)*eta**2/(-1+2*eta))-np.divide(3353747,2129400)*vomega**4+(np.divide(190606537999247,11957879934000)-np.divide(1546,429)*eulerlog5)*vomega**6+(-np.divide(1213641959949291437,118143853747920000)+np.divide(376451,83655)*eulerlog5)*vomega**8+(-np.divide(150082616449726042201261,4837990810977324000000)+np.divide(2592446431,456756300)*eulerlog5)*vomega**10+eta*(-2.61*vomega**4+1.25*vomega**6+-35.7*vomega**8)
    rhoNS41=1+((288*eta**2-1385*eta+602)/(528*(2*eta-1)))*vomega**2-(np.divide(7775491,21141120))*vomega**4+(np.divide(1227423222031,1758095539200)-np.divide(1571,6930)*eulerlog1)*vomega**6
    rho42=(1+vomega**2*((1146.0-3530.0*eta+285.0*eta**2)/(1320.0*(-1+3*eta)))+vomega**3*((chi_A*(10.0-21.0*eta)*delta+chi_S*(10.0-59.0*eta+78.0*eta**2))/(15.0*(-1+3*eta)))+vomega**4*((-114859044.0+295834536.0*eta+1204388696.0*eta**2-3047981160.0*eta**3-379526805.0*eta**4)/(317116800.0*(-1+3*eta)**2))+vomega**6*(848238724511.0/219761942400.-(3142.0/3465.0)*eulerlog2))
    rhoNS43=1+vomega**2/(1-2*eta)*(-np.divide(10,11)*eta**2+np.divide(547,176)*eta-np.divide(111,88))-np.divide(6894273,7047040)*vomega**4+vomega**6*(np.divide(1664224207351,195343948800)-np.divide(1571,770)*eulerlog3)+vomega**8*(-np.divide(2465107182496333,460490801971200)+np.divide(174381,67760)*eulerlog3)+eta*(-0.654*vomega**4+-3.69*vomega**6+18.5*vomega**8)
    rho44=(1+vomega**2*((1614.0-5870.0*eta+2625.0*eta**2)/(1320.0*(-1+3*eta)))+vomega**3*(chi_A*(10.0-39.0*eta)*delta+chi_S*(10.0-41.0*eta+42.0*eta**2))/(15.0*(-1+3*eta))+vomega**4*((-511573572.0+2338945704.0*eta-313857376.0*eta**2-6733146000.0*eta**3+1252563795.0*eta**4)/(317116800.0*(-1+3*eta)**2)+chi_S**2/2.0+delta*chi_S*chi_A+delta**2*chi_A**2/2.0)+vomega**5*(chi_A*delta*(-8280.0+42716.0*eta-57990.0*eta**2+8955*eta**3)/(6600.0*(-1+3*eta)**2)+chi_S*(-8280.0+66284.0*eta-176418.0*eta**2+128085.0*eta**3+88650*eta**2*eta**2)/(6600.0*(-1+3*eta)**2))+vomega**6*(np.divide(16600939332793,1098809712000)-np.divide(12568.0,3465.0)*eulerlog4)+vomega**8*(-np.divide(172066910136202271,19426955708160000)+np.divide(845198.0,190575.0)*eulerlog4)+vomega**10*(-np.divide(17154485653213713419357,568432724020761600000)+np.divide(22324502267,3815311500)*eulerlog4)+eta*(-3.56*vomega**6+15.6*vomega**8+-216*vomega**10))
    rhoNS31=(1-vomega**2*(np.divide(2,9)*eta+np.divide(13,18))+vomega**4*(-np.divide(829,1782)*eta**2-np.divide(1685,1782)*eta+np.divide(101,7128))+vomega**6*(np.divide(11706720301,6129723600)-np.divide(26,63)*eulerlog1)+vomega**8*(np.divide(169,567)*eulerlog1+np.divide(2606097992581,4854741091200)))
    rho32=1+vomega*((4*eta*chi_S)/(3*(1-3*eta)))+vomega**2*((np.divide(-32,27)*eta**2+np.divide(223,54)*eta-np.divide(164,135))/(1-3*eta)-(16*eta**2*chi_S**2)/(9*(1-3*eta)**2))+vomega**3*((np.divide(13,9)*eta+np.divide(2,9))*(delta*chi_A)/(1-3*eta)+(np.divide(607,81)*eta**3+np.divide(503,81)*eta**2-np.divide(1478,405)*eta+np.divide(2,9))*chi_S/(1-3*eta)**2+(320*eta**3*chi_S**3)/(81*(1-3*eta)**3))+vomega**4*((np.divide(77141,40095)*eta**4-np.divide(508474,40095)*eta**3-np.divide(945121,320760)*eta**2+np.divide(1610009,320760)*eta-np.divide(180566,200475))/(1-3*eta)**2+(4*eta**2-3*eta+np.divide(1,3))*(chi_A**2)/(1-3*eta)+(np.divide(-50,27)*eta**2-np.divide(88,27)*eta+np.divide(2,3))*(delta*chi_A*chi_S)/(1-3*eta)**2+(np.divide(-2452,243)*eta**4-np.divide(1997,243)*eta**3+np.divide(1435,243)*eta**2-np.divide(43,27)*eta+np.divide(1,3))*(chi_S**2)/((1-3*eta)**3))+vomega**5*((np.divide(-1184225,96228)*eta**5-np.divide(40204523,962280)*eta**4+np.divide(101706029,962280)*eta**3-np.divide(14103833,192456)*eta**2+np.divide(20471053,962280)*eta-np.divide(2788,1215))*chi_S/(1-3*eta)**3+(np.divide(608,81)*eta**3+np.divide(736,81)*eta**2-np.divide(16,9)*eta)*(delta*chi_A*chi_S**2)/(1-3*eta)**3+(np.divide(889673,106920)*eta**3-np.divide(75737,5346)*eta**2+np.divide(376177,35640)*eta-np.divide(2788,1215))*(delta*chi_A)/(1-3*eta)**2+(np.divide(96176,2187)*eta**5+np.divide(43528,2187)*eta**4-np.divide(40232,2187)*eta**3+np.divide(376,81)*eta**2-np.divide(8,9)*eta)*(chi_S**3)/(1-3*eta)**4+(np.divide(-32,3)*eta**3+8*eta**2-np.divide(8,9)*eta)*(chi_A**2*chi_S)/((1-3*eta)**2))+vomega**6*(np.divide(5849948554,940355325)-np.divide(104,63)*eulerlog2)+vomega**8*(np.divide(17056,8505)*eulerlog2-np.divide(10607269449358,3072140846775))+vomega**10*(-np.divide(1312549797426453052,176264081083715625)+np.divide(18778864,12629925)*eulerlog2)+eta*(+.333*vomega**6-6.5*vomega**8+98*vomega**10)
    rhoNS33=1+vomega**2*(np.divide(2,3)*eta-np.divide(7,6))+vomega**4*(-np.divide(6719,3960)-np.divide(1861,990)*eta+np.divide(149,330)*eta**2)+vomega**6*(np.divide(3203101567,227026800)+(-np.divide(129509,25740)+np.divide(41,192)*np.pi**2)*eta-np.divide(274621,154440)*eta**2+np.divide(12011,46332)*eta**3-np.divide(26,7)*eulerlog3)+vomega**8*(-np.divide(57566572157,8562153600)+np.divide(13,3)*eulerlog3)+vomega**10*(-np.divide(903823148417327,30566888352000)+np.divide(87347,13860)*eulerlog3)+eta*(12*vomega**8+-215*vomega**10)
    rhoNS21=1+vomega**2*(np.divide(23,84)*eta-np.divide(59,56))+vomega**4*(np.divide(617,4704)*eta**2-np.divide(10993,14112)*eta-np.divide(47009,56448))+vomega**6*(np.divide(7613184941,2607897600)-np.divide(107,105)*eulerlog1)+vomega**8*(-np.divide(1168617463883,911303737344)+np.divide(6313,5880)*eulerlog1)+vomega**10*(-np.divide(63735873771463,16569158860800)+np.divide(5029963,5927040)*eulerlog1)+eta*(1.65*vomega**6+26.5*vomega**8+80*vomega**10)
    rho22=1+vomega**2*(np.divide(55,84)*eta-np.divide(43,42))+vomega**3*((-2.0*(chi_S+chi_A*delta-chi_S*eta))/3.0)+vomega**4*(np.divide(19583,42336)*eta**2-np.divide(33025,21168)*eta-np.divide(20555,10584)+(np.divide(1,2)-2*eta)*chi_A**2+delta*chi_A*chi_S+np.divide(1,2)*chi_S**2)+vomega**5*(delta*(-np.divide(19,42)*eta-np.divide(34,21))*chi_A+(np.divide(209,126)*eta**2+np.divide(49,18)*eta-np.divide(34,21))*chi_S)+vomega**6*(np.divide(10620745,39118464)*eta**3-np.divide(6292061,3259872)*eta**2+np.divide(41,192)*np.pi**2*eta-np.divide(48993925,9779616)*eta-np.divide(428,105)*eulerlog2+np.divide(1556919113,122245200)+delta*(np.divide(89,126)-np.divide(781,252)*eta)*chi_A*chi_S+(-np.divide(27,14)*eta**2-np.divide(457,504)*eta+np.divide(89,252))*chi_A**2+(np.divide(10,9)*eta**2-np.divide(1817,504)*eta+np.divide(89,252))*chi_S**2)+vomega**7*(delta*(np.divide(97865,63504)*eta**2+np.divide(50140,3969)*eta+np.divide(18733,15876))*chi_A+(np.divide(50803,63504)*eta**3-np.divide(245717,63504)*eta**2+np.divide(74749,5292)*eta+np.divide(18733,15876))*chi_S+delta*chi_A**3*(np.divide(1,3)-np.divide(4,3)*eta)+delta*(2*eta+1)*chi_A*chi_S**2+(np.divide(-4,1)*eta**2-np.divide(3,1)*eta+np.divide(1,1))*chi_A**2*chi_S+(eta+np.divide(1,3))*chi_S**3)+vomega**8*(np.divide(9202,2205)*eulerlog2-np.divide(387216563023,160190110080))+vomega**10*(np.divide(439877,55566)*eulerlog2-np.divide(16094530514677,533967033600))+eta*(21.2*vomega**8+-411*vomega**10)
    f81amp = (rhoNS81)**8
    f82amp = (rho82)**8
    f83amp = (rhoNS83)**8
    f84amp = (rho84)**8
    f85amp = (rhoNS85)**8
    f86amp = (rho86)**8
    f87amp = (rhoNS87)**8
    f88amp = (rho88)**8
    f71amp = (rhoNS71)**7
    f72amp = (rho72)**7
    f73amp = (rhoNS73)**7
    f74amp = (rho74)**7
    f75amp = (rhoNS75)**7
    f76amp = (rho76)**7
    f77amp = (rhoNS77)**7
    f61amp = (rhoNS61)**6
    f62amp = (rho62)**6
    f63amp = (rhoNS63)**6
    f64amp = (rho64)**6
    f65amp = (rhoNS65)**6
    f66amp = (rho66)**6
    f51amp = ((rhoNS51)**5)
    f52amp = (rho52)**5
    f53amp = ((rhoNS53)**5)
    f54amp = (rho54)**5
    f55 = ((rhoNS55)**5 + fspin55)*noneqcond + fspin55_limit*eqcond
    f55amp = f55
    f41amp = ((rhoNS41)**4 + fspin41)*noneqcond + fspin41_limit*eqcond
    f42amp = (rho42)**4
    f43 = ((rhoNS43)**4 + fspin43)*noneqcond + fspin43_limit*eqcond
    f43amp = f43
    f44 = (rho44)**4
    f44amp = f44
    f31amp = ((rhoNS31)**3 + fspin31)*noneqcond + fspin31_limit*eqcond
    f32 = (rho32)**3
    f32amp = f32
    f33 = ((rhoNS33)**3 + fspin33)*noneqcond + fspin33_limit*eqcond
    f33amp = ((rhoNS33)**3 + fspin33amp)*noneqcond + fspin33amp_limit*eqcond
    f21 = ((rhoNS21)**2 + fspin21)*noneqcond + fspin21_limit*eqcond
    f21amp = f21
    f22 = (rho22)**2
    f22amp = f22
    Y82amp = 0.32254835519288305
    Y84amp = 0.3382915688890245
    Y86amp = 0.3764161087284946
    Y88amp = 0.5154289843972844
    Y71amp = 0.31937046138540076
    Y73amp = 0.331899519333737
    Y75amp = 0.3669287245764378
    Y77amp = 0.5000395635705508
    Y62amp = 0.32569524293385776
    Y64amp = 0.3567812628539981
    Y66amp = 0.48308411358006625
    Y51amp = 0.32028164857621516
    Y53amp = 0.34594371914684025
    Y55amp = 0.46413220344085826
    Y55 = 0.46413220344085826*np.exp(-5*1j*phi)
    Y42amp = 0.33452327177864466
    Y44amp = 0.4425326924449826
    Y44 = 0.4425326924449826*np.exp(-4*1j*phi)
    Y31amp = 0.3231801841141506
    Y33amp = 0.4172238236327842
    Y33 = 0.4172238236327842*np.exp(-3*1j*phi)
    Y11amp = 0.3454941494713355
    Y11 = 0.3454941494713355*np.exp(-1j*phi)
    Y22amp = 0.3862742020231896
    Y22 = 0.3862742020231896*np.exp(-2*1j*phi)
    c9 = ((1 - delta)/2)**8 + ((-1)**9)*((1 + delta)/2)**8
    c8 = ((1 - delta)/2)**7 + ((-1)**8)*((1 + delta)/2)**7
    c7 = ((1 - delta)/2)**6 + ((-1)**7)*((1 + delta)/2)**6
    c6 = ((1 - delta)/2)**5 + ((-1)**6)*((1 + delta)/2)**5
    c5 = (((1 - delta)/2)**4 + ((-1)**5)*((1 + delta)/2)**4)*noneqcond + (-np.divide(1,2))*eqcond
    c4 = ((1 - delta)/2)**3 + ((-1)**4)*((1 + delta)/2)**3
    c3 = (((1 - delta)/2)**2 + ((-1)**3)*((1 + delta)/2)**2)*noneqcond + (-1)*eqcond
    c2 = ((1 - delta)/2)**1 + ((-1)**2)*((1 + delta)/2)**1
    N81 = (-16*1j*np.pi*(1j*1)**8/factorial2(2*8+1))*np.sqrt( ( (2*8+1)*(8+2)*(8**2 - 1**2) ) / ( (2*8-1)*(8+1)*(8)*(8-1) ) )
    N82 = (8*1j*np.pi*(1j*2)**8/factorial2(2*8 + 1))*np.sqrt( ( (8+1)*(8+2) ) / ( (8)*(8-1) ) )
    N83 = (-16*1j*np.pi*(1j*3)**8/factorial2(2*8+1))*np.sqrt( ( (2*8+1)*(8+2)*(8**2 - 3**2) ) / ( (2*8-1)*(8+1)*(8)*(8-1) ) )
    N84 = (8*1j*np.pi*(1j*4)**8/factorial2(2*8 + 1))*np.sqrt( ( (8+1)*(8+2) ) / ( (8)*(8-1) ) )
    N85 = (-16*1j*np.pi*(1j*5)**8/factorial2(2*8+1))*np.sqrt( ( (2*8+1)*(8+2)*(8**2 - 5**2) ) / ( (2*8-1)*(8+1)*(8)*(8-1) ) )
    N86 = (8*1j*np.pi*(1j*6)**8/factorial2(2*8 + 1))*np.sqrt( ( (8+1)*(8+2) ) / ( (8)*(8-1) ) )
    N87 = (-16*1j*np.pi*(1j*7)**8/factorial2(2*8+1))*np.sqrt( ( (2*8+1)*(8+2)*(8**2 - 7**2) ) / ( (2*8-1)*(8+1)*(8)*(8-1) ) )
    N88 = (8*1j*np.pi*(1j*8)**8/factorial2(2*8 + 1))*np.sqrt( ( (8+1)*(8+2) ) / ( (8)*(8-1) ) )
    N71 = (8*1j*np.pi*(1j*1)**7/factorial2(2*7 + 1))*np.sqrt( ( (7+1)*(7+2) ) / ( (7)*(7-1) ) )
    N72 = (-16*1j*np.pi*(1j*2)**7/factorial2(2*7+1))*np.sqrt( ( (2*7+1)*(7+2)*(7**2 - 2**2) ) / ( (2*7-1)*(7+1)*(7)*(7-1) ) )
    N73 = (8*1j*np.pi*(1j*3)**7/factorial2(2*7 + 1))*np.sqrt( ( (7+1)*(7+2) ) / ( (7)*(7-1) ) )
    N74 = (-16*1j*np.pi*(1j*4)**7/factorial2(2*7+1))*np.sqrt( ( (2*7+1)*(7+2)*(7**2 - 4**2) ) / ( (2*7-1)*(7+1)*(7)*(7-1) ) )
    N75 = (8*1j*np.pi*(1j*5)**7/factorial2(2*7 + 1))*np.sqrt( ( (7+1)*(7+2) ) / ( (7)*(7-1) ) )
    N76 = (-16*1j*np.pi*(1j*6)**7/factorial2(2*7+1))*np.sqrt( ( (2*7+1)*(7+2)*(7**2 - 6**2) ) / ( (2*7-1)*(7+1)*(7)*(7-1) ) )
    N77 = (8*1j*np.pi*(1j*7)**7/factorial2(2*7 + 1))*np.sqrt( ( (7+1)*(7+2) ) / ( (7)*(7-1) ) )
    N61 = (-16*1j*np.pi*(1j*1)**6/factorial2(2*6+1))*np.sqrt( ( (2*6+1)*(6+2)*(6**2 - 1**2) ) / ( (2*6-1)*(6+1)*(6)*(6-1) ) )
    N62 = (8*1j*np.pi*(1j*2)**6/factorial2(2*6 + 1))*np.sqrt( ( (6+1)*(6+2) ) / ( (6)*(6-1) ) )
    N63 = (-16*1j*np.pi*(1j*3)**6/factorial2(2*6+1))*np.sqrt( ( (2*6+1)*(6+2)*(6**2 - 3**2) ) / ( (2*6-1)*(6+1)*(6)*(6-1) ) )
    N64 = (8*1j*np.pi*(1j*4)**6/factorial2(2*6 + 1))*np.sqrt( ( (6+1)*(6+2) ) / ( (6)*(6-1) ) )
    N65 = (-16*1j*np.pi*(1j*5)**6/factorial2(2*6+1))*np.sqrt( ( (2*6+1)*(6+2)*(6**2 - 5**2) ) / ( (2*6-1)*(6+1)*(6)*(6-1) ) )
    N66 = (8*1j*np.pi*(1j*6)**6/factorial2(2*6 + 1))*np.sqrt( ( (6+1)*(6+2) ) / ( (6)*(6-1) ) )
    N51 = (8*1j*np.pi*(1j*1)**5/factorial2(2*5 + 1))*np.sqrt( ( (5+1)*(5+2) ) / ( (5)*(5-1) ) )
    N52 = (-16*1j*np.pi*(1j*2)**5/factorial2(2*5+1))*np.sqrt( ( (2*5+1)*(5+2)*(5**2 - 2**2) ) / ( (2*5-1)*(5+1)*(5)*(5-1) ) )
    N53 = (8*1j*np.pi*(1j*3)**5/factorial2(2*5 + 1))*np.sqrt( ( (5+1)*(5+2) ) / ( (5)*(5-1) ) )
    N54 = (-16*1j*np.pi*(1j*4)**5/factorial2(2*5+1))*np.sqrt( ( (2*5+1)*(5+2)*(5**2 - 4**2) ) / ( (2*5-1)*(5+1)*(5)*(5-1) ) )
    N55 = (8*1j*np.pi*(1j*5)**5/factorial2(2*5 + 1))*np.sqrt( ( (5+1)*(5+2) ) / ( (5)*(5-1) ) )
    N41 = (-16*1j*np.pi*(1j*1)**4/factorial2(2*4+1))*np.sqrt( ( (2*4+1)*(4+2)*(4**2 - 1**2) ) / ( (2*4-1)*(4+1)*(4)*(4-1) ) )
    N42 = (8*1j*np.pi*(1j*2)**4/factorial2(2*4 + 1))*np.sqrt( ( (4+1)*(4+2) ) / ( (4)*(4-1) ) )
    N43 = (-16*1j*np.pi*(1j*3)**4/factorial2(2*4+1))*np.sqrt( ( (2*4+1)*(4+2)*(4**2 - 3**2) ) / ( (2*4-1)*(4+1)*(4)*(4-1) ) )
    N44 = (8*1j*np.pi*(1j*4)**4/factorial2(2*4 + 1))*np.sqrt( ( (4+1)*(4+2) ) / ( (4)*(4-1) ) )
    N31 = (8*1j*np.pi*(1j*1)**3/factorial2(2*3 + 1))*np.sqrt( ( (3+1)*(3+2) ) / ( (3)*(3-1) ) )
    N32 = (-16*1j*np.pi*(1j*2)**3/factorial2(2*3+1))*np.sqrt( ( (2*3+1)*(3+2)*(3**2 - 2**2) ) / ( (2*3-1)*(3+1)*(3)*(3-1) ) )
    N33 = (8*1j*np.pi*(1j*3)**3/factorial2(2*3 + 1))*np.sqrt( ( (3+1)*(3+2) ) / ( (3)*(3-1) ) )
    N21 = (-16*1j*np.pi*(1j*1)**2/factorial2(2*2+1))*np.sqrt( ( (2*2+1)*(2+2)*(2**2 - 1**2) ) / ( (2*2-1)*(2+1)*(2)*(2-1) ) )
    N22 = (8*1j*np.pi*(1j*2)**2/factorial2(2*2 + 1))*np.sqrt( ( (2+1)*(2+2) ) / ( (2)*(2-1) ) )
    N81amp = np.absolute(N81)
    N82amp = np.absolute(N82)
    N83amp = np.absolute(N83)
    N84amp = np.absolute(N84)
    N85amp = np.absolute(N85)
    N86amp = np.absolute(N86)
    N87amp = np.absolute(N87)
    N88amp = np.absolute(N88)
    N71amp = np.absolute(N71)
    N72amp = np.absolute(N72)
    N73amp = np.absolute(N73)
    N74amp = np.absolute(N74)
    N75amp = np.absolute(N75)
    N76amp = np.absolute(N76)
    N77amp = np.absolute(N77)
    N61amp = np.absolute(N61)
    N62amp = np.absolute(N62)
    N63amp = np.absolute(N63)
    N64amp = np.absolute(N64)
    N65amp = np.absolute(N65)
    N66amp = np.absolute(N66)
    N51amp = np.absolute(N51)
    N52amp = np.absolute(N52)
    N53amp = np.absolute(N53)
    N54amp = np.absolute(N54)
    N55amp = np.absolute(N55)
    N41amp = np.absolute(N41)
    N42amp = np.absolute(N42)
    N43amp = np.absolute(N43)
    N44amp = np.absolute(N44)
    N31amp = np.absolute(N31)
    N32amp = np.absolute(N32)
    N33amp = np.absolute(N33)
    N21amp = np.absolute(N21)
    N22amp = np.absolute(N22)
    r0 = 2/np.sqrt(np.exp(1))
    b8 = -2*khat8
    b7 = -2*khat7
    b6 = -2*khat6
    b5 = -2*khat5
    b4 = -2*khat4
    b3 = -2*khat3
    b2 = -2*khat2
    b1 = -2*khat1
    T81prodfac = (1**2 + b1**2)*(2**2 + b1**2)*(3**2 + b1**2)*(4**2 + b1**2)*(5**2 + b1**2)*(6**2 + b1**2)*(7**2 + b1**2)*(8**2 + b1**2)
    T82prodfac = (1**2 + b2**2)*(2**2 + b2**2)*(3**2 + b2**2)*(4**2 + b2**2)*(5**2 + b2**2)*(6**2 + b2**2)*(7**2 + b2**2)*(8**2 + b2**2)
    T83prodfac = (1**2 + b3**2)*(2**2 + b3**2)*(3**2 + b3**2)*(4**2 + b3**2)*(5**2 + b3**2)*(6**2 + b3**2)*(7**2 + b3**2)*(8**2 + b3**2)
    T84prodfac = (1**2 + b4**2)*(2**2 + b4**2)*(3**2 + b4**2)*(4**2 + b4**2)*(5**2 + b4**2)*(6**2 + b4**2)*(7**2 + b4**2)*(8**2 + b4**2)
    T85prodfac = (1**2 + b5**2)*(2**2 + b5**2)*(3**2 + b5**2)*(4**2 + b5**2)*(5**2 + b5**2)*(6**2 + b5**2)*(7**2 + b5**2)*(8**2 + b5**2)
    T86prodfac = (1**2 + b6**2)*(2**2 + b6**2)*(3**2 + b6**2)*(4**2 + b6**2)*(5**2 + b6**2)*(6**2 + b6**2)*(7**2 + b6**2)*(8**2 + b6**2)
    T87prodfac = (1**2 + b7**2)*(2**2 + b7**2)*(3**2 + b7**2)*(4**2 + b7**2)*(5**2 + b7**2)*(6**2 + b7**2)*(7**2 + b7**2)*(8**2 + b7**2)
    T88prodfac = (1**2 + b8**2)*(2**2 + b8**2)*(3**2 + b8**2)*(4**2 + b8**2)*(5**2 + b8**2)*(6**2 + b8**2)*(7**2 + b8**2)*(8**2 + b8**2)
    T71prodfac = (1**2 + b1**2)*(2**2 + b1**2)*(3**2 + b1**2)*(4**2 + b1**2)*(5**2 + b1**2)*(6**2 + b1**2)*(7**2 + b1**2)
    T72prodfac = (1**2 + b2**2)*(2**2 + b2**2)*(3**2 + b2**2)*(4**2 + b2**2)*(5**2 + b2**2)*(6**2 + b2**2)*(7**2 + b2**2)
    T73prodfac = (1**2 + b3**2)*(2**2 + b3**2)*(3**2 + b3**2)*(4**2 + b3**2)*(5**2 + b3**2)*(6**2 + b3**2)*(7**2 + b3**2)
    T74prodfac = (1**2 + b4**2)*(2**2 + b4**2)*(3**2 + b4**2)*(4**2 + b4**2)*(5**2 + b4**2)*(6**2 + b4**2)*(7**2 + b4**2)
    T75prodfac = (1**2 + b5**2)*(2**2 + b5**2)*(3**2 + b5**2)*(4**2 + b5**2)*(5**2 + b5**2)*(6**2 + b5**2)*(7**2 + b5**2)
    T76prodfac = (1**2 + b6**2)*(2**2 + b6**2)*(3**2 + b6**2)*(4**2 + b6**2)*(5**2 + b6**2)*(6**2 + b6**2)*(7**2 + b6**2)
    T77prodfac = (1**2 + b7**2)*(2**2 + b7**2)*(3**2 + b7**2)*(4**2 + b7**2)*(5**2 + b7**2)*(6**2 + b7**2)*(7**2 + b7**2)
    T61prodfac = (1**2 + b1**2)*(2**2 + b1**2)*(3**2 + b1**2)*(4**2 + b1**2)*(5**2 + b1**2)*(6**2 + b1**2)
    T62prodfac = (1**2 + b2**2)*(2**2 + b2**2)*(3**2 + b2**2)*(4**2 + b2**2)*(5**2 + b2**2)*(6**2 + b2**2)
    T63prodfac = (1**2 + b3**2)*(2**2 + b3**2)*(3**2 + b3**2)*(4**2 + b3**2)*(5**2 + b3**2)*(6**2 + b3**2)
    T64prodfac = (1**2 + b4**2)*(2**2 + b4**2)*(3**2 + b4**2)*(4**2 + b4**2)*(5**2 + b4**2)*(6**2 + b4**2)
    T65prodfac = (1**2 + b5**2)*(2**2 + b5**2)*(3**2 + b5**2)*(4**2 + b5**2)*(5**2 + b5**2)*(6**2 + b5**2)
    T66prodfac = (1**2 + b6**2)*(2**2 + b6**2)*(3**2 + b6**2)*(4**2 + b6**2)*(5**2 + b6**2)*(6**2 + b6**2)
    T51prodfac = (1**2 + b1**2)*(2**2 + b1**2)*(3**2 + b1**2)*(4**2 + b1**2)*(5**2 + b1**2)
    T52prodfac = (1**2 + b2**2)*(2**2 + b2**2)*(3**2 + b2**2)*(4**2 + b2**2)*(5**2 + b2**2)
    T53prodfac = (1**2 + b3**2)*(2**2 + b3**2)*(3**2 + b3**2)*(4**2 + b3**2)*(5**2 + b3**2)
    T54prodfac = (1**2 + b4**2)*(2**2 + b4**2)*(3**2 + b4**2)*(4**2 + b4**2)*(5**2 + b4**2)
    T55prodfac = (1**2 + b5**2)*(2**2 + b5**2)*(3**2 + b5**2)*(4**2 + b5**2)*(5**2 + b5**2)
    T41prodfac = (1**2 + b1**2)*(2**2 + b1**2)*(3**2 + b1**2)*(4**2 + b1**2)
    T42prodfac = (1**2 + b2**2)*(2**2 + b2**2)*(3**2 + b2**2)*(4**2 + b2**2)
    T43prodfac = (1**2 + b3**2)*(2**2 + b3**2)*(3**2 + b3**2)*(4**2 + b3**2)
    T44prodfac = (1**2 + b4**2)*(2**2 + b4**2)*(3**2 + b4**2)*(4**2 + b4**2)
    T31prodfac = (1**2 + b1**2)*(2**2 + b1**2)*(3**2 + b1**2)
    T32prodfac = (1**2 + b2**2)*(2**2 + b2**2)*(3**2 + b2**2)
    T33prodfac = (1**2 + b3**2)*(2**2 + b3**2)*(3**2 + b3**2)
    T21prodfac = (1**2 + b1**2)*(2**2 + b1**2)
    T22prodfac = (1**2 + b2**2)*(2**2 + b2**2)
    T8prefac = 2*np.pi*b8/(np.exp(np.pi*b8) - np.exp(-np.pi*b8))
    T7prefac = 2*np.pi*b7/(np.exp(np.pi*b7) - np.exp(-np.pi*b7))
    T6prefac = 2*np.pi*b6/(np.exp(np.pi*b6) - np.exp(-np.pi*b6))
    T5prefac = 2*np.pi*b5/(np.exp(np.pi*b5) - np.exp(-np.pi*b5))
    T4prefac = 2*np.pi*b4/(np.exp(np.pi*b4) - np.exp(-np.pi*b4))
    T3prefac = 2*np.pi*b3/(np.exp(np.pi*b3) - np.exp(-np.pi*b3))
    T2prefac = 2*np.pi*b2/(np.exp(np.pi*b2) - np.exp(-np.pi*b2))
    T1prefac = 2*np.pi*b1/(np.exp(np.pi*b1) - np.exp(-np.pi*b1))
    gamma_amp_81 = np.sqrt(T1prefac*T81prodfac)
    gamma_amp_82 = np.sqrt(T2prefac*T82prodfac)
    gamma_amp_83 = np.sqrt(T3prefac*T83prodfac)
    gamma_amp_84 = np.sqrt(T4prefac*T84prodfac)
    gamma_amp_85 = np.sqrt(T5prefac*T85prodfac)
    gamma_amp_86 = np.sqrt(T6prefac*T86prodfac)
    gamma_amp_87 = np.sqrt(T7prefac*T87prodfac)
    gamma_amp_88 = np.sqrt(T8prefac*T88prodfac)
    gamma_amp_71 = np.sqrt(T1prefac*T71prodfac)
    gamma_amp_72 = np.sqrt(T2prefac*T72prodfac)
    gamma_amp_73 = np.sqrt(T3prefac*T73prodfac)
    gamma_amp_74 = np.sqrt(T4prefac*T74prodfac)
    gamma_amp_75 = np.sqrt(T5prefac*T75prodfac)
    gamma_amp_76 = np.sqrt(T6prefac*T76prodfac)
    gamma_amp_77 = np.sqrt(T7prefac*T77prodfac)
    gamma_amp_61 = np.sqrt(T1prefac*T61prodfac)
    gamma_amp_62 = np.sqrt(T2prefac*T62prodfac)
    gamma_amp_63 = np.sqrt(T3prefac*T63prodfac)
    gamma_amp_64 = np.sqrt(T4prefac*T64prodfac)
    gamma_amp_65 = np.sqrt(T5prefac*T65prodfac)
    gamma_amp_66 = np.sqrt(T6prefac*T66prodfac)
    gamma_amp_51 = np.sqrt(T1prefac*T51prodfac)
    gamma_amp_52 = np.sqrt(T2prefac*T52prodfac)
    gamma_amp_53 = np.sqrt(T3prefac*T53prodfac)
    gamma_amp_54 = np.sqrt(T4prefac*T54prodfac)
    gamma_amp_55 = np.sqrt(T5prefac*T55prodfac)
    gamma_amp_41 = np.sqrt(T1prefac*T41prodfac)
    gamma_amp_42 = np.sqrt(T2prefac*T42prodfac)
    gamma_amp_43 = np.sqrt(T3prefac*T43prodfac)
    gamma_amp_44 = np.sqrt(T4prefac*T44prodfac)
    gamma_amp_31 = np.sqrt(T1prefac*T31prodfac)
    gamma_amp_32 = np.sqrt(T2prefac*T32prodfac)
    gamma_amp_33 = np.sqrt(T3prefac*T33prodfac)
    gamma_amp_21 = np.sqrt(T1prefac*T21prodfac)
    gamma_amp_22 = np.sqrt(T2prefac*T22prodfac)
    T81amp = gamma_amp_81*np.exp(np.pi*khat1)/factorial(8)
    T82amp = gamma_amp_82*np.exp(np.pi*khat2)/factorial(8)
    T83amp = gamma_amp_83*np.exp(np.pi*khat3)/factorial(8)
    T84amp = gamma_amp_84*np.exp(np.pi*khat4)/factorial(8)
    T85amp = gamma_amp_85*np.exp(np.pi*khat5)/factorial(8)
    T86amp = gamma_amp_86*np.exp(np.pi*khat6)/factorial(8)
    T87amp = gamma_amp_87*np.exp(np.pi*khat7)/factorial(8)
    T88amp = gamma_amp_88*np.exp(np.pi*khat8)/factorial(8)
    T71amp = gamma_amp_71*np.exp(np.pi*khat1)/factorial(7)
    T72amp = gamma_amp_72*np.exp(np.pi*khat2)/factorial(7)
    T73amp = gamma_amp_73*np.exp(np.pi*khat3)/factorial(7)
    T74amp = gamma_amp_74*np.exp(np.pi*khat4)/factorial(7)
    T75amp = gamma_amp_75*np.exp(np.pi*khat5)/factorial(7)
    T76amp = gamma_amp_76*np.exp(np.pi*khat6)/factorial(7)
    T77amp = gamma_amp_77*np.exp(np.pi*khat7)/factorial(7)
    T61amp = gamma_amp_61*np.exp(np.pi*khat1)/factorial(6)
    T62amp = gamma_amp_62*np.exp(np.pi*khat2)/factorial(6)
    T63amp = gamma_amp_63*np.exp(np.pi*khat3)/factorial(6)
    T64amp = gamma_amp_64*np.exp(np.pi*khat4)/factorial(6)
    T65amp = gamma_amp_65*np.exp(np.pi*khat5)/factorial(6)
    T66amp = gamma_amp_66*np.exp(np.pi*khat6)/factorial(6)
    T51amp = gamma_amp_51*np.exp(np.pi*khat1)/factorial(5)
    T52amp = gamma_amp_52*np.exp(np.pi*khat2)/factorial(5)
    T53amp = gamma_amp_53*np.exp(np.pi*khat3)/factorial(5)
    T54amp = gamma_amp_54*np.exp(np.pi*khat4)/factorial(5)
    T55amp = gamma_amp_55*np.exp(np.pi*khat5)/factorial(5)
    T55 = gamma(5 + 1 - 2*1j*khat5)*np.exp(np.pi*khat5)*(np.exp(2*1j*khat5*np.log(2*5*Omega*r0)))/factorial(5)
    T41amp = gamma_amp_41*np.exp(np.pi*khat1)/factorial(4)
    T42amp = gamma_amp_42*np.exp(np.pi*khat2)/factorial(4)
    T43amp = gamma_amp_43*np.exp(np.pi*khat3)/factorial(4)
    T43 = gamma(4 + 1 - 2*1j*khat3)*np.exp(np.pi*khat3)*(np.exp(2*1j*khat3*np.log(2*3*Omega*r0)))/factorial(4)
    T44amp = gamma_amp_44*np.exp(np.pi*khat4)/factorial(4)
    T44 = gamma(4 + 1 - 2*1j*khat4)*np.exp(np.pi*khat4)*(np.exp(2*1j*khat4*np.log(2*4*Omega*r0)))/factorial(4)
    T31amp = gamma_amp_31*np.exp(np.pi*khat1)/factorial(3)
    T32amp = gamma_amp_32*np.exp(np.pi*khat2)/factorial(3)
    T32 = gamma(3 + 1 - 2*1j*khat2)*np.exp(np.pi*khat2)*(np.exp(2*1j*khat2*np.log(2*2*Omega*r0)))/factorial(3)
    T33amp = gamma_amp_33*np.exp(np.pi*khat3)/factorial(3)
    T33 = gamma(3 + 1 - 2*1j*khat3)*np.exp(np.pi*khat3)*(np.exp(2*1j*khat3*np.log(2*3*Omega*r0)))/factorial(3)
    T21amp = gamma_amp_21*np.exp(np.pi*khat1)/factorial(2)
    T21 = gamma(2 + 1 - 2*1j*khat1)*np.exp(np.pi*khat1)*(np.exp(2*1j*khat1*np.log(2*1*Omega*r0)))/factorial(2)
    T22amp = gamma_amp_22*np.exp(np.pi*khat2)/factorial(2)
    T22 = gamma(2 + 1 - 2*1j*khat2)*np.exp(np.pi*khat2)*(np.exp(2*1j*khat2*np.log(2*2*Omega*r0)))/factorial(2)
    source_odd = vomega*pphi
    source_even = Heff
    hN81amp = eta*N81amp*np.absolute(c9)*vphi**9*Y71amp
    hN82amp = eta*N82amp*np.absolute(c8)*vphi**8*Y82amp
    hN83amp = eta*N83amp*np.absolute(c9)*vphi**9*Y73amp
    hN84amp = eta*N84amp*np.absolute(c8)*vphi**8*Y84amp
    hN85amp = eta*N85amp*np.absolute(c9)*vphi**9*Y75amp
    hN86amp = eta*N86amp*np.absolute(c8)*vphi**8*Y86amp
    hN87amp = eta*N87amp*np.absolute(c9)*vphi**9*Y77amp
    hN88amp = eta*N88amp*np.absolute(c8)*vphi**8*Y88amp
    hN71amp = eta*N71amp*np.absolute(c7)*vphi**7*Y71amp
    hN72amp = eta*N72amp*np.absolute(c8)*vphi**8*Y62amp
    hN73amp = eta*N73amp*np.absolute(c7)*vphi**7*Y73amp
    hN74amp = eta*N74amp*np.absolute(c8)*vphi**8*Y64amp
    hN75amp = eta*N75amp*np.absolute(c7)*vphi**7*Y75amp
    hN76amp = eta*N76amp*np.absolute(c8)*vphi**8*Y66amp
    hN77amp = eta*N77amp*np.absolute(c7)*vphi**7*Y77amp
    hN61amp = eta*N61amp*np.absolute(c7)*vphi**7*Y51amp
    hN62amp = eta*N62amp*np.absolute(c6)*vphi**6*Y62amp
    hN63amp = eta*N63amp*np.absolute(c7)*vphi**7*Y53amp
    hN64amp = eta*N64amp*np.absolute(c6)*vphi**6*Y64amp
    hN65amp = eta*N65amp*np.absolute(c7)*vphi**7*Y55amp
    hN66amp = eta*N66amp*np.absolute(c6)*vphi**6*Y66amp
    hN51amp = eta*N51amp*np.absolute(c5)*vphi**5*Y51amp
    hN52amp = eta*N52amp*np.absolute(c6)*vphi**6*Y42amp
    hN53amp = eta*N53amp*np.absolute(c5)*vphi**5*Y53amp
    hN54amp = eta*N54amp*np.absolute(c6)*vphi**6*Y44amp
    hN55amp = eta*N55amp*np.absolute(c5)*vphi**5*Y55amp
    hN55 = eta*N55*c5*vphi**5*Y55
    hN41amp = eta*N41amp*np.absolute(c5)*vphi**5*Y31amp
    hN42amp = eta*N42amp*np.absolute(c4)*vphi**4*Y42amp
    hN43amp = eta*N43amp*np.absolute(c5)*vphi**5*Y33amp
    hN43 = eta*N43*c5*vphi**5*Y33
    hN44amp = eta*N44amp*np.absolute(c4)*vphi**4*Y44amp
    hN44 = eta*N44*c4*vphi**4*Y44
    hN31amp = eta*N31amp*np.absolute(c3)*vphi**3*Y31amp
    hN32amp = eta*N32amp*np.absolute(c4)*vphi**4*Y22amp
    hN32 = eta*N32*c4*vphi**4*Y22
    hN33amp = eta*N33amp*np.absolute(c3)*vphi**3*Y33amp
    hN33 = eta*N33*c3*vphi**3*Y33
    hN21amp = eta*N21amp*np.absolute(c3)*vphi**3*Y11amp
    hN21 = eta*N21*c3*vphi**3*Y11
    hN22amp = eta*N22amp*np.absolute(c2)*vphi**2*Y22amp
    hN22 = eta*N22*c2*vphi**2*Y22
    h81amp = hN81amp*source_odd*T81amp*f81amp
    h82amp = hN82amp*source_even*T82amp*f82amp
    h83amp = hN83amp*source_odd*T83amp*f83amp
    h84amp = hN84amp*source_even*T84amp*f84amp
    h85amp = hN85amp*source_odd*T85amp*f85amp
    h86amp = hN86amp*source_even*T86amp*f86amp
    h87amp = hN87amp*source_odd*T87amp*f87amp
    h88amp = hN88amp*source_even*T88amp*f88amp
    h71amp = hN71amp*source_even*T71amp*f71amp
    h72amp = hN72amp*source_odd*T72amp*f72amp
    h73amp = hN73amp*source_even*T73amp*f73amp
    h74amp = hN74amp*source_odd*T74amp*f74amp
    h75amp = hN75amp*source_even*T75amp*f75amp
    h76amp = hN76amp*source_odd*T76amp*f76amp
    h77amp = hN77amp*source_even*T77amp*f77amp
    h61amp = hN61amp*source_odd*T61amp*f61amp
    h62amp = hN62amp*source_even*T62amp*f62amp
    h63amp = hN63amp*source_odd*T63amp*f63amp
    h64amp = hN64amp*source_even*T64amp*f64amp
    h65amp = hN65amp*source_odd*T65amp*f65amp
    h66amp = hN66amp*source_even*T66amp*f66amp
    h51amp = hN51amp*source_even*T51amp*f51amp
    h52amp = hN52amp*source_odd*T52amp*f52amp
    h53amp = hN53amp*source_even*T53amp*f53amp
    h54amp = hN54amp*source_odd*T54amp*f54amp
    h55amp = hN55amp*source_even*T55amp*f55amp
    h55 = hN55*source_even*T55*f55*np.exp(1j*delta55)
    h41amp = hN41amp*source_odd*T41amp*f41amp
    h42amp = hN42amp*source_even*T42amp*f42amp
    h43amp = hN43amp*source_odd*T43amp*f43amp
    h43 = hN43*source_odd*T43*f43*np.exp(1j*delta43)
    h44amp = hN44amp*source_even*T44amp*f44amp
    h44 = hN44*source_even*T44*f44*np.exp(1j*delta44)
    h31amp = hN31amp*source_even*T31amp*f31amp
    h32amp = hN32amp*source_odd*T32amp*f32amp
    h32 = hN32*source_odd*T32*f32*np.exp(1j*delta32)
    h33amp = hN33amp*source_even*T33amp*f33amp
    h33 = hN33*source_even*T33*f33*np.exp(1j*delta33)
    h21amp = hN21amp*source_odd*T21amp*f21amp
    h21 = hN21*source_odd*T21*f21*np.exp(1j*delta21)
    h22amp = hN22amp*source_even*T22amp*f22amp
    h22 = hN22*source_even*T22*f22*np.exp(1j*delta22)
    flux=-(Omega**2/8/np.pi)*(2*2*h22amp**2+1*1*h21amp**2+3*3*h33amp**2+2*2*h32amp**2+1*1*h31amp**2+4*4*h44amp**2+3*3*h43amp**2+2*2*h42amp**2+1*1*h41amp**2+5*5*h55amp**2+4*4*h54amp**2+3*3*h53amp**2+2*2*h52amp**2+1*1*h51amp**2+6*6*h66amp**2+5*5*h65amp**2+4*4*h64amp**2+3*3*h63amp**2+2*2*h62amp**2+1*1*h61amp**2+7*7*h77amp**2+6*6*h76amp**2+5*5*h75amp**2+4*4*h74amp**2+3*3*h73amp**2+2*2*h72amp**2+1*1*h71amp**2+8*8*h88amp**2+7*7*h87amp**2+6*6*h86amp**2+5*5*h85amp**2+4*4*h84amp**2+3*3*h83amp**2+2*2*h82amp**2+1*1*h81amp**2)
    if not verbose:
        return flux
    else:
        return np.real(flux), np.array([[0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0],[0,f21amp,f22amp,0,0,0,0,0,0],[0,f31amp,f32amp,f33amp,0,0,0,0,0],[0,f41amp,f42amp,f43amp,f44amp,0,0,0,0],[0,f51amp,f52amp,f53amp,f54amp,f55amp,0,0,0],[0,f61amp,f62amp,f63amp,f64amp,f65amp,f66amp,0,0],[0,f71amp,f72amp,f73amp,f74amp,f75amp,f76amp,f77amp,0],[0,f81amp,f82amp,f83amp,f84amp,f85amp,f86amp,f87amp,f88amp]]),np.array([[0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0],[0,rhoNS21,rho22,0,0,0,0,0,0],[0,rhoNS31,rho32,rhoNS33,0,0,0,0,0],[0,rhoNS41,rho42,rhoNS43,rho44,0,0,0,0],[0,rhoNS51,rho52,rhoNS53,rho54,rhoNS55,0,0,0],[0,rhoNS61,rho62,rhoNS63,rho64,rhoNS65,rho66,0,0],[0,rhoNS71,rho72,rhoNS73,rho74,rhoNS75,rho76,rhoNS77,0],[0,rhoNS81,rho82,rhoNS83,rho84,rhoNS85,rho86,rhoNS87,rho88]]), np.array([[0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0],[0,fspin21,0,0,0,0,0,0,0],[0,fspin31,0,fspin33amp,0,0,0,0,0],[0,fspin41,0,fspin43,0,0,0,0,0],[0,0,0,0,0,fspin55,0,0,0],[0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0]]), np.array([[0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0],[0,T21amp,T22amp,0,0,0,0,0,0],[0,T31amp,T32amp,T33amp,0,0,0,0,0],[0,T41amp,T42amp,T43amp,T44amp,0,0,0,0],[0,T51amp,T52amp,T53amp,T54amp,T55amp,0,0,0],[0,T61amp,T62amp,T63amp,T64amp,T65amp,T66amp,0,0],[0,T71amp,T72amp,T73amp,T74amp,T75amp,T76amp,T77amp,0],[0,T81amp,T82amp,T83amp,T84amp,T85amp,T86amp,T87amp,T88amp]]), source_even, source_odd, delta, deltainv
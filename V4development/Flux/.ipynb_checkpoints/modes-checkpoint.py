import numpy as np
def compute_modes(m1=23., m2=10., a=-9.086823027135883e-03, chiA=-4.516044029838846e-02, chiS=1.506880542362110e-02):
    modes={}
    eta = (m1*m2)/(m1+m2)/(m1+m2)
    deltam = (m1 - m2)/(m1 + m2)
    deltam2=deltam*deltam
    a2 = a*a
    a3 = a2*a
    a4 = a3*a
    a5 = a4*a
    a6 = a5*a
    a7=a6*a
    eta2 = eta*eta
    eta3 = eta2*eta
    eta4=eta3*eta
    m1p3eta=-1+3*eta
    m1p3eta2 = m1p3eta*m1p3eta
    m1p3eta3 = m1p3eta2*m1p3eta
    
    modes['delta22vh3'] = 7./3.
    modes['delta22vh6'] = -4.*(deltam*chiA + chiS*(1. - 2.*eta))/3. + (428.*np.pi)/105.
    modes['delta22v8'] = 0.
    modes['delta22vh9'] = -2203./81. + 1712.*np.pi*np.pi/ 315.
    modes['delta22v5'] = -24.*eta
    modes['delta22v6'] = 0.
    modes['rho22v2'] = -43./42. + (55.*eta)/84.
    modes['rho22v3'] = (-2.*(chiS + chiA*deltam - chiS*eta))/3.
    modes['rho22v4'] = -20555./10584. + 0.5*(chiS + chiA*deltam)*(chiS + chiA*deltam) - (33025.*eta)/21168. + (19583.*eta2)/42336.
    modes['rho22v5'] = (-34./21. + 49.*eta/18. + 209.*eta2/126.)*chiS + (-34./21. - 19.*eta/42.)*deltam*chiA
    modes['rho22v6'] = 1556919113./122245200.+(89.*a2)/252.-(48993925.*eta)/9779616.-(6292061.*eta2)/3259872.+(10620745.*eta3)/39118464.+(41.*eta*np.pi*np.pi)/192.
    modes['rho22v6l'] = -428./105.
    modes['rho22v7'] = a3/3. + chiA*deltam*(18733./15876. + (50140.*eta)/3969. + (97865.*eta2)/63504.) + chiS*(18733./15876. + (74749.*eta)/5292. - (245717.*eta2)/63504. + (50803.*eta3)/63504.)
    modes['rho22v8'] = -387216563023./160190110080. + (18353.*a2)/21168. - a2*a2/8.
    modes['rho22v8l'] =  9202./2205.
    modes['rho22v10'] =  -16094530514677./533967033600.
    modes['rho22v10l'] = 439877./55566.
    
    modes['delta21vh3'] = 2./3.
    modes['delta21vh6'] = 107.*np.pi/105.
    modes['delta21ch7'] = 0.
    modes['delta21vh9'] = -272./81. + (214.*np.pi*np.pi)/315.
    modes['delta21v5'] = -493.*eta/42.
    modes['delta21v7'] = 0.0    
    modes['rho21v1'] = 0.    
    if(deltam2):
        modes['rho21v2'] =  -59./56. + (23.*eta)/84.
        modes['rho21v4'] =  -47009./56448. - (865.*a2)/1792. - (405.*a4)/2048. - (10993.*eta)/14112. + (617.*eta2)/4704.
        modes['rho21v5'] =  (-98635.*a)/75264. + (2031.*a3)/7168. - (1701.*a5)/8192.
        modes['rho21v6'] =  7613184941./2607897600. + (9032393.*a2)/1806336. + (3897.*a4)/16384. - (15309.*a6)/65536.
        modes['rho21v6l'] =  -107./105.
        modes['rho21v7'] =  (-3859374457.*a)/1159065600. - (55169.*a3)/16384. + (18603.*a5)/65536. - (72171.*a7)/262144.
        modes['rho21v7l'] =  107.*a/140.
        modes['rho21v8'] =  -1168617463883./911303737344.
        modes['rho21v8l'] =  6313./5880.
        modes['rho21v10'] =  -63735873771463./16569158860800.
        modes['rho21v10l'] = 5029963./5927040.
        modes['f21v1'] =  (-3.*(chiS + chiA/deltam))/2.
        modes['f21v3'] = (chiS*deltam*(427.+79.*eta)+chiA*(147.+280.*deltam2+1251.*eta))/84./deltam
    else:
        modes['rho21v2'] = 0.
        modes['rho21v3'] = 0.
        modes['rho21v4'] = 0.
        modes['rho21v5'] = 0.
        modes['rho21v6'] = 0.
        modes['rho21v6l'] = 0.
        modes['rho21v7'] = 0.
        modes['rho21v7l'] = 0.
        modes['rho21v8'] = 0.
        modes['rho21v8l'] = 0.
        modes['rho21v10'] = 0.
        modes['rho21v10l'] = 0.
        modes['f21v1'] = -3.*chiA/2.
        modes['f21v3'] = (chiS*deltam*(427. + 79.*eta) + chiA*(147. + 280.*deltam2 + 1251.*eta))/84.
     
    modes['delta33vh3'] = 13./10.
    modes['delta33vh6'] = (39.*np.pi)/7.
    modes['delta33vh9'] = -227827./3000. + (78.*np.pi*np.pi)/7.
    modes['delta33v5'] = -80897.*eta/2430.; 
    if(deltam2):
        modes['rho33v2'] = -7./6. + (2.*eta)/3.
        modes['rho33v3'] = 0.
        modes['rho33v4'] = -6719./3960. + a2/2. - (1861.*eta)/990. + (149.*eta2)/330.
        modes['rho33v5'] = (-4.*a)/3.
        modes['rho33v6'] = 3203101567./227026800. + (5.*a2)/36.
        modes['rho33v6l'] = -26./7.
        modes['rho33v7'] = (5297.*a)/2970. + a3/3.
        modes['rho33v8'] = -57566572157./8562153600.
        modes['rho33v8l'] = 13./3.
        modes['f33v3'] = (chiS*deltam*(-4.+5.*eta)+chiA*(-4.+19.*eta))/(2.*deltam)
    else:
        modes['rho33v2'] = 0.
        modes['rho33v3'] = 0.
        modes['rho33v4'] = 0.
        modes['rho33v5'] = 0.
        modes['rho33v6'] = 0.
        modes['rho33v6l'] = 0.
        modes['rho33v7'] = 0.
        modes['rho33v8'] = 0.
        modes['rho33v8l'] = 0.
        modes['f33v3'] = chiA*3./8.
    
    modes['delta32vh3'] = (10. + 33.*eta)/(-15.*m1p3eta)
    modes['delta32vh4'] = 0.
    modes['delta32vh6'] = (52.*np.pi)/21.
    modes['delta32vh9'] = -9112./405. + (208.*np.pi*np.pi)/63.
    modes['rho32v1'] = (4.*chiS*eta)/(-3.*m1p3eta)
    modes['rho32v2'] =  (328. - 1115.*eta + 320.*eta2)/(270.*m1p3eta)
    modes['rho32v3'] =  2./9.*a
    modes['rho32v4'] = a2/3.+(-1444528.+8050045.*eta-4725605.*eta2-20338960.*eta3+3085640.*eta4)/(1603800.*m1p3eta*m1p3eta)
    modes['rho32v5'] =  (-2788.*a)/1215.
    modes['rho32v6'] =  5849948554./940355325. + (488.*a2)/405.
    modes['rho32v6l'] =  -104./63.
    modes['rho32v8'] =  -10607269449358./3072140846775.
    modes['rho32v8l'] = 17056./8505.    
    
    if(deltam2):
        modes['delta31vh3'] = 13./30.
        modes['delta31vh6'] = (13.*np.pi)/21.
        modes['delta31vh7'] = 0.
        modes['delta31vh9'] = -227827./81000. + (26.*np.pi*np.pi)/63.
        modes['delta31v5'] = -17.*eta/10.
        modes['rho31v2'] =  -13./18. - (2.*eta)/9.
        modes['rho31v3'] =  0.0
        modes['rho31v4'] =  101./7128. - (5.*a2)/6. - (1685.*eta)/1782. - (829.*eta2)/1782.
        modes['rho31v5'] =  (4.*a)/9.
        modes['rho31v6'] =  11706720301./6129723600. - (49.*a2)/108.
        modes['rho31v6l'] =  -26./63.
        modes['rho31v7'] =  (-2579.*a)/5346. + a3/9.
        modes['rho31v8'] =  2606097992581./4854741091200.
        modes['rho31v8l'] = 169./567.
        modes['f31v3'] = (chiA*(-4.+11.*eta)+chiS*deltam*(-4.+13.*eta))/(2.*deltam)
    else:
        modes['delta31vh3'] = 0.
        modes['delta31vh6'] = 0.
        modes['delta31vh7'] = 0.
        modes['delta31vh9'] = 0.
        modes['delta31v5'] = 0.        
        modes['rho31v2'] =  0.
        modes['rho31v3'] =  0.
        modes['rho31v4'] =  0.
        modes['rho31v5'] =  0.
        modes['rho31v6'] =  0.
        modes['rho31v6l'] = 0.
        modes['rho31v7'] = 0.
        modes['rho31v8'] = 0.
        modes['rho31v8l'] = 0.
        modes['f31v3'] = -chiA*5./8.    
    
    modes['delta44vh3'] = (112. + 219.*eta)/(-120.*m1p3eta)
    modes['delta44vh6'] = (25136.*np.pi)/3465.
    modes['delta44vh9'] = 0.    
    modes['rho44v2'] =  (1614. - 5870.*eta + 2625.*eta2)/(1320.*m1p3eta)
    modes['rho44v3'] =  (chiA * (10. - 39.*eta)*deltam + chiS * (10. - 41.*eta + 42.*eta2)) / (15.*m1p3eta)
    modes['rho44v4'] = a2/2.+(-511573572.+2338945704.*eta-313857376.*eta2-6733146000.*eta3+1252563795.*eta4)/(317116800.*np.power(m1p3eta,2))
    modes['rho44v5'] =  (-69.*a)/55.
    modes['rho44v6'] =  16600939332793./1098809712000. + (217.*a2)/3960.
    modes['rho44v6l'] = -12568./3465.    
    
    if(deltam2):
        modes['delta43vh3'] = (486. + 4961.*eta)/(810.*(1. - 2.*eta))
        modes['delta43vh4'] = 0.
        modes['delta43vh6'] = 1571.*np.pi/385.        
        modes['rho43v1'] =  0.
        modes['rho43v2'] =  (222. - 547.*eta + 160.*eta2)/(176.*(-1. + 2.*eta))
        modes['rho43v4'] =  -6894273./7047040. + (3.*a2)/8.
        modes['rho43v5'] =  (-12113.*a)/6160.
        modes['rho43v6'] =  1664224207351./195343948800.
        modes['rho43v6l'] = -1571./770.    
        modes['f43v1'] = (5.*(chiA-chiS*deltam)*eta)/(2.*deltam*(-1.+2.*eta))
    else:
        modes['delta43vh3'] = 0.
        modes['delta43vh4'] = 0.
        modes['delta43vh6'] = 0.        
        modes['rho43v1'] =  0.
        modes['rho43v2'] =  0.
        modes['rho43v4'] =  0.
        modes['rho43v5'] =  0.
        modes['rho43v6'] =  0.
        modes['rho43v6l'] = 0.    
        modes['f43v1'] = 5.*chiA/4.    
    
    modes['delta42vh3'] = (7.*(1. + 6.*eta))/(-15.*m1p3eta)
    modes['delta42vh6'] = (6284.*np.pi)/3465.
    modes['rho42v2'] =  (1146. - 3530.*eta + 285.*eta2)/(1320.*m1p3eta)
    modes['rho42v3'] =  (chiA * (10. - 21.*eta)*deltam + chiS*(10. - 59.*eta + 78.*eta2))/(15.*m1p3eta)
    modes['rho42v4'] = a2/2.+(-114859044.+295834536.*eta+1204388696.*eta2-3047981160.*eta3-379526805.*eta4)/(317116800.*m1p3eta*m1p3eta)
    modes['rho42v5'] =  (-7.*a)/110.
    modes['rho42v6'] =  848238724511./219761942400. + (2323.*a2)/3960.
    modes['rho42v6l'] = -3142./3465.    
    
    if(deltam2):
        modes['delta41vh3'] = (2. + 507.*eta)/(10.*(1. - 2.*eta))
        modes['delta41vh4'] = 0.
        modes['delta41vh6'] = 1571.*np.pi/3465.        
        modes['rho41v1'] =  0.0
        modes['rho41v2'] =  (602. - 1385.*eta + 288.*eta2)/(528.*(-1. + 2.*eta))
        modes['rho41v4'] =  -7775491./21141120. + (3.* a2)/8.
        modes['rho41v5'] =  (-20033.*a)/55440. - (5*a3)/6.
        modes['rho41v6'] =  1227423222031./1758095539200.
        modes['rho41v6l'] = -1571./6930.    
        modes['f41v1'] = (5.*(chiA-chiS*deltam)*eta)/(2.*deltam*(-1.+2.*eta))
    else:
        modes['delta41vh3'] = 0.
        modes['delta41vh4'] = 0.
        modes['delta41vh6'] = 0.        
        modes['rho41v1'] =  0.
        modes['rho41v2'] =  0.
        modes['rho41v4'] =  0.
        modes['rho41v5'] =  0.
        modes['rho41v6'] =  0.
        modes['rho41v6l'] = 0.    
        modes['f41v1'] = -5.*chiA/4.        
    
    modes['delta55vh3'] = (96875. + 857528.*eta)/(131250.*(1. - 2.*eta))
    modes['delta55vh6'] = 0.
    modes['delta55vh9'] = 0.
    
    if(deltam2):
        modes['rho55v2'] =  (487. - 1298.*eta + 512.*eta2)/(390.*(-1. + 2.*eta))
        modes['rho55v3'] =  (-2.*a)/3.
        modes['rho55v4'] =  -3353747./2129400. + a2/2.
        modes['rho55v5'] = -241.*a/195.
        modes['rho55v6'] = 0.
    else:
        modes['rho55v2'] = 0.
        modes['rho55v3'] = 0.
        modes['rho55v4'] = 0.
        modes['rho55v5'] = 0.
        modes['rho55v6'] = 0.
    
    modes['delta54vh3'] = 8. / 15.
    modes['delta54vh4'] = 0.
    modes['rho54v2'] =  (-17448. + 96019.*eta - 127610.*eta2 + 33320.*eta3)/(13650.*(1. - 5.*eta + 5.*eta2))
    modes['rho54v3'] =  (-2.*a)/15.
    modes['rho54v4'] = -16213384./15526875.+(2.*a2)/5.
        
    if(deltam2):
        modes['delta53vh3'] = 31. / 70.
        modes['rho53v2'] =  (375. - 850.*eta + 176.*eta2)/(390.*(-1. + 2.*eta))
        modes['rho53v3'] =  (-2.*a)/3.
        modes['rho53v4'] =  -410833./709800. + a2/2.
        modes['rho53v5'] = -103.*a/325.
    else:
        modes['delta53vh3'] = 0.
        modes['rho53v2'] = 0.
        modes['rho53v3'] = 0.
        modes['rho53v4'] = 0.
        modes['rho53v5'] = 0.
        
    modes['delta52vh3'] = 4. / 15.
    modes['delta52vh4'] = 0.
    modes['rho52v2'] =  (-15828. + 84679.*eta - 104930.*eta2 + 21980.*eta3)/(13650.*(1. - 5.*eta + 5.*eta2))
    modes['rho52v3'] =  (-2.*a)/15.
    modes['rho52v4'] = -7187914./15526875.+(2.*a2)/5.
        
    if(deltam2):
        modes['delta51vh3'] = 31. / 210.
        modes['rho51v2'] =  (319. - 626.*eta + 8.*eta2)/(390.*(-1. + 2.*eta))
        modes['rho51v3'] =  (-2.*a)/3.
        modes['rho51v4'] =  -31877./304200. + a2/2.
        modes['rho51v5'] = 139.*a/975.
    else:
        modes['delta51vh3'] = 0.
        modes['rho51v2'] = 0.
        modes['rho51v3'] = 0.
        modes['rho51v4'] = 0.
        modes['rho51v5'] = 0.
        
    modes['delta66vh3'] = 43. / 70.
    modes['rho66v2'] =  (-106. + 602.*eta - 861.*eta2 + 273.*eta3)/(84.*(1. - 5.*eta + 5.*eta2))
    modes['rho66v3'] =  (-2.*a)/3.
    modes['rho66v4'] = -1025435./659736.+a2/2.
        
    if(deltam2):
        modes['delta65vh3'] = 10. / 21.
        modes['rho65v2'] =  (-185. + 838.*eta - 910.*eta2 + 220.*eta3)/(144.*(deltam2 + 3.*eta2))
        modes['rho65v3'] = -2.*a/9.
    else:
        modes['delta65vh3'] = 0.
        modes['rho65v2'] = 0.
        modes['rho65v3'] = 0.
        
    modes['delta64vh3'] = 43. / 105.
    modes['rho64v2'] =  (-86. + 462.*eta - 581.*eta2 + 133.*eta3)/(84.*(1. - 5.*eta + 5.*eta2))
    modes['rho64v3'] =  (-2.*a)/3.
    modes['rho64v4'] = -476887./659736.+a2/2.
        
    if(deltam2):
        modes['delta63vh3'] = 2. / 7.
        modes['rho63v2'] =  (-169. + 742.*eta - 750.*eta2 + 156.*eta3)/(144.*(deltam2 + 3.*eta2))
        modes['rho63v3'] = -2.*a/9.
    else:
        modes['delta63vh3'] = 0.
        modes['rho63v2'] = 0.
        modes['rho63v3'] = 0.
        
    modes['delta62vh3'] = 43. / 210.
    modes['rho62v2'] =  (-74. + 378.*eta - 413.*eta2 + 49.*eta3)/(84.*(1. - 5.*eta + 5.*eta2))
    modes['rho62v3'] =  (-2.*a)/3.
    modes['rho62v4'] = -817991./3298680.+a2/2.
        
    if(deltam2):
        modes['delta61vh3'] = 2. / 21.
        modes['rho61v2'] =  (-161. + 694.*eta - 670.*eta2 + 124.*eta3)/(144.*(deltam2 + 3.*eta2))
        modes['rho61v3'] = -2.*a/9.
    else:
        modes['delta61vh3'] = 0.
        modes['rho61v2'] = 0.
        modes['rho61v3'] = 0.
        
    if(deltam2):
        modes['delta77vh3'] = 19. / 36.
        modes['rho77v2'] =  (-906. + 4246.*eta - 4963.*eta2 + 1380.*eta3)/(714.*(deltam2 + 3.*eta2))
        modes['rho77v3'] = -2.*a/3.
    else:
        modes['delta77vh3'] = 0.
        modes['rho77v2'] =  0.
        modes['rho77v3'] = 0.
        
    modes['rho76v2'] = (2144.-16185.*eta+37828.*eta2-29351.*eta3+6104.*eta4)/(1666.*(-1+7*eta-14*eta2+7*eta3))
        
    if(deltam2):
        modes['delta75vh3'] = 95. / 252.
        modes['rho75v2'] =  (-762. + 3382.*eta - 3523.*eta2 + 804.*eta3)/(714.*(deltam2 + 3.*eta2))
        modes['rho75v3'] = -2.*a/3.
    else:
        modes['delta75vh3'] = 0.;
        modes['rho75v2'] = 0.
        modes['rho75v3'] = 0.
            
    modes['rho74v2'] = (17756.-131805.*eta+298872.*eta2-217959.*eta3+41076.*eta4)/(14994.*(-1.+7.*eta-14.*eta2+7.*eta3))
        
    if(deltam2):
        modes['delta73vh3'] = 19. / 84.
        modes['rho73v2'] =  (-666. + 2806.*eta - 2563.*eta2 + 420.*eta3)/(714.*(deltam2 + 3.*eta2))
        modes['rho73v3'] = -2.*a/3.
    else:
        modes['delta73vh3'] = 0.
        modes['rho73v2'] = 0.
        modes['rho73v3'] = 0.
        
    modes['rho72v2'] = (16832.-123489.*eta+273924.*eta2-190239.*eta3+32760.*eta4)/(14994.*(-1.+7.*eta-14.*eta2+7.*eta3))
        
    if(deltam2):
        modes['delta71vh3'] = 19. / 252.
        modes['rho71v2'] =  (-618. + 2518.*eta - 2083.*eta2 + 228.*eta3)/(714.*(deltam*deltam + 3.*eta2))
        modes['rho71v3'] = -2.*a/3.
    else:
        modes['delta71vh3'] = 0.
        modes['rho71v2'] = 0.
        modes['rho71v3'] = 0.
            
    modes['rho88v2'] = (3482.-26778.*eta+64659.*eta2-53445.*eta3+12243.*eta4)/(2736.*(-1.+7.*eta-14.*eta2+7.*eta3))
    modes['rho86v2'] = (1002.-7498.*eta+17269.*eta2-13055.*eta3+2653.*eta4)/(912.*(-1.+7.*eta-14.*eta2+7.*eta3))
    modes['rho84v2'] = (2666.-19434.*eta+42627.*eta2-28965.*eta3+4899.*eta4)/(2736.*(-1.+7.*eta-14.*eta2+7.*eta3))
    modes['rho82v2'] = (2462.-17598.*eta+37119.*eta2-22845.*eta3+3063.*eta4)/(2736.*(-1.+7.*eta-14.*eta2+7.*eta3))
    
    if(deltam2):        
        modes['rho81v2'] =  (20022.-126451.*eta+236922.*eta2 - 138430.*eta3 + 21640.*eta4)/(18240.*(-1. + 6.*eta - 10.*eta2 + 4.*eta3))
        modes['rho83v2'] = (20598.-131059.*eta+249018.*eta2-149950.*eta3+24520.*eta4)/(18240.*(-1.+6.*eta-10.*eta2+4.*eta3))        
        modes['rho85v2'] = (4350.-28055.*eta+54642.*eta2-34598.*eta3+6056.*eta4)/(3648.*(-1.+6.*eta-10.*eta2+4.*eta3))        
        modes['rho87v2'] = (23478.-154099.*eta+309498.*eta2-207550.*eta3+38920*eta4)/(18240.*(-1+6*eta-10*eta2+4*eta3))
    else:
        modes['rho81v2'] = 0.
        modes['rho83v2'] = 0.
        modes['rho85v2'] = 0.
        modes['rho87v2'] = 0.
                    
    return modes

## Un-Comment for Validation
#test=compute_modes(m1=6.969697e-01,m2=3.0303030303030304e-01,a=-8.9973238670953502e-03,chiA = -4.5068658220324853e-02, chiS = 1.5161158482118485e-02)
#for key, value in test.items():
#    print("%s = %.15e" % (key,value))
#    if key == 'rho22v10l':
#        break

#print("lalvalues")
#lalvalues = {'v2' : -8.8552188552188538e-01, 'v3' : 3.8634950966511088e-03, 'v4' : -2.2509514149671328e+00, 'v5' : 1.5733331438572230e-02, 'v6' : 1.2039564081033429e+01, 'v6l' : -4.0761904761904759e+00, 'v7' : -8.9243738723327609e-03, 'v8' : -2.4171612078198277e+00, 'v8l' : 4.1732426303854879e+00, 'v10' : -3.0141431028368640e+01, 'v10l' : 7.9162977360256273e+00}
#for key, value in lalvalues.items():
#    print("%s = %.15e"% (key,value))
    
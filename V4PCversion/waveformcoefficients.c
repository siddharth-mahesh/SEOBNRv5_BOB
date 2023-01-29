#include "mainv4pheader.h"
#include <string.h>

int compute_modes(WaveformCoefficients * const coeffs, double m1,double  m2,double eta, double  a,double  chiA,double chiS){
    
    double deltam = (m1 - m2)/(m1 + m2);
    double deltam2=deltam*deltam;
    double a2 = a*a;
    double a3 = a2*a;
    double a4 = a3*a;
    double a5 = a4*a;
    double a6 = a5*a;
    double a7=a6*a;
    double eta2 = eta*eta;
    double eta3 = eta2*eta;
    double eta4=eta3*eta;
    double m1p3eta=-1+3*eta;
    double m1p3eta2 = m1p3eta*m1p3eta;
    double m1p3eta3 = m1p3eta2*m1p3eta;
    memset(coeffs,0, sizeof(WaveformCoefficients));
    
    coeffs->delta22vh3 = 7./3.;
    coeffs->delta22vh6 = -4.*(deltam*chiA + chiS*(1. - 2.*eta))/3. + (428.*M_PI)/105.;
    coeffs->delta22v8 = 0.;
    coeffs->delta22vh9 = -2203./81. + 1712.*M_PI*M_PI/ 315.;
    coeffs->delta22v5 = -24.*eta;
    coeffs->delta22v6 = 0.;
    coeffs->rho22v2 = -43./42. + (55.*eta)/84.;
    coeffs->rho22v3 = (-2.*(chiS + chiA*deltam - chiS*eta))/3.;
    coeffs->rho22v4 = -20555./10584. + 0.5*(chiS + chiA*deltam)*(chiS + chiA*deltam) - (33025.*eta)/21168. + (19583.*eta2)/42336.;
    coeffs->rho22v5 = (-34./21. + 49.*eta/18. + 209.*eta2/126.)*chiS + (-34./21. - 19.*eta/42.)*deltam*chiA;
    coeffs->rho22v6 = 1556919113./122245200.+(89.*a2)/252.-(48993925.*eta)/9779616.-(6292061.*eta2)/3259872.+(10620745.*eta3)/39118464.+(41.*eta*M_PI*M_PI)/192.;
    coeffs->rho22v6l = -428./105.;
    coeffs->rho22v7 = a3/3. + chiA*deltam*(18733./15876. + (50140.*eta)/3969. + (97865.*eta2)/63504.) + chiS*(18733./15876. + (74749.*eta)/5292. - (245717.*eta2)/63504. + (50803.*eta3)/63504.);
    coeffs->rho22v8 = -387216563023./160190110080. + (18353.*a2)/21168. - a2*a2/8.;
    coeffs->rho22v8l =  9202./2205.;
    coeffs->rho22v10 =  -16094530514677./533967033600.;
    coeffs->rho22v10l = 439877./55566.;
    
    coeffs->delta21vh3 = 2./3.;
    coeffs->delta21vh6 = 107.*M_PI/105.;
    coeffs->delta21ch7 = 0.;
    coeffs->delta21vh9 = -272./81. + (214.*M_PI*M_PI)/315.;
    coeffs->delta21v5 = -493.*eta/42.;
    coeffs->delta21v7 = 0.;
    coeffs->rho21v1 = 0.;    
    if(deltam2){
        coeffs->rho21v2 =  -59./56. + (23.*eta)/84.;
        coeffs->rho21v4 =  -47009./56448. - (865.*a2)/1792. - (405.*a4)/2048. - (10993.*eta)/14112. + (617.*eta2)/4704.;
        coeffs->rho21v5 =  (-98635.*a)/75264. + (2031.*a3)/7168. - (1701.*a5)/8192.;
        coeffs->rho21v6 =  7613184941./2607897600. + (9032393.*a2)/1806336. + (3897.*a4)/16384. - (15309.*a6)/65536.;
        coeffs->rho21v6l =  -107./105.;
        coeffs->rho21v7 =  (-3859374457.*a)/1159065600. - (55169.*a3)/16384. + (18603.*a5)/65536. - (72171.*a7)/262144.;
        coeffs->rho21v7l =  107.*a/140.;
        coeffs->rho21v8 =  -1168617463883./911303737344.;
        coeffs->rho21v8l =  6313./5880.;
        coeffs->rho21v10 =  -63735873771463./16569158860800.;
        coeffs->rho21v10l = 5029963./5927040.;
        coeffs->f21v1 =  (-3.*(chiS + chiA/deltam))/2.;
        coeffs->f21v3 = (chiS*deltam*(427.+79.*eta)+chiA*(147.+280.*deltam2+1251.*eta))/84./deltam;
    }else{
        coeffs->f21v1 = -3.*chiA/2.;
        coeffs->f21v3 = (chiS*deltam*(427. + 79.*eta) + chiA*(147. + 280.*deltam2 + 1251.*eta))/84.;
    }
    
    coeffs->delta33vh3 = 13./10.;
    coeffs->delta33vh6 = (39.*M_PI)/7.;
    coeffs->delta33vh9 = -227827./3000. + (78.*M_PI*M_PI)/7.;
    coeffs->delta33v5 = -80897.*eta/2430.; 
    if(deltam2){
        coeffs->rho33v2 = -7./6. + (2.*eta)/3.;
        coeffs->rho33v3 = 0.;
        coeffs->rho33v4 = -6719./3960. + a2/2. - (1861.*eta)/990. + (149.*eta2)/330.;
        coeffs->rho33v5 = (-4.*a)/3.;
        coeffs->rho33v6 = 3203101567./227026800. + (5.*a2)/36.;
        coeffs->rho33v6l = -26./7.;
        coeffs->rho33v7 = (5297.*a)/2970. + a3/3.;
        coeffs->rho33v8 = -57566572157./8562153600.;
        coeffs->rho33v8l = 13./3.;
        coeffs->f33v3 = (chiS*deltam*(-4.+5.*eta)+chiA*(-4.+19.*eta))/(2.*deltam);
    }else{
        coeffs->f33v3 = chiA*3./8.;
    }
    
    coeffs->delta32vh3 = (10. + 33.*eta)/(-15.*m1p3eta);
    coeffs->delta32vh4 = 0.;
    coeffs->delta32vh6 = (52.*M_PI)/21.;
    coeffs->delta32vh9 = -9112./405. + (208.*M_PI*M_PI)/63.;
    coeffs->rho32v1 = (4.*chiS*eta)/(-3.*m1p3eta);
    coeffs->rho32v2 =  (328. - 1115.*eta + 320.*eta2)/(270.*m1p3eta);
    coeffs->rho32v3 =  2./9.*a;
    coeffs->rho32v4 = a2/3.+(-1444528.+8050045.*eta-4725605.*eta2-20338960.*eta3+3085640.*eta4)/(1603800.*m1p3eta*m1p3eta);
    coeffs->rho32v5 =  (-2788.*a)/1215.;
    coeffs->rho32v6 =  5849948554./940355325. + (488.*a2)/405.;
    coeffs->rho32v6l =  -104./63.;
    coeffs->rho32v8 =  -10607269449358./3072140846775.;
    coeffs->rho32v8l = 17056./8505.;    
    
    if(deltam2){
        coeffs->delta31vh3 = 13./30.;
        coeffs->delta31vh6 = (13.*M_PI)/21.;
        coeffs->delta31vh7 = 0.;
        coeffs->delta31vh9 = -227827./81000. + (26.*M_PI*M_PI)/63.;
        coeffs->delta31v5 = -17.*eta/10.;
        coeffs->rho31v2 =  -13./18. - (2.*eta)/9.;
        coeffs->rho31v3 =  0.;
        coeffs->rho31v4 =  101./7128. - (5.*a2)/6. - (1685.*eta)/1782. - (829.*eta2)/1782.;
        coeffs->rho31v5 =  (4.*a)/9.;
        coeffs->rho31v6 =  11706720301./6129723600. - (49.*a2)/108.;
        coeffs->rho31v6l =  -26./63.;
        coeffs->rho31v7 =  (-2579.*a)/5346. + a3/9.;
        coeffs->rho31v8 =  2606097992581./4854741091200.;
        coeffs->rho31v8l = 169./567.;
        coeffs->f31v3 = (chiA*(-4.+11.*eta)+chiS*deltam*(-4.+13.*eta))/(2.*deltam);
    }else{
        coeffs->f31v3 = -chiA*5./8.;    
    }
    
    coeffs->delta44vh3 = (112. + 219.*eta)/(-120.*m1p3eta);
    coeffs->delta44vh6 = (25136.*M_PI)/3465.;
    coeffs->delta44vh9 = 0.;    
    coeffs->rho44v2 =  (1614. - 5870.*eta + 2625.*eta2)/(1320.*m1p3eta);
    coeffs->rho44v3 =  (chiA * (10. - 39.*eta)*deltam + chiS * (10. - 41.*eta + 42.*eta2)) / (15.*m1p3eta);
    coeffs->rho44v4 = a2/2.+(-511573572.+2338945704.*eta-313857376.*eta2-6733146000.*eta3+1252563795.*eta4)/(317116800.*np.power(m1p3eta,2));
    coeffs->rho44v5 =  (-69.*a)/55.;
    coeffs->rho44v6 =  16600939332793./1098809712000. + (217.*a2)/3960.;
    coeffs->rho44v6l = -12568./3465.;    
    
    if(deltam2){
        coeffs->delta43vh3 = (486. + 4961.*eta)/(810.*(1. - 2.*eta));
        coeffs->delta43vh4 = 0.;
        coeffs->delta43vh6 = 1571.*M_PI/385.;
        coeffs->rho43v1 =  0.;
        coeffs->rho43v2 =  (222. - 547.*eta + 160.*eta2)/(176.*(-1. + 2.*eta));
        coeffs->rho43v4 =  -6894273./7047040. + (3.*a2)/8.;
        coeffs->rho43v5 =  (-12113.*a)/6160.;
        coeffs->rho43v6 =  1664224207351./195343948800.;
        coeffs->rho43v6l = -1571./770.;
        coeffs->f43v1 = (5.*(chiA-chiS*deltam)*eta)/(2.*deltam*(-1.+2.*eta));
    }else{
        coeffs->f43v1 = 5.*chiA/4.;
    }
    
    coeffs->delta42vh3 = (7.*(1. + 6.*eta))/(-15.*m1p3eta);
    coeffs->delta42vh6 = (6284.*M_PI)/3465.;
    coeffs->rho42v2 =  (1146. - 3530.*eta + 285.*eta2)/(1320.*m1p3eta);
    coeffs->rho42v3 =  (chiA * (10. - 21.*eta)*deltam + chiS*(10. - 59.*eta + 78.*eta2))/(15.*m1p3eta);
    coeffs->rho42v4 = a2/2.+(-114859044.+295834536.*eta+1204388696.*eta2-3047981160.*eta3-379526805.*eta4)/(317116800.*m1p3eta*m1p3eta);
    coeffs->rho42v5 =  (-7.*a)/110.;
    coeffs->rho42v6 =  848238724511./219761942400. + (2323.*a2)/3960.;
    coeffs->rho42v6l = -3142./3465.;
    
    if(deltam2){
        coeffs->delta41vh3 = (2. + 507.*eta)/(10.*(1. - 2.*eta));
        coeffs->delta41vh4 = 0.;
        coeffs->delta41vh6 = 1571.*M_PI/3465.;        
        coeffs->rho41v1 =  0.0;
        coeffs->rho41v2 =  (602. - 1385.*eta + 288.*eta2)/(528.*(-1. + 2.*eta));
        coeffs->rho41v4 =  -7775491./21141120. + (3.* a2)/8.;
        coeffs->rho41v5 =  (-20033.*a)/55440. - (5*a3)/6.;
        coeffs->rho41v6 =  1227423222031./1758095539200.;
        coeffs->rho41v6l = -1571./6930.;
        coeffs->f41v1 = (5.*(chiA-chiS*deltam)*eta)/(2.*deltam*(-1.+2.*eta));
    }else{
        coeffs->f41v1 = -5.*chiA/4.;
    }
    
    coeffs->delta55vh3 = (96875. + 857528.*eta)/(131250.*(1. - 2.*eta));
    coeffs->delta55vh6 = 0.;
    coeffs->delta55vh9 = 0.;
    
    if(deltam2){
        coeffs->rho55v2 =  (487. - 1298.*eta + 512.*eta2)/(390.*(-1. + 2.*eta));
        coeffs->rho55v3 =  (-2.*a)/3.;
        coeffs->rho55v4 =  -3353747./2129400. + a2/2.;
        coeffs->rho55v5 = -241.*a/195.;
        coeffs->rho55v6 = 0.;
    }
    
    coeffs->delta54vh3 = 8. / 15.;
    coeffs->delta54vh4 = 0.;
    coeffs->rho54v2 =  (-17448. + 96019.*eta - 127610.*eta2 + 33320.*eta3)/(13650.*(1. - 5.*eta + 5.*eta2));
    coeffs->rho54v3 =  (-2.*a)/15.;
    coeffs->rho54v4 = -16213384./15526875.+(2.*a2)/5.;
        
    if(deltam2){
        coeffs->delta53vh3 = 31. / 70.;
        coeffs->rho53v2 =  (375. - 850.*eta + 176.*eta2)/(390.*(-1. + 2.*eta));
        coeffs->rho53v3 =  (-2.*a)/3.;
        coeffs->rho53v4 =  -410833./709800. + a2/2.;
        coeffs->rho53v5 = -103.*a/325.;
    }
    
    coeffs->delta52vh3 = 4. / 15.;
    coeffs->delta52vh4 = 0.;
    coeffs->rho52v2 =  (-15828. + 84679.*eta - 104930.*eta2 + 21980.*eta3)/(13650.*(1. - 5.*eta + 5.*eta2));
    coeffs->rho52v3 =  (-2.*a)/15.;
    coeffs->rho52v4 = -7187914./15526875.+(2.*a2)/5.;
        
    if(deltam2){
        coeffs->delta51vh3 = 31. / 210.;
        coeffs->rho51v2 =  (319. - 626.*eta + 8.*eta2)/(390.*(-1. + 2.*eta));
        coeffs->rho51v3 =  (-2.*a)/3.;
        coeffs->rho51v4 =  -31877./304200. + a2/2.;
        coeffs->rho51v5 = 139.*a/975.;
    }
    
    coeffs->delta66vh3 = 43. / 70.;
    coeffs->rho66v2 =  (-106. + 602.*eta - 861.*eta2 + 273.*eta3)/(84.*(1. - 5.*eta + 5.*eta2));
    coeffs->rho66v3 =  (-2.*a)/3.;
    coeffs->rho66v4 = -1025435./659736.+a2/2.;
        
    if(deltam2){
        coeffs->delta65vh3 = 10. / 21.;
        coeffs->rho65v2 =  (-185. + 838.*eta - 910.*eta2 + 220.*eta3)/(144.*(deltam2 + 3.*eta2));
        coeffs->rho65v3 = -2.*a/9.;
        
    coeffs->delta64vh3 = 43. / 105.;
    coeffs->rho64v2 =  (-86. + 462.*eta - 581.*eta2 + 133.*eta3)/(84.*(1. - 5.*eta + 5.*eta2));
    coeffs->rho64v3 =  (-2.*a)/3.;
    coeffs->rho64v4 = -476887./659736.+a2/2.;
        
    if(deltam2){
        coeffs->delta63vh3 = 2. / 7.;
        coeffs->rho63v2 =  (-169. + 742.*eta - 750.*eta2 + 156.*eta3)/(144.*(deltam2 + 3.*eta2));
        coeffs->rho63v3 = -2.*a/9.;
    }
        
    coeffs->delta62vh3 = 43. / 210.;
    coeffs->rho62v2 =  (-74. + 378.*eta - 413.*eta2 + 49.*eta3)/(84.*(1. - 5.*eta + 5.*eta2));
    coeffs->rho62v3 =  (-2.*a)/3.;
    coeffs->rho62v4 = -817991./3298680.+a2/2.;
        
    if(deltam2){
        coeffs->delta61vh3 = 2. / 21.;
        coeffs->rho61v2 =  (-161. + 694.*eta - 670.*eta2 + 124.*eta3)/(144.*(deltam2 + 3.*eta2));
        coeffs->rho61v3 = -2.*a/9.;
    }
        
    if(deltam2){
        coeffs->delta77vh3 = 19. / 36.;
        coeffs->rho77v2 =  (-906. + 4246.*eta - 4963.*eta2 + 1380.*eta3)/(714.*(deltam2 + 3.*eta2));
        coeffs->rho77v3 = -2.*a/3.;
    }
    
    coeffs->rho76v2 = (2144.-16185.*eta+37828.*eta2-29351.*eta3+6104.*eta4)/(1666.*(-1+7*eta-14*eta2+7*eta3));
        
    if(deltam2){
        coeffs->delta75vh3 = 95. / 252.;
        coeffs->rho75v2 =  (-762. + 3382.*eta - 3523.*eta2 + 804.*eta3)/(714.*(deltam2 + 3.*eta2));
        coeffs->rho75v3 = -2.*a/3.;
    }
            
    coeffs->rho74v2 = (17756.-131805.*eta+298872.*eta2-217959.*eta3+41076.*eta4)/(14994.*(-1.+7.*eta-14.*eta2+7.*eta3));
        
    if(deltam2){
        coeffs->delta73vh3 = 19. / 84.;
        coeffs->rho73v2 =  (-666. + 2806.*eta - 2563.*eta2 + 420.*eta3)/(714.*(deltam2 + 3.*eta2));
        coeffs->rho73v3 = -2.*a/3.;
        
    coeffs->rho72v2 = (16832.-123489.*eta+273924.*eta2-190239.*eta3+32760.*eta4)/(14994.*(-1.+7.*eta-14.*eta2+7.*eta3));
        
    if(deltam2){
        coeffs->delta71vh3 = 19. / 252.;
        coeffs->rho71v2 =  (-618. + 2518.*eta - 2083.*eta2 + 228.*eta3)/(714.*(deltam*deltam + 3.*eta2));
        coeffs->rho71v3 = -2.*a/3.;
    }
            
    coeffs->rho88v2 = (3482.-26778.*eta+64659.*eta2-53445.*eta3+12243.*eta4)/(2736.*(-1.+7.*eta-14.*eta2+7.*eta3));
    coeffs->rho86v2 = (1002.-7498.*eta+17269.*eta2-13055.*eta3+2653.*eta4)/(912.*(-1.+7.*eta-14.*eta2+7.*eta3));
    coeffs->rho84v2 = (2666.-19434.*eta+42627.*eta2-28965.*eta3+4899.*eta4)/(2736.*(-1.+7.*eta-14.*eta2+7.*eta3));
    coeffs->rho82v2 = (2462.-17598.*eta+37119.*eta2-22845.*eta3+3063.*eta4)/(2736.*(-1.+7.*eta-14.*eta2+7.*eta3));
    
    if(deltam2){
        coeffs->rho81v2 =  (20022.-126451.*eta+236922.*eta2 - 138430.*eta3 + 21640.*eta4)/(18240.*(-1. + 6.*eta - 10.*eta2 + 4.*eta3));
        coeffs->rho83v2 = (20598.-131059.*eta+249018.*eta2-149950.*eta3+24520.*eta4)/(18240.*(-1.+6.*eta-10.*eta2+4.*eta3));        
        coeffs->rho85v2 = (4350.-28055.*eta+54642.*eta2-34598.*eta3+6056.*eta4)/(3648.*(-1.+6.*eta-10.*eta2+4.*eta3));        
        coeffs->rho87v2 = (23478.-154099.*eta+309498.*eta2-207550.*eta3+38920*eta4)/(18240.*(-1+6*eta-10*eta2+4*eta3));
    }
    
    return 1
 }   
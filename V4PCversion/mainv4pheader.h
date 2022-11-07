# include<math.h>

typedef struct WaveformCoefficients
{
  REAL8 delta22vh3;
  REAL8 delta22vh6;
  REAL8 delta22vh6S;
  REAL8 delta22v8;
  REAL8 delta22v8S;
  REAL8 delta22vh9;
  REAL8 delta22v5;
  REAL8 delta22v6;
  REAL8 delta22v6S;

  REAL8 rho22v2;
  REAL8 rho22v3;
  REAL8 rho22v3S;
  REAL8 rho22v4;
  REAL8 rho22v4S;
  REAL8 rho22v5;
  REAL8 rho22v5S;
  REAL8 rho22v6;
  REAL8 rho22v6S;
  REAL8 rho22v6l;
  REAL8 rho22v7;
  REAL8 rho22v7S;
  REAL8 rho22v8;
  REAL8 rho22v8S;
  REAL8 rho22v8l;
  REAL8 rho22v10;
  REAL8 rho22v10l;

  REAL8 delta21vh3;
  REAL8 delta21vh6;
  REAL8 delta21vh6S;
  REAL8 delta21vh7;
  REAL8 delta21vh7S;
  REAL8 delta21vh9;
  REAL8 delta21v5;
  REAL8 delta21v7;

  REAL8 rho21v1;
  REAL8 rho21v2;
  REAL8 rho21v2S;
  REAL8 rho21v3;
  REAL8 rho21v3S;
  REAL8 rho21v4;
  REAL8 rho21v4S;
  REAL8 rho21v5;
  REAL8 rho21v5S;
  REAL8 rho21v6;
  REAL8 rho21v6S;
  REAL8 rho21v6l;
  REAL8 rho21v7;
  REAL8 rho21v7S;
  REAL8 rho21v7l;
  REAL8 rho21v7lS;
  REAL8 rho21v8;
  REAL8 rho21v8l;
  REAL8 rho21v10;
  REAL8 rho21v10l;

  REAL8 f21v1;
  REAL8 f21v1S;
  REAL8 f21v3;
  REAL8 f21v3S;
  REAL8 f21v4;
  REAL8 f21v5;
  REAL8 f21v6;
  REAL8 f21v7c;

  REAL8 delta33vh3;
  REAL8 delta33vh6;
  REAL8 delta33vh6S;
  REAL8 delta33vh9;
  REAL8 delta33v5;
  REAL8 delta33v7;

  REAL8 rho33v2;
  REAL8 rho33v3;
  REAL8 rho33v4;
  REAL8 rho33v4S;
  REAL8 rho33v5;
  REAL8 rho33v5S;
  REAL8 rho33v6;
  REAL8 rho33v6S;
  REAL8 rho33v6l;
  REAL8 rho33v7;
  REAL8 rho33v7S;
  REAL8 rho33v8;
  REAL8 rho33v8l;
  REAL8 rho33v10;
  REAL8 rho33v10l;

  REAL8 f33v3;
  REAL8 f33v4;
  REAL8 f33v5;
  REAL8 f33v6;
  REAL8 f33v3S;
  REAL8 f33vh6;

  REAL8 delta32vh3;
  REAL8 delta32vh4;
  REAL8 delta32vh4S;
  REAL8 delta32vh6;
  REAL8 delta32vh6S;
  REAL8 delta32vh9;

  REAL8 rho32v;
  REAL8 rho32vS;
  REAL8 rho32v2;
  REAL8 rho32v2S;
  REAL8 rho32v3;
  REAL8 rho32v3S;
  REAL8 rho32v4;
  REAL8 rho32v4S;
  REAL8 rho32v5;
  REAL8 rho32v5S;
  REAL8 rho32v6;
  REAL8 rho32v6S;
  REAL8 rho32v6l;
  REAL8 rho32v8;
  REAL8 rho32v8l;

  REAL8 delta31vh3;
  REAL8 delta31vh6;
  REAL8 delta31vh6S;
  REAL8 delta31vh7;
  REAL8 delta31vh7S;
  REAL8 delta31vh9;
  REAL8 delta31v5;

  REAL8 rho31v2;
  REAL8 rho31v3;
  REAL8 rho31v4;
  REAL8 rho31v4S;
  REAL8 rho31v5;
  REAL8 rho31v5S;
  REAL8 rho31v6;
  REAL8 rho31v6S;
  REAL8 rho31v6l;
  REAL8 rho31v7;
  REAL8 rho31v7S;
  REAL8 rho31v8;
  REAL8 rho31v8l;

  REAL8 f31v3;
  REAL8 f31v3S;

  REAL8 delta44vh3;
  REAL8 delta44vh6;
  REAL8 delta44vh6S;
  REAL8 delta44v5;
  REAL8 delta44vh9;

  REAL8 rho44v2;
  REAL8 rho44v3;
  REAL8 rho44v3S;
  REAL8 rho44v4;
  REAL8 rho44v4S;
  REAL8 rho44v5;
  REAL8 rho44v5S;
  REAL8 rho44v6;
  REAL8 rho44v6S;
  REAL8 rho44v6l;
  REAL8 rho44v8;
  REAL8 rho44v8l;
  REAL8 rho44v10;
  REAL8 rho44v10l;

  REAL8 delta43vh3;
  REAL8 delta43vh4;
  REAL8 delta43vh4S;
  REAL8 delta43vh6;

  REAL8 rho43v;
  REAL8 rho43v2;
  REAL8 rho43v4;
  REAL8 rho43v4S;
  REAL8 rho43v5;
  REAL8 rho43v5S;
  REAL8 rho43v6;
  REAL8 rho43v6l;

  REAL8 f43v;
  REAL8 f43vS;

  REAL8 delta42vh3;
  REAL8 delta42vh6;
  REAL8 delta42vh6S;

  REAL8 rho42v2;
  REAL8 rho42v3;
  REAL8 rho42v3S;
  REAL8 rho42v4;
  REAL8 rho42v4S;
  REAL8 rho42v5;
  REAL8 rho42v5S;
  REAL8 rho42v6;
  REAL8 rho42v6S;
  REAL8 rho42v6l;

  REAL8 delta41vh3;
  REAL8 delta41vh4;
  REAL8 delta41vh4S;
  REAL8 delta41vh6;

  REAL8 rho41v;
  REAL8 rho41v2;
  REAL8 rho41v4;
  REAL8 rho41v4S;
  REAL8 rho41v5;
  REAL8 rho41v5S;
  REAL8 rho41v6;
  REAL8 rho41v6l;

  REAL8 f41v;
  REAL8 f41vS;

  REAL8 delta55vh3;
  REAL8 delta55vh6;
  REAL8 delta55vh9;

  REAL8 delta55v5;
  REAL8 rho55v2;
  REAL8 rho55v3;
  REAL8 rho55v3S;
  REAL8 rho55v4;
  REAL8 rho55v4S;
  REAL8 rho55v5;
  REAL8 rho55v5S;
  REAL8 rho55v6;
  REAL8 rho55v6l;
  REAL8 rho55v8;
  REAL8 rho55v8l;
  REAL8 rho55v10;
  REAL8 rho55v10l;
  REAL8 f55v3;
  REAL8 f55v4;
  REAL8 f55v5c;


  REAL8 delta54vh3;
  REAL8 delta54vh4;
  REAL8 delta54vh4S;
  REAL8 rho54v2;
  REAL8 rho54v3;
  REAL8 rho54v3S;
  REAL8 rho54v4;
  REAL8 rho54v4S;

  REAL8 delta53vh3;
  REAL8 rho53v2;
  REAL8 rho53v3;
  REAL8 rho53v3S;
  REAL8 rho53v4;
  REAL8 rho53v4S;
  REAL8 rho53v5;
  REAL8 rho53v5S;

  REAL8 delta52vh3;
  REAL8 delta52vh4;
  REAL8 delta52vh4S;
  REAL8 rho52v2;
  REAL8 rho52v3;
  REAL8 rho52v3S;
  REAL8 rho52v4;
  REAL8 rho52v4S;

  REAL8 delta51vh3;
  REAL8 rho51v2;
  REAL8 rho51v3;
  REAL8 rho51v3S;
  REAL8 rho51v4;
  REAL8 rho51v4S;
  REAL8 rho51v5;
  REAL8 rho51v5S;

  REAL8 delta66vh3;
  REAL8 rho66v2;
  REAL8 rho66v3;
  REAL8 rho66v3S;
  REAL8 rho66v4;
  REAL8 rho66v4S;

  REAL8 delta65vh3;
  REAL8 rho65v2;
  REAL8 rho65v3;
  REAL8 rho65v3S;

  REAL8 delta64vh3;
  REAL8 rho64v2;
  REAL8 rho64v3;
  REAL8 rho64v3S;
  REAL8 rho64v4;
  REAL8 rho64v4S;

  REAL8 delta63vh3;
  REAL8 rho63v2;
  REAL8 rho63v3;
  REAL8 rho63v3S;

  REAL8 delta62vh3;
  REAL8 rho62v2;
  REAL8 rho62v3;
  REAL8 rho62v3S;
  REAL8 rho62v4;
  REAL8 rho62v4S;

  REAL8 delta61vh3;
  REAL8 rho61v2;
  REAL8 rho61v3;
  REAL8 rho61v3S;

  REAL8 delta77vh3;
  REAL8 rho77v2;
  REAL8 rho77v3;
  REAL8 rho77v3S;

  REAL8 rho76v2;

  REAL8 delta75vh3;
  REAL8 rho75v2;
  REAL8 rho75v3;
  REAL8 rho75v3S;

  REAL8 rho74v2;

  REAL8 delta73vh3;
  REAL8 rho73v2;
  REAL8 rho73v3;
  REAL8 rho73v3S;

  REAL8 rho72v2;

  REAL8 delta71vh3;
  REAL8 rho71v2;
  REAL8 rho71v3;
  REAL8 rho71v3S;

  REAL8 rho88v2;
  REAL8 rho87v2;
  REAL8 rho86v2;
  REAL8 rho85v2;
  REAL8 rho84v2;
  REAL8 rho83v2;
  REAL8 rho82v2;
  REAL8 rho81v2;
}
FacWaveformCoeffs;
  // ========================= Begin Sid Hamiltonian Verbose 1 ====================== //

  #include<stdio.h>
  FILE *fp;
  fp = fopen("CI-HamiltonianTestMain.dat","rb");
  double values[15];
  fseek(fp,0,SEEK_SET);
  fread(&values,sizeof(double),15,fp);
  fclose(fp);
  for(int sidi = 0; sidi<3; sidi++){
        x->data[sidi] = (REAL8)values[2 + sidi];
        p->data[sidi] = (REAL8)values[5 + sidi];
        s1Vec->data[sidi] = (REAL8)values[8 + sidi];
        s2Vec->data[sidi] = (REAL8)values[11 + sidi];
  }
  double sidmass1 = values[0];
  double sidmass2 = values[1];
  eta = (REAL8) (sidmass1*sidmass2/((sidmass1 + sidmass2)*(sidmass1 + sidmass2)));
  tortoise = 2;
  XLALSimIMRSpinEOBCalculateSigmaKerr(sigmaKerr,sidmass1, sidmass2, s1Vec, s2Vec);
  XLALSimIMRSpinEOBCalculateSigmaStar(sigmaStar,sidmass1, sidmass2, s1Vec, s2Vec);
  coeffs->updateHCoeffs = 1;

  // ========================= End   Sid Hamiltonian Verbose 1 ====================== //

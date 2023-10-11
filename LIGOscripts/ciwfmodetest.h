        //================= Begin Sid Waveform Coefficients Intensive Debug =======================//

        #include<stdio.h>
        #include<time.h>
        FILE *fpmain, *fppert, *fpciwfmodes, *fpciwfmodespert;
        fpmain = fopen("CI-InputsWFModesMain.dat","wb");
        fppert = fopen("CI-InputsWFModesPert.dat","wb");
        fpciwfmodes = fopen("CI-WFModesMain.dat","wb");
        fpciwfmodespert = fopen("CI-WFModesPert.dat","wb");
for (int sidi = 0; sidi < 10000; sidi++){
        srand(time(NULL));
        REAL8 sidmass1 = 20*((REAL8)rand())/((REAL8)RAND_MAX);
        REAL8 sidmass2 = 10*((REAL8)rand())/((REAL8)RAND_MAX);
        REAL8 sideta   = (sidmass1*sidmass2/((sidmass1+sidmass2)*(sidmass1+sidmass2)));
        REAL8 sida     = 0.9*((REAL8)rand())/((REAL8)RAND_MAX) + 0.1;
        REAL8 sidchiS  = ((REAL8) (2*(rand() % 2) - 1))*(0.9*((REAL8)rand())/((REAL8)RAND_MAX) + 0.1);
        REAL8 sidchiA  = ((REAL8) (2*(rand() % 2) - 1))*(0.9*((REAL8)rand())/((REAL8)RAND_MAX) + 0.1);
        UINT4 sidSpinAlignedEOBversion = 4;
        FacWaveformCoeffs *sidhCoeffs = params.params->eobParams->hCoeffs;
        XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(sidhCoeffs, sidmass1, sidmass2, sideta, sida,sidchiS, sidchiA, sidSpinAlignedEOBversion);
        double indata[7] = {0.};
        indata[0] = (double)sidmass1;
        indata[1] = (double)sidmass2;
        indata[2] = (double)sideta;
        indata[3] = (double)sida;
        indata[4] = (double)sidchiS;
        indata[5] = (double)sidchiA;
        indata[6] = (double)(sizeof(FacWaveformCoeffs));
        fwrite(indata, sizeof(double), 7, fpmain);
        fwrite(sidhCoeffs,sizeof(FacWaveformCoeffs),1,fpciwfmodes);

        REAL8 pert_exponent = 1.e-14;
        REAL8 pert_sign = (REAL8) (2*(rand() % 2) - 1);
        REAL8 pert_mantissa = ( 3.*(REAL8)rand() )/((REAL8)RAND_MAX) + 1. ;
        REAL8 pert = 1. + pert_sign*pert_mantissa*pert_exponent;
        REAL8 sidmass1pert = pert*sidmass1;
        REAL8 sidmass2pert = pert*sidmass2;
        REAL8 sidetapert   = (sidmass1pert*sidmass2pert/((sidmass1pert+sidmass2pert)*(sidmass1pert+sidmass2pert)));
        REAL8 sidapert     = pert*sida;
        REAL8 sidchiSpert  = pert*sidchiS;
        REAL8 sidchiApert  = pert*sidchiS;
        FacWaveformCoeffs *sidhCoeffspert = sidhCoeffs;
        XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(sidhCoeffspert, sidmass1pert, sidmass2pert, sidetapert, sidapert,sidchiSpert, sidchiApert, sidSpinAlignedEOBversion);
        double indatapert[7] = {0.};
        indatapert[0] = (double)sidmass1pert;
        indatapert[1] = (double)sidmass2pert;
        indatapert[2] = (double)sidetapert;
        indatapert[3] = (double)sidapert;
        indatapert[4] = (double)sidchiSpert;
        indatapert[5] = (double)sidchiApert;
        indatapert[6] = (double)(sizeof(FacWaveformCoeffs));
        fwrite(indatapert, sizeof(double), 7, fppert);
        fwrite(sidhCoeffspert,sizeof(FacWaveformCoeffs),1,fpciwfmodespert);

        
}
fclose(fpmain);
fclose(fppert);
fclose(fpciwfmodes);
fclose(fpciwfmodespert);
exit(0); 
        //================= End Sid Waveform Coefficients Intensive Debug =======================//

        //================= Begin Sid Hamiltonian Intensive Debug =======================//

        #include<stdio.h>
        #include<time.h>
        FILE *fpmain, *fppert, *fpciham;
        fpmain = fopen("CI-InputsMain.dat","wb");
        fppert = fopen("CI-InputsPert.dat","wb");
        fpciham = fopen("CI-HamiltonianMain.dat","wb");
for (int sidi = 0; sidi < 100000; sidi++){
        srand(time(NULL));
        SpinEOBHCoeffs *sidHCoeffs;
        REAL8 sidKerrdata[3] = {0.};
        REAL8 sidStardata[3] = {0.};
        REAL8Vector     sidx, sidp, sids1, sids2, sidKerr, sidStar;
        REAL8           sidxdata[3], sidpdata[3], sids1Vecdata[3], sids2Vecdata[3];
        sidx.length = sidp.length = sids1.length = sids2.length = sidKerr.length = sidStar.length = 3;
        int sidtortoise = 2;
        REAL8 sidmass1 = 20*((REAL8)rand())/((REAL8)RAND_MAX);
        REAL8 sidmass2 = 10*((REAL8)rand())/((REAL8)RAND_MAX);
        REAL8 sideta = (sidmass1*sidmass2/((sidmass1+sidmass2)*(sidmass1+sidmass2)));
        sidxdata[0]         = 20*((REAL8)rand())/((REAL8)RAND_MAX) - 10;
        sidxdata[1]         = 20*((REAL8)rand())/((REAL8)RAND_MAX) - 10;
        sidxdata[2]         = 20*((REAL8)rand())/((REAL8)RAND_MAX) - 10;
        sidpdata[0]         = 20*((REAL8)rand())/((REAL8)RAND_MAX) - 10;
        sidpdata[1]         = 20*((REAL8)rand())/((REAL8)RAND_MAX) - 10;
        sidpdata[2]         = 20*((REAL8)rand())/((REAL8)RAND_MAX) - 10;
        sids1Vecdata[0]     = 2*((REAL8)rand())/((REAL8)RAND_MAX) - 1;
        sids1Vecdata[1]     = 2*((REAL8)rand())/((REAL8)RAND_MAX) - 1;
        sids1Vecdata[2]     = 2*((REAL8)rand())/((REAL8)RAND_MAX) - 1;
        sids2Vecdata[0]     = 2*((REAL8)rand())/((REAL8)RAND_MAX) - 1;
        sids2Vecdata[1]     = 2*((REAL8)rand())/((REAL8)RAND_MAX) - 1;
        sids2Vecdata[2]     = 2*((REAL8)rand())/((REAL8)RAND_MAX) - 1;

        sidx.data = sidxdata;
        sidp.data = sidpdata;
        sids1.data = sids1Vecdata;
        sids2.data = sids2Vecdata;
        sidKerr.data = sidKerrdata;
        sidStar.data = sidStardata;
        XLALSimIMRSpinEOBCalculateSigmaKerr(&sidKerr, sidmass1, sidmass2, &sids1, &sids2);
        XLALSimIMRSpinEOBCalculateSigmaStar(&sidStar, sidmass1, sidmass2, &sids1, &sids2);
	sidHCoeffs = params.params->seobCoeffs;
        sidHCoeffs->updateHCoeffs = 1;
	sidHCoeffs->SpinAlignedEOBversion = 4;
        REAL8 sidHam = XLALSimIMRSpinPrecEOBHamiltonian(sideta, &sidx, &sidp, &sids1, &sids2,
                                             &sidKerr, &sidStar, sidtortoise,sidHCoeffs);
        double indata[14] = {0.};
        double hamdata[2] = {0.};
        indata[0] = (double)sidmass1;
        indata[1] = (double)sidmass2;
        for (int indatai = 0; indatai<3; indatai++){
                indata[indatai + 2] = (double)sidxdata[indatai];
                indata[indatai + 5] = (double)sidpdata[indatai];
                indata[indatai + 8] = (double)sids1Vecdata[indatai];
                indata[indatai + 11] = (double)sids2Vecdata[indatai];
        }
        hamdata[0] = (double)sidHam;
        fwrite(indata, sizeof(double), 14, fpmain);

        REAL8 sidKerrdatapert[3] = {0.};
        REAL8 sidStardatapert[3] = {0.};
        REAL8 pert_exponent = 1.e-14;
        REAL8 pert_sign = (REAL8) (2*(rand() % 2) - 1);
        REAL8 pert_mantissa = ( 3.*(REAL8)rand() )/((REAL8)RAND_MAX) + 1. ;
        REAL8 pert = 1. + pert_sign*pert_mantissa*pert_exponent;
        REAL8Vector     sidxpert, sidppert, sids1pert, sids2pert, sidKerrpert, sidStarpert;
        REAL8           sidxdatapert[3], sidpdatapert[3], sids1Vecdatapert[3], sids2Vecdatapert[3];
        sidxpert.length = sidppert.length = sids1pert.length = sids2pert.length = sidKerrpert.length = sidStarpert.length = 3;
        REAL8 sidmass1pert = pert*sidmass1;
        REAL8 sidmass2pert = pert*sidmass2;
        REAL8 sidetapert = (sidmass1pert*sidmass2pert/((sidmass1pert+sidmass2pert)*(sidmass1pert+sidmass2pert)));
        sidxdatapert[0]         = sidxdata[0]*pert;
        sidxdatapert[1]         = sidxdata[1]*pert;
        sidxdatapert[2]         = sidxdata[2]*pert;
        sidpdatapert[0]         = sidpdata[0]*pert;
        sidpdatapert[1]         = sidpdata[1]*pert;
        sidpdatapert[2]         = sidpdata[2]*pert;
        sids1Vecdatapert[0]     = sids1Vecdata[0]*pert;
        sids1Vecdatapert[1]     = sids1Vecdata[1]*pert;
        sids1Vecdatapert[2]     = sids1Vecdata[2]*pert;
        sids2Vecdatapert[0]     = sids2Vecdata[0]*pert;
        sids2Vecdatapert[1]     = sids2Vecdata[1]*pert;
        sids2Vecdatapert[2]     = sids2Vecdata[2]*pert;

        sidxpert.data = sidxdatapert;
        sidppert.data = sidpdatapert;
        sids1pert.data = sids1Vecdatapert;
        sids2pert.data = sids2Vecdatapert;
        sidKerrpert.data = sidKerrdatapert;
        sidStarpert.data = sidStardatapert;
        XLALSimIMRSpinEOBCalculateSigmaKerr(&sidKerrpert, sidmass1pert, sidmass2pert, &sids1pert, &sids2pert);
        XLALSimIMRSpinEOBCalculateSigmaStar(&sidStarpert, sidmass1pert, sidmass2pert, &sids1pert, &sids2pert);
        sidHCoeffs->updateHCoeffs = 1;
        REAL8 sidHampert = XLALSimIMRSpinPrecEOBHamiltonian(sidetapert, &sidxpert, &sidppert, &sids1pert, &sids2pert,
                                             &sidKerrpert, &sidStarpert, sidtortoise, sidHCoeffs);
        indata[0] = (double)sidmass1pert;
        indata[1] = (double)sidmass2pert;
        for (int indatai = 0; indatai<3; indatai++){
                indata[indatai + 2] = (double)sidxdatapert[indatai];
                indata[indatai + 5] = (double)sidpdatapert[indatai];
                indata[indatai + 8] = (double)sids1Vecdatapert[indatai];
                indata[indatai + 11] =(double) sids2Vecdatapert[indatai];
        }
        hamdata[1] = (double) sidHampert;
        fwrite(hamdata, sizeof(double),2,fpciham);
        fwrite(indata, sizeof(double), 14, fppert);

}
fclose(fpmain);
fclose(fppert);
fclose(fpciham);
exit(0); 

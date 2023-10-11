        //================= Begin Sid Hamiltonian Circular Omega Intensive Debug =======================//

        //declare and open files
        #include<stdio.h>
        FILE *fpmain, *fppert, *fpciomegacirc;
        fpmain = fopen("CI-InputsMain.dat","rb");
        fppert = fopen("CI-InputsPert.dat","rb");
        fpciomegacirc = fopen("CI-OmegaCircular.dat","wb");
        
        //declare the gsl and lalsuite structures
        SpinEOBParams *sidparams = funcParams
        double indata[14] = {0.};
        double indatapert[14] = {0.};
        REAL8 sidvalues[12] = {0.};
        REAL8 sidvaluespert[12] = {0.};
        REAL8 omegacircvalues[2] = {0.};

        // start the loop
for (int sidi = 0; sidi < 100000; sidi++){
        fread(indata, sizeof(double), 14, fpmain);

        for (int indatai = 0; indatai<12; indatai++){
                sidvalues[indatai] = (REAL8)indata[indatai + 2];
        }

        // populate lalsuite parameter structure
        
        sidparams->seobCoeffs->updateHCoeffs = 1;
        sidparams->seobCoeffs->SpinAlignedEOBversion = 4;
        sidparams->eobParams->m1 = (REAL8)indata[0];
        REAL8 sidmass1 = (REAL8)indata[0];
        sidparams->eobParams->m2 = (REAL8)indata[1];
        REAL8 sidmass2 = (REAL8)indata[1];
        sidderivparams.params->eobParams->eta = (sidmass1*sidmass2/((sidmass1+sidmass2)*(sidmass1+sidmass2)));
        
        // start computing the partials
        REAL8 omegacircmain = XLALSimIMRSpinPrecEOBCalcOmega( sidvalues, sidparams );
        // do the same for the perturbed values

        fread(indatapert, sizeof(double), 14, fppert);

        for (int indatai = 0; indatai<12; indatai++){
                sidvaluespert[indatai] = (REAL8)indatapert[indatai + 2];
        }

        // populate lalsuite parameter structure

        sidparams->seobCoeffs->updateHCoeffs = 1;
        sidparams->seobCoeffs->SpinAlignedEOBversion = 4;
        sidparams->eobParams->m1 = (REAL8)indatapert[0];
        REAL8 sidmass1pert = (REAL8)indatapert[0];
        sidparams->eobParams->m2 = (REAL8)indatapert[1];
        REAL8 sidmass2pert = (REAL8)indatapert[1];
        sidderivparams.params->eobParams->eta = (sidmass1pert*sidmass2pert/((sidmass1pert+sidmass2pert)*(sidmass1pert+sidmass2pert)));
        
        // start computing the partials
        REAL8 omegacircpert = XLALSimIMRSpinPrecEOBCalcOmega( sidvaluespert, sidparams );
        
        omegacircvalues[0] = (double)omegacircmain;
        omegacircvalues[1] = (double)omegacircpert;

        fwrite(omegacircvalues,sizeof(double),2,fpciomegacirc);

}
fclose(fpmain);
fclose(fppert);
fclose(fpciomegacirc);
int siddebug = 1;
if (siddebug){
exit(0);
}

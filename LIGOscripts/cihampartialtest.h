        //================= Begin Sid Hamiltonian Partial Derivative Intensive Debug =======================//

        //declare and open files
        #include<stdio.h>
        FILE *fpmain, *fppert, *fpcihampartial, *fpcihampartialpert;
        fpmain = fopen("CI-InputsMain.dat","rb");
        fppert = fopen("CI-InputsPert.dat","rb");
        fpcihampartial = fopen("CI-HamiltonianPartialsMain.dat","wb");
        fpcihampartialpert = fopen("CI-HamiltonianPartialsPert.dat","wb");
        
        //declare the gsl and lalsuite structures
        gsl_function    sidF;
        HcapDerivParams sidderivparams;
        sidF.function = &GSLSpinPrecHamiltonianWrapper;
        F.params = &sidderivparams;
        REAL8 v4pstepsize = 2.0e-3;
        REAL8 sidabserror = 0.;
        double indata[14] = {0.};
        double indatapert[14] = {0.};
        REAL8 sidvalues[12] = {0.};

        // start the loop
for (int sidi = 0; sidi < 100000; sidi++){
        fread(&indata, sizeof(double), 14, fpmain);

        for (int indatai = 0; indatai<12; indatai++){
                sidvalues[indatai] = (REAL8)indata[indatai + 2];
        }

        // populate lalsuite derivative structure
        sidderivparams.values = sidvalues;
        sidderivparams.params = (SpinEOBParams *) funcParams;
        sidderivparams.params->tortoise = 2;
        sidderivparams.params->seobCoeffs->updateHCoeffs = 1;
        sidderivparams.params->seobCoeffs->SpinAlignedEOBversion = 4;
        sidderivparams.params->eobParams->m1 = (REAL8)indata[0];
        REAL8 sidmass1 = (REAL8)indata[0];
        sidderivparams.params->eobParams->m2 = (REAL8)indata[1];
        REAL8 sidmass2 = (REAL8)indata[1];
        sidderivparams.params->eobParams->eta = (sidmass1*sidmass2/((sidmass1+sidmass2)*(sidmass1+sidmass2)));
        
        // start computing the partials
        REAL8 HamiltonianPartials[12] = {0.};
        double hampartials[12];
        int cigslStatus = 0.;
        for (int vari = 0; vari < 12; vari++){
            sidderivparams.varyParam = vari;
            REAL8 actual_stepsize = v4pstepsize;
            if (vari > 5 && vari < 9){
                actual_stepsize *= ((sidmass1/(sidmass1 + sidmass2))*(sidmass1/(sidmass1 + sidmass2)));
            }
            if (vari > 8 && vari < 12){
                actual_stepsize *= ((sidmass2/(sidmass1 + sidmass2))*(sidmass2/(sidmass1 + sidmass2)));
            }
            XLAL_CALLGSL(cigslStatus = gsl_deriv_central(&sidF, sidvalues[vari], actual_stepsize, &HamiltonianPartials[vari], &sidabserror));
            if (cigslStatus != GSL_SUCCESS){
                printf("it's okay\n");
            }
            hampartials[vari] = (double) HamiltonianPartials[vari]; 
        }
        // write the partials into the file
        
        fwrite(hampartials,sizeof(double),12,fpcihampartial);

        // do the same for the perturbed values

        fread(&indatapert, sizeof(double), 14, fppert);
        sidderivparams.params->tortoise = 2;
        sidderivparams.params->seobCoeffs->updateHCoeffs = 1;
        sidderivparams.params->seobCoeffs->SpinAlignedEOBversion = 4;
        REAL8 sidvaluespert[12];
        sidderivparams.params->eobParams->m1 = (REAL8)indatapert[0];
         REAL8 sidmass1pert = (REAL8)indatapert[0];
        sidderivparams.params->eobParams->m2 = (REAL8)indatapert[1];
        REAL8 sidmass2pert = (REAL8)indatapert[1];
        sidderivparams.params->eobParams->eta = (sidmass1pert*sidmass2pert/((sidmass1pert+sidmass2pert)*(sidmass1pert+sidmass2pert)));
        for (int indatai = 0; indatai<12; indatai++){
                sidvaluespert[indatai] = (REAL8)indatapert[indatai + 2];
        }
        sidderivparams.values = sidvaluespert;
        REAL8 HamiltonianPartialsPert[12] = {0.};
        double hampartialspert[12] = {0.};
        for (int vari = 0; vari < 12; vari++){
            sidderivparams.varyParam = vari;
            REAL8 actual_stepsize = v4pstepsize;
            if (vari > 5 && vari < 9){
                actual_stepsize *= ((sidmass1pert/(sidmass1pert + sidmass2pert))*(sidmass1pert/(sidmass1pert + sidmass2pert)));
            }
            if (vari > 8 && vari < 12){
                actual_stepsize *= ((sidmass2pert/(sidmass1pert + sidmass2pert))*(sidmass2pert/(sidmass1pert + sidmass2pert)));
            }
        
            XLAL_CALLGSL(cigslStatus = gsl_deriv_central(&sidF, sidvaluespert[vari], actual_stepsize, &HamiltonianPartialsPert[vari], &sidabserror));
            hampartialspert[vari] = (double)HamiltonianPartialsPert[vari];
        
        }

        fwrite(hampartialspert,sizeof(double),12,fpcihampartialpert);

}
fclose(fpmain);
fclose(fppert);
fclose(fpcihampartial);
fclose(fpcihampartialpert);
int siddebug = 1;
if (siddebug){
exit(0);
}

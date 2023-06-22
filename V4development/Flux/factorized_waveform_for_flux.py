import numpy as np
import Flux.auxiliaryfunctions as aux

def factorial(n):
    if n == 1:
        return 1
    return n*factorial(n-1)

def compute_hFlm_for_flux(m1, m2, tortoise, x, y, z, p1, p2, p3, S1x, S1y, S1z, S2x, S2y, S2z, l, m, Hreal, v, rho,f,a,chia,chis):
    debugwaveform = 0
    eta = m1*m2/(m1+m2)/(m1+m2)
    r = sqrt(x*x + y*y + z*z)
    v2 = v*v
    Omega = v2*v
    vPhi = r*Omega*cbrt(vPhiNonKeplerian(m1,m2,tortoise,x,y,z,p1,p2,p3,S1x,S1y,S1z,S2x,S2y,S2z))
    L1 = y*p3 - z*p2
    L2 = z*p1 - x*p3
    L3 = x*p2 - y*p1
    pp = sqrt(L1*L1 + L2*L2 + L3*L3)
    EMgamma = 0.577215664901532860606512090082402431
    lMax = 8
    hlm = np.zeros([lMax+1,lMax+1])
    for l in range(2,lMax+1):

        for m in range(1,l+1):
            epsilon = (l+m)%2
            eulerlogxabs = EMgamma + log(2.0* m * v)
            Vl_Phi = np.pow(vPhi,l+epsilon)
            Ylminuse_minusm = AbsSphericalHarmonicAtPiOver2(l-epsilon,-m)
            prefixstatus = Newtonian_flux_prefix(&prefix,m1,m2,l,m)
            prefixflux = sqrt(creal(prefix)*creal(prefix) + cimag(prefix)*cimag(prefix))
            hlm = prefix*Vl_Phi*Ylminuse_minusm
            
            if (epsilon == 0):
                Slm = (Hreal * Hreal - 1.) / (2. * eta) + 1.
            else:
                Slm = v * pp
            
            hathatk = m*Hreal*Omega
            hathatksq4 = 4. * hathatk * hathatk
            hathatk4pi = 4. * np.pi * hathatk
            exppihathatk4m1 = 1 - exp( -hathatk4pi )
            lfactorialinv = 1./factorial(l)
            Tlmprefac = sqrt(hathatk4pi / (1. - exp(-hathatk4pi)))*lfactorialinv
            Tlmprodfac = 1.
            for i in range(1,l+1):
                Tlmprodfac *= (hathatksq4 + i * i)    
            Tlm = Tlmprefac*sqrt(Tlmprodfac)
            rholm = 1.
            vp = 1.
            auxflm = v*(f[l,m,1] + v2*f[l,m,3])
            for p in range(1,10):
                vp *=v
                rholm += vp*rho[l,m,p,0]
                if p in [6,8,10]:
                    rholm += vp*rho[l,m,p,1]
            
            switch (l) {
			case 2:
				switch (abs(m)) {
				case 2:
					rholm = 1. + v2 * (hCoeffs->rho22v2 + v * (hCoeffs->rho22v3
								     + v * (hCoeffs->rho22v4
					     + v * (hCoeffs->rho22v5 + v * (hCoeffs->rho22v6
									    + hCoeffs->rho22v6l * eulerlogxabs + v * (hCoeffs->rho22v7
														      + v * (hCoeffs->rho22v8 + hCoeffs->rho22v8l * eulerlogxabs
															     + (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs) * v2)))))))
					break
				case 1:
					{
						rholm = 1. + v * (hCoeffs->rho21v1
								  + v * (hCoeffs->rho21v2 + v * (hCoeffs->rho21v3 + v * (hCoeffs->rho21v4
															 + v * (hCoeffs->rho21v5 + v * (hCoeffs->rho21v6 + hCoeffs->rho21v6l * eulerlogxabs
																			+ v * (hCoeffs->rho21v7 + hCoeffs->rho21v7l * eulerlogxabs
																			       + v * (hCoeffs->rho21v8 + hCoeffs->rho21v8l * eulerlogxabs
																				      + (hCoeffs->rho21v10 + hCoeffs->rho21v10l * eulerlogxabs) * v2))))))))
						auxflm = v * hCoeffs->f21v1 + v2 * v * hCoeffs->f21v3
					}
					break
				default:
					exit(0)
					break
				}
				break
			case 3:
				switch (m) {
				case 3:
					rholm = 1. + v2 * (hCoeffs->rho33v2 + v * (hCoeffs->rho33v3 + v * (hCoeffs->rho33v4
													   + v * (hCoeffs->rho33v5 + v * (hCoeffs->rho33v6 + hCoeffs->rho33v6l * eulerlogxabs
																	  + v * (hCoeffs->rho33v7 + (hCoeffs->rho33v8 + hCoeffs->rho33v8l * eulerlogxabs) * v))))))
					auxflm = v * v2 * hCoeffs->f33v3
					break
				case 2:
					rholm = 1. + v * (hCoeffs->rho32v1
							  + v * (hCoeffs->rho32v2 + v * (hCoeffs->rho32v3 + v * (hCoeffs->rho32v4 + v * (hCoeffs->rho32v5
																	 + v * (hCoeffs->rho32v6 + hCoeffs->rho32v6l * eulerlogxabs
																		+ (hCoeffs->rho32v8 + hCoeffs->rho32v8l * eulerlogxabs) * v2))))))
					break
				case 1:
					rholm = 1. + v2 * (hCoeffs->rho31v2 + v * (hCoeffs->rho31v3 + v * (hCoeffs->rho31v4
													   + v * (hCoeffs->rho31v5 + v * (hCoeffs->rho31v6 + hCoeffs->rho31v6l * eulerlogxabs
																	  + v * (hCoeffs->rho31v7 + (hCoeffs->rho31v8 + hCoeffs->rho31v8l * eulerlogxabs) * v))))))
					auxflm = v * v2 * hCoeffs->f31v3
					break
				default:
					exit(0)
					break
				}
				break
			case 4:
				switch (m) {
				case 4:

					rholm = 1. + v2 * (hCoeffs->rho44v2
					     + v * (hCoeffs->rho44v3 + v * (hCoeffs->rho44v4
						 + v * (hCoeffs->rho44v5 + (hCoeffs->rho44v6
						+ hCoeffs->rho44v6l * eulerlogxabs) * v))))
					break
				case 3:
					rholm = 1. + v * (hCoeffs->rho43v1
							  + v * (hCoeffs->rho43v2
					    + v2 * (hCoeffs->rho43v4 + v * (hCoeffs->rho43v5
									    + (hCoeffs->rho43v6 + hCoeffs->rho43v6l * eulerlogxabs) * v))))
					auxflm = v * hCoeffs->f43v1
					break
				case 2:
					rholm = 1. + v2 * (hCoeffs->rho42v2
							   + v * (hCoeffs->rho42v3 + v * (hCoeffs->rho42v4 + v * (hCoeffs->rho42v5
														  + (hCoeffs->rho42v6 + hCoeffs->rho42v6l * eulerlogxabs) * v))))
					break
				case 1:
					rholm = 1. + v * (hCoeffs->rho41v1
							  + v * (hCoeffs->rho41v2
					    + v2 * (hCoeffs->rho41v4 + v * (hCoeffs->rho41v5
									    + (hCoeffs->rho41v6 + hCoeffs->rho41v6l * eulerlogxabs) * v))))
					auxflm = v * hCoeffs->f41v1
					break
				default:
					exit(0)
					break
				}
				break
			case 5:
				switch (m) {
				case 5:
					rholm = 1. + v2 * (hCoeffs->rho55v2
					     + v * (hCoeffs->rho55v3 + v * (hCoeffs->rho55v4
					 + v * (hCoeffs->rho55v5 + hCoeffs->rho55v6 * v))))
					break
				case 4:
					rholm = 1. + v2 * (hCoeffs->rho54v2 + v * (hCoeffs->rho54v3
								   + hCoeffs->rho54v4 * v))
					break
				case 3:
					rholm = 1. + v2 * (hCoeffs->rho53v2
							   + v * (hCoeffs->rho53v3 + v * (hCoeffs->rho53v4 + hCoeffs->rho53v5 * v)))
					break
				case 2:
					rholm = 1. + v2 * (hCoeffs->rho52v2 + v * (hCoeffs->rho52v3
								   + hCoeffs->rho52v4 * v))
					break
				case 1:
					rholm = 1. + v2 * (hCoeffs->rho51v2
							   + v * (hCoeffs->rho51v3 + v * (hCoeffs->rho51v4 + hCoeffs->rho51v5 * v)))
					break
				default:
					exit(0)
					break
				}
				break
			case 6:
				switch (m) {
				case 6:
					rholm = 1. + v2 * (hCoeffs->rho66v2 + v * (hCoeffs->rho66v3
								   + hCoeffs->rho66v4 * v))
					break
				case 5:
					rholm = 1. + v2 * (hCoeffs->rho65v2 + hCoeffs->rho65v3 * v)
					break
				case 4:
					rholm = 1. + v2 * (hCoeffs->rho64v2 + v * (hCoeffs->rho64v3
								   + hCoeffs->rho64v4 * v))
					break
				case 3:
					rholm = 1. + v2 * (hCoeffs->rho63v2 + hCoeffs->rho63v3 * v)
					break
				case 2:
					rholm = 1. + v2 * (hCoeffs->rho62v2 + v * (hCoeffs->rho62v3
								   + hCoeffs->rho62v4 * v))
					break
				case 1:
					rholm = 1. + v2 * (hCoeffs->rho61v2 + hCoeffs->rho61v3 * v)
					break
				default:
					exit(0)
					break
				}
				break
			case 7:
				switch (m) {
				case 7:
					rholm = 1. + v2 * (hCoeffs->rho77v2 + hCoeffs->rho77v3 * v)
					break
				case 6:
					rholm = 1. + hCoeffs->rho76v2 * v2
					break
				case 5:
					rholm = 1. + v2 * (hCoeffs->rho75v2 + hCoeffs->rho75v3 * v)
					break
				case 4:
					rholm = 1. + hCoeffs->rho74v2 * v2
					break
				case 3:
					rholm = 1. + v2 * (hCoeffs->rho73v2 + hCoeffs->rho73v3 * v)
					break
				case 2:
					rholm = 1. + hCoeffs->rho72v2 * v2
					break
				case 1:
					rholm = 1. + v2 * (hCoeffs->rho71v2 + hCoeffs->rho71v3 * v)
					break
				default:
					exit(0)
					break
				}
				break
			case 8:
				switch (m) {
				case 8:
					rholm = 1. + hCoeffs->rho88v2 * v2
					break
				case 7:
					rholm = 1. + hCoeffs->rho87v2 * v2
					break
				case 6:
					rholm = 1. + hCoeffs->rho86v2 * v2
					break
				case 5:
					rholm = 1. + hCoeffs->rho85v2 * v2
					break
				case 4:
					rholm = 1. + hCoeffs->rho84v2 * v2
					break
				case 3:
					rholm = 1. + hCoeffs->rho83v2 * v2
					break
				case 2:
					rholm = 1. + hCoeffs->rho82v2 * v2
					break
				case 1:
					rholm = 1. + hCoeffs->rho81v2 * v2
					break
				default:
					exit(0)
					break
				}
				break
			default:
				exit(0)
				break
			}

            rholmPwrl = 1.0
            /*
             * In the equal-mass odd m case, there is no contribution from
             * nonspin terms,  and the only contribution comes from the auxflm
			 * term that is proportional to chiA (asymmetric spins). In this
			 * case, we must ignore the nonspin terms directly, since the leading
			 * term defined by CalculateThisMultipolePrefix in
			 * LALSimIMREOBNewtonianMultipole.c is not zero (see comments there).
			 */
            if (eta == 0.25 && m % 2) {
                rholmPwrl = auxflm
            }
            else {
                i = l
                while (i--) {
                    rholmPwrl *= rholm
                }
                rholmPwrl += auxflm
            }
            
            if (debugwaveform){
                printf("\n Checking l,m \n l = %d , m = %d \n",l,m)
                printf("\n checking functions \n")
                printf("eta = %.4e",eta)
                printf("Slm = %.4e \n",Slm)
                printf("lfactorialinv = %.4e \n",lfactorialinv)
                printf("Vl_Phi = %.4e \n",Vl_Phi)
                printf("Ylminuse_minusm = %.4e \n",Ylminuse_minusm)
                printf("Newtonian prefix = %.4e \n",prefix)
                printf("hathatk = %.4e \n", hathatk)
                printf("Tlmprefac = %.4e \n",Tlmprefac)
                printf("Tlmprodfac = %.4e \n",Tlmprodfac)
                printf("Tlm = %.4e \n",Tlm)
                printf("rholmPwrl = %.4e \n",rholmPwrl)
                exit(1)
            }
            
            *hlm *= Slm*Tlm*rholmPwrl
        }
    }
    return 1
}
    
    eta = m1*m2/(m1+m2)/(m1+m2)
    
    epsilon = (l+m)%2
    
    r = np.linalg.norm(q)
    
    rholmpowl = aux.rholmpowl(m1,m2,l,m,chiA,chiS,v)
    
    Se_eff = aux.Se_eff(l,m,Hreal,v,q,p,eta)
    
    lfactorialinv = 1/np.math.factorial(l)
    
    Omega = np.power(v,3)
    
    hathatk = m*Hreal*Omega
    
    pihathatk4 = 4*np.pi*hathatk
    
    Tlmprefac = np.sqrt( pihathatk4 / ( 1 - np.exp( -pihathatk4 ) ) )
    
    Tlmprodfac = aux.Tlmprodfac(l,hathatk)
    
    T_lm = lfactorialinv*Tlmprefac*np.sqrt(Tlmprodfac)
    
    vPhi = r*Omega*np.cbrt(aux.vPhiNonKeplerian(m1,m2,tortoise,q,p,S1,S2))
    
    Vl_Phi = np.power(vPhi,l+epsilon)
    
    Ylminuse_minusm = aux.AbsSphericalHarmonicAtPiOver2(l-epsilon,-m)
    
    c_lpluse = aux.Newtonian_c(m1,m2,l,m)
    
    ne_lm = aux.Newtonian_n(m1,m2,l,m)
    
    hNe_lm = eta*ne_lm*c_lpluse*Vl_Phi*Ylminuse_minusm
    
    hF_lm = hNe_lm*Se_eff*T_lm*rholmpowl
    return hF_lm
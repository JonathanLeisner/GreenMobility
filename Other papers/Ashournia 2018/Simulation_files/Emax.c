//------------------------------------------------------------------------------
// Emax.c
//
// By Damoun Ashournia, University of Oxford
//
//------------------------------------------------------------------------------
// Functions for computing expected values under static and forward looking
// expectations with unemployment probabilities fixed over time.
//------------------------------------------------------------------------------

#define iEpsEmax	(iAge - 1) + iDraw * (nAge - 1) + (iSec-1) * nDraws * (nAge - 1)
#define iPI1		iReg + iAge * nReg + iLagSec * (nAge - 1) * nReg + iFemale * (nSector+1) * (nAge - 1) * nReg + iEduc * CatGender * (nSector+1) * (nAge - 1) * nReg + iElig * CatEduc * CatGender * (nSector+1) * (nAge - 1) * nReg + iType * CatElig * CatEduc * CatGender * (nSector+1) * (nAge - 1) * nReg
#define iPIEmax		iReg + (iAge - 1) * nReg + iLagSec * (nAge - 1) * nReg + iFemale * (nSector+1) * (nAge - 1) * nReg + iEduc * CatGender * (nSector+1) * (nAge - 1) * nReg + iElig * CatEduc * CatGender * (nSector+1) * (nAge - 1) * nReg + iType * CatElig * CatEduc * CatGender * (nSector+1) * (nAge - 1) * nReg
#define iRsq		(iAge - 1) + iLagSec * (nAge - 1) + iFemale * (nSector+1) * (nAge - 1) + iEduc * CatGender * (nSector+1) * (nAge - 1) + iElig * CatEduc * CatGender * (nSector+1) * (nAge - 1) + iType * CatElig * CatEduc * CatGender * (nSector+1) * (nAge - 1)
#define iDeltaSE          0 + iAge * (LastYr_sim - FirstYr + 1) + female * nAge * (LastYr_sim - FirstYr + 1) + educ * CatGender * nAge * (LastYr_sim - FirstYr + 1) + iLagSec * CatEduc * CatGender * nAge * (LastYr_sim - FirstYr + 1)
#define iDeltaRE    iYear-1 + iAge * (LastYr_sim - FirstYr + 1) + female * nAge * (LastYr_sim - FirstYr + 1) + educ * CatGender * nAge * (LastYr_sim - FirstYr + 1) + iLagSec * CatEduc * CatGender * nAge * (LastYr_sim - FirstYr + 1)

//-------------------------------------------------------------------------------------
// This function computes the EMAX coefficients when agents have static expectations.
//-------------------------------------------------------------------------------------
void Emax_SE(int tid, double *param, int iFemale, int iEduc, int iElig, int iType, double *Low_RSk, double *Up_RSk, double Low_Wage, double Up_Wage, double *PI_SE, double *R_sq_SE)
{
	//---------------------------------------------------------------------------------
	// Initialize
	//---------------------------------------------------------------------------------
	double *beta, *sigma, *gamma, *kappa, *xi, *nu, *lambda, *phi, *eps, *rsk, *PI_Coef, *Emax, *w, *V, *Cost, SumExpV, MeanEmax, VM, benefit, lastWage, Emax0, exper, experChg, ExperTomorrow;
	int female, educ, elig, typ;
	int iAge, iLagSec, iSec, iN, iDraw, iReg;

	// Emax regression
	double *X, *Y, *coef, SST, SSE;


	// GSL Random Number Generator declarations
	double unif;
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng * s;


	//---------------------------------------------------------------------------------
	// Allocate memory
	//---------------------------------------------------------------------------------
	eps     = calloc((nAge-1)*nDraws*nSector,sizeof(double));
	rsk     = calloc(nSector,sizeof(double));
	PI_Coef = calloc(nReg,sizeof(double));
	Emax    = calloc(nSector,sizeof(double));
	w       = calloc(nSector,sizeof(double));
	V       = calloc(nSector,sizeof(double));
	Cost    = calloc(nSector,sizeof(double));
	beta    = calloc(nSector*nParamBeta,sizeof(double));
	sigma   = calloc(nSector,sizeof(double));
	gamma   = calloc(nSector+1,sizeof(double));
	kappa   = calloc(4,sizeof(double));
	xi      = calloc(2*nSector,sizeof(double));
	nu      = calloc(1,sizeof(double));
	lambda  = calloc(nSector,sizeof(double));
	phi     = calloc(4,sizeof(double));

	// Emax regressions
	X    = calloc(Intp*nReg,sizeof(double));
	Y    = calloc(Intp,sizeof(double));
	coef = calloc(Intp,sizeof(double));

	// Initialize Emax0
	Emax0 = 0.0;


	//--------------------------------------------------------------------
	// Unpack parameters
	//--------------------------------------------------------------------
	unpackPars(param, beta, sigma, gamma, xi, kappa, nu, lambda, phi);


	//---------------------------------------------------------------------------------
	// Random numbers for Monte Carlo integration over wage shocks. Draw only nDraws
	// for each age to save computational time.
	//---------------------------------------------------------------------------------
	T = gsl_rng_mt19937; // Algorithm
	r = gsl_rng_alloc(T); // Allocate memory
	gsl_rng_set(r,12345678); // Set seed

	for (iAge = 1; iAge < nAge; iAge++)
	{
		for (iDraw = 0; iDraw < nDraws; iDraw++)
		{
			for (iSec = 1; iSec < nSector + 1; iSec++)
			{
				eps[iEpsEmax] = gsl_ran_gaussian(r, sigma[iSec-1]);
			}
		}
	}

	// Free allocated memory for random number generator
	gsl_rng_free (r);

	// Random number generator for experience and tenure draws
	T = gsl_rng_mt19937; // Algorithm
	s = gsl_rng_alloc(T); // Allocate memory
	gsl_rng_set(s,987656); // Set seed

	//---------------------------------------------------------------------------------
	// Compute Emax_SE
	//---------------------------------------------------------------------------------
	for (iAge = nAge - 1; iAge > 0; iAge--)
	{
		for (iLagSec = 0; iLagSec < nSector + 1; iLagSec++)
		{
			for (iN = 0; iN < Intp; iN++)
			{
				//-------------------------------------------------------------
				// Random number generation - skill returns, exper and tenure
				//-------------------------------------------------------------
				exper  = gsl_rng_uniform_int(s, iAge + 5) + 2.0; // Draw exper from 2 to iAge + 5

				for (iSec = 1; iSec < nSector + 1; iSec++)
				{
					unif = gsl_rng_uniform_pos(s);
					rsk[iSec-1] = Low_RSk[iSec-1] + (Up_RSk[iSec-1] - Low_RSk[iSec-1])*unif;
				}

				unif = gsl_rng_uniform_pos(s);
				lastWage = Low_Wage + (Up_Wage - Low_Wage) * unif;


				//-------------------------------------------------------------
				// Compute Emax for iAge < LastAge
				//-------------------------------------------------------------
				for (iSec = 1; iSec < nSector + 1; iSec++)
				{
					if (iLagSec == iSec)
					{
						ExperTomorrow = exper + 1;
					}
					else
					{
						ExperTomorrow = gamma[iSec]*exper;
					}

					if (iAge < nAge - 1)
					{
						for (iReg = 0; iReg < nReg; iReg++)
						{
							PI_Coef[iReg] = PI_SE[iPI1];
						}

						Emax[iSec-1] = Emax_Hat(PI_Coef,rsk,ExperTomorrow,lastWage);
						Emax0 = Emax_Hat(PI_Coef,rsk,gamma[0]*exper,lastWage);
					}
					else
					{
						Emax[iSec-1] = 0.0;
						Emax0 = 0.0;
					}
				}


				//-------------------------------------------------------------
				// Wage and and value function
				//-------------------------------------------------------------
				female = iFemale;
				educ = iEduc;
				elig = iElig;
				typ = iType + 1;

				for (iSec = 1; iSec < nSector + 1; iSec++)
				{
					if (iLagSec == iSec)
					{
						if (typ == 1)
						{
							w[iSec-1] = rsk[iSec-1] * exp( beta[0+(iSec-1)*nParamBeta]*female + beta[1+(iSec-1)*nParamBeta]*educ
											         + beta[2+(iSec-1)*nParamBeta]*exper  + beta[3+(iSec-1)*nParamBeta]*exper*exper );
						}
						else if (typ == 2)
						{
							w[iSec-1] = rsk[iSec-1] * exp( beta[0+(iSec-1)*nParamBeta]*female + beta[1+(iSec-1)*nParamBeta]*educ
											         + beta[2+(iSec-1)*nParamBeta]*exper  + beta[3+(iSec-1)*nParamBeta]*exper*exper + lambda[iSec-1] );
						}

						Cost[iSec-1] = 0.0;
					}
					else if (iLagSec == 0)
					{
						if (typ == 1)
						{
							w[iSec-1] = rsk[iSec-1] * exp( beta[0+(iSec-1)*nParamBeta]*female + beta[1+(iSec-1)*nParamBeta]*educ
											         + beta[2+(iSec-1)*nParamBeta]*exper  + beta[3+(iSec-1)*nParamBeta]*exper*exper );
						}
						else if (typ == 2)
						{
							w[iSec-1] = rsk[iSec-1] * exp( beta[0+(iSec-1)*nParamBeta]*female + beta[1+(iSec-1)*nParamBeta]*educ
											         + beta[2+(iSec-1)*nParamBeta]*exper  + beta[3+(iSec-1)*nParamBeta]*exper*exper + lambda[iSec-1] );
						}

						Cost[iSec-1] = 0.0;
					}
					else
					{
						experChg = gamma[iSec]*exper;

						if (typ == 1)
						{
							w[iSec-1] = rsk[iSec-1] * exp( beta[0+(iSec-1)*nParamBeta]*female + beta[1+(iSec-1)*nParamBeta]*educ
											         + beta[2+(iSec-1)*nParamBeta]*experChg  + beta[3+(iSec-1)*nParamBeta]*experChg*experChg );
						}
						else if (typ == 2)
						{
							w[iSec-1] = rsk[iSec-1] * exp( beta[0+(iSec-1)*nParamBeta]*female + beta[1+(iSec-1)*nParamBeta]*educ
											         + beta[2+(iSec-1)*nParamBeta]*experChg  + beta[3+(iSec-1)*nParamBeta]*experChg*experChg + lambda[iSec-1] );
						}

						Cost[iSec-1] = exp( xi[iLagSec-1 + 0 * nSector] + xi[iSec-1 + 1 * nSector] + kappa[0]*female +
										  kappa[1]*educ + kappa[2]*iAge + kappa[3]*iAge*iAge );
					}
				}

				if (elig == 0)
				{
					benefit = WA;
				}
				else
				{
					benefit = fmin(compRate * lastWage, UIbar);
				}


				// Compute value functions and integrate ove eta and eps
				MeanEmax = 0.0;
				for (iDraw = 0; iDraw < nDraws; iDraw++)
				{
					VM = 0.0;
					for (iSec = 1; iSec < nSector + 1; iSec++)
					{
						V[iSec-1] = (1 - Delta[iDeltaSE])*(w[iSec-1]*exp(eps[iEpsEmax]) + rho*Emax[iSec-1] - Cost[iSec-1]) + Delta[iDeltaSE]*(benefit + rho*Emax0);

						if (V[iSec-1] < 0)
						{
							V[iSec-1] = 0.00001;
						}

						// Find max
						if (V[iSec-1] > VM)
						{
							VM = V[iSec-1];
						}

					}

					SumExpV = 0.0;
					for (iSec = 1; iSec < nSector + 1; iSec++)
					{
						SumExpV += exp( (V[iSec-1] - VM)/nu[0] );
					}

					MeanEmax += VM + nu[0]*log(SumExpV);
				}
				MeanEmax = MeanEmax / nDraws;


				// Fill regression matrix and vector
				X[iN + 0  * Intp] = 1.0;
				X[iN + 1  * Intp] = rsk[0];
				X[iN + 2  * Intp] = rsk[1];
				X[iN + 3  * Intp] = rsk[2];
				X[iN + 4  * Intp] = rsk[3];
				X[iN + 5  * Intp] = rsk[4];
				X[iN + 6  * Intp] = lastWage;
				X[iN + 7  * Intp] = exper;
				X[iN + 8  * Intp] = rsk[0]*rsk[0];
				X[iN + 9  * Intp] = rsk[1]*rsk[1];
				X[iN + 10 * Intp] = rsk[2]*rsk[2];
				X[iN + 11 * Intp] = rsk[3]*rsk[3];
				X[iN + 12 * Intp] = rsk[4]*rsk[4];
				X[iN + 13 * Intp] = lastWage*lastWage;
				X[iN + 14 * Intp] = exper*exper;
				X[iN + 15 * Intp] = rsk[0]*rsk[1];
				X[iN + 16 * Intp] = rsk[0]*rsk[2];
				X[iN + 17 * Intp] = rsk[0]*rsk[3];
				X[iN + 18 * Intp] = rsk[0]*rsk[4];
				X[iN + 19 * Intp] = rsk[0]*lastWage;
				X[iN + 20 * Intp] = rsk[0]*exper;
				X[iN + 21 * Intp] = rsk[1]*rsk[2];
				X[iN + 22 * Intp] = rsk[1]*rsk[3];
				X[iN + 23 * Intp] = rsk[1]*rsk[4];
				X[iN + 24 * Intp] = rsk[1]*lastWage;
				X[iN + 25 * Intp] = rsk[1]*exper;
				X[iN + 26 * Intp] = rsk[2]*rsk[3];
				X[iN + 27 * Intp] = rsk[2]*rsk[4];
				X[iN + 28 * Intp] = rsk[2]*lastWage;
				X[iN + 29 * Intp] = rsk[2]*exper;
				X[iN + 30 * Intp] = rsk[3]*rsk[4];
				X[iN + 31 * Intp] = rsk[3]*lastWage;
				X[iN + 32 * Intp] = rsk[3]*exper;
				X[iN + 33 * Intp] = rsk[4]*lastWage;
				X[iN + 34 * Intp] = rsk[4]*exper;
				X[iN + 35 * Intp] = lastWage*exper;

				Y[iN] = MeanEmax;

			} // End iN loop


			// OLS
			fcn_linreg(Y, X, Intp, nReg, coef, &SST, &SSE);

			// Store coefficients in PI_SE
			for (iReg = 0; iReg < nReg; iReg++)
			{
				PI_SE[iPIEmax] = coef[iReg];
			}

			// Store R-squared
			R_sq_SE[iRsq] = 1.0 - SSE / SST;

		} // End iLagSec loop
	} // End iAge loop


	// Free allocated memory for random number Generator
	gsl_rng_free (s);

	// Free memory allocated to matrices and vectors
	free(eps);
	free(rsk);
	free(PI_Coef);
	free(Emax);
	free(w);
	free(V);
	free(Cost);
	free(beta);
	free(sigma);
	free(gamma);
	free(kappa);
	free(xi);
	free(nu);
	free(lambda);
	free(phi);

	free(X);
	free(Y);
	free(coef);

} // End PI_SE



//-------------------------------------------------------------------------------------
// This function computes the EMAX coefficients when agents have forward looking
// expectations.
//-------------------------------------------------------------------------------------
void Emax_RE(double *param, int iFemale, int iEduc, int iElig, int iType, double *avec, double *Low_RSk, double *Up_RSk, double Low_Wage, double Up_Wage, double *PI_SE, double *PI_RE,
	         double *R_sq_RE)
{
	//---------------------------------------------------------------------------------
	// Initialize
	//---------------------------------------------------------------------------------
	double *eps, *rsk, *beta, *sigma, *gamma, *kappa, *xi, *nu, *lambda, *phi, *PI_Coef, *Emax, *rsk_tomorrow, *w, *V, *Cost, SumExpV, MeanEmax, VM, lastWage, benefit, Emax0, exper, ExperTomorrow, experChg;
	int iYear, iAge, iDraw, iSec, iLagSec, iN, iReg;
	int female, educ, elig, typ;

	// Emax regression
	double *X, *Y, *coef, SST, SSE;

	// GSL Random Number Generator declarations
	double unif;
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng * s;

	//---------------------------------------------------------------------------------
	// Allocate memory
	//---------------------------------------------------------------------------------
	eps     = calloc((nAge-1)*nDraws*nSector,sizeof(double));
	rsk     = calloc(nSector,sizeof(double));
	beta    = calloc(nParamBeta*nSector,sizeof(double));
    sigma   = calloc(nSector,sizeof(double));
    gamma   = calloc(nSector+1,sizeof(double));
    kappa   = calloc(4,sizeof(double));
    xi      = calloc(2*nSector,sizeof(double));
    nu      = calloc(1,sizeof(double));
    lambda  = calloc(nSector,sizeof(double));
	phi     = calloc(4,sizeof(double));
	PI_Coef = calloc(nReg,sizeof(double));
	Emax    = calloc(nSector,sizeof(double));
	w       = calloc(nSector,sizeof(double));
	V       = calloc(nSector,sizeof(double));
	Cost    = calloc(nSector,sizeof(double));
	rsk_tomorrow = calloc(nSector,sizeof(double));

	// Emax regressions
	X    = calloc(Intp*nReg,sizeof(double));
	Y    = calloc(Intp,sizeof(double));
	coef = calloc(Intp,sizeof(double));

	// Initialize Emax0
	Emax0 = 0.0;

	// Definitions
	#define iAvec 		(iSec-1) + (iYear - 1) * nSector
	#define iPIRE1		iReg + iYear     * nReg + iAge     * (LastYr_sim-(Yr_shock+1)-1) * nReg + iLagSec * (nAge-1) * (LastYr_sim-(Yr_shock+1)-1) * nReg + iFemale * (nSector+1) * (nAge-1) * (LastYr_sim-(Yr_shock+1)-1) * nReg + iEduc * CatGender * (nSector+1) * (nAge-1) * (LastYr_sim-(Yr_shock+1)-1) * nReg + iElig * CatEduc * CatGender * (nSector+1) * (nAge-1) * (LastYr_sim-(Yr_shock+1)-1) * nReg + iType * CatElig * CatEduc * CatGender * (nSector+1) * (nAge-1) * (LastYr_sim-(Yr_shock+1)-1) * nReg
	#define iPIREEmax	iReg + (iYear-1) * nReg + (iAge-1) * (LastYr_sim-(Yr_shock+1)-1) * nReg + iLagSec * (nAge-1) * (LastYr_sim-(Yr_shock+1)-1) * nReg + iFemale * (nSector+1) * (nAge-1) * (LastYr_sim-(Yr_shock+1)-1) * nReg + iEduc * CatGender * (nSector+1) * (nAge-1) * (LastYr_sim-(Yr_shock+1)-1) * nReg + iElig * CatEduc * CatGender * (nSector+1) * (nAge-1) * (LastYr_sim-(Yr_shock+1)-1) * nReg + iType * CatElig * CatEduc * CatGender * (nSector+1) * (nAge-1) * (LastYr_sim-(Yr_shock+1)-1) * nReg
	#define iRsqRE		(iYear-1) + (iAge-1) * (LastYr_sim-(Yr_shock+1)-1) + iLagSec * (nAge-1) * (LastYr_sim-(Yr_shock+1)-1) + iFemale * (nSector+1) * (nAge-1) * (LastYr_sim-(Yr_shock+1)-1) + iEduc * CatGender * (nSector+1) * (nAge-1) * (LastYr_sim-(Yr_shock+1)-1) + iElig * CatEduc * CatGender * (nSector+1) * (nAge-1) * (LastYr_sim-(Yr_shock+1)-1) + iType * CatElig * CatEduc * CatGender * (nSector+1) * (nAge-1) * (LastYr_sim-(Yr_shock+1)-1)

	//--------------------------------------------------------------------
	// Unpack parameters
	//--------------------------------------------------------------------
	unpackPars(param, beta, sigma, gamma, xi, kappa, nu, lambda, phi);


	//---------------------------------------------------------------------------------
	// Algorithm goes from 2307 to 2160 as workers have static expectations in
	// last year. I subtract 2 as indexation starts at 0.
	//---------------------------------------------------------------------------------
	for (iYear = LastYr_sim - (Yr_shock+1) - 2; iYear > 0; iYear--)
	{
		//-----------------------------------------------------------------------------
		// Random numbers for Monte Carlo integration over wage shocks. Draw only
		// nDraws for each age to save computational time.
		//-----------------------------------------------------------------------------
		T = gsl_rng_mt19937; // Algorithm
		r = gsl_rng_alloc(T); // Allocate memory
		gsl_rng_set(r,12345678); // Set seed

		for (iAge = 1; iAge < nAge; iAge++)
		{
			for (iDraw = 0; iDraw < nDraws; iDraw++)
			{
				for (iSec = 1; iSec < nSector + 1; iSec++)
				{
					eps[iEpsEmax] = gsl_ran_gaussian(r, sigma[iSec-1]);
				}
			}
		}


		// Free allocated memory for random number generator
		gsl_rng_free(r);

		// Random number generator for experience and tenure draws
		T = gsl_rng_mt19937; // Algorithm
		s = gsl_rng_alloc(T); // Allocate memory
		gsl_rng_set(s,987656); // Set seed


		for (iAge = nAge - 1; iAge > 0; iAge--)
		{
			for (iLagSec = 0; iLagSec < nSector + 1; iLagSec++)
			{
				for (iN = 0; iN < Intp; iN++)
				{
					//-------------------------------------------------------------
					// Random number generation - skill returns, exper and tenure
					//-------------------------------------------------------------
					exper  = gsl_rng_uniform_int(s, iAge + 5) + 2.0; // Draw exper from 2 to iAge + 5

					for (iSec = 1; iSec < nSector + 1; iSec++)
					{
						unif = gsl_rng_uniform_pos(s);
						rsk[iSec-1] = Low_RSk[iSec-1] + (Up_RSk[iSec-1] - Low_RSk[iSec-1])*unif;
					}

					unif = gsl_rng_uniform_pos(s);
					lastWage = Low_Wage + (Up_Wage - Low_Wage) * unif;

					//-------------------------------------------------------------
					// Compute Emax for iAge < LastAge
					//-------------------------------------------------------------
					for (iSec = 1; iSec < nSector + 1; iSec++)
					{
						if (iAge < nAge - 1)
						{
							if (iLagSec == iSec)
							{
								ExperTomorrow = exper + 1;
							}
							else
							{
								ExperTomorrow = gamma[iSec]*exper;
							}

							rsk_tomorrow[iSec-1] = avec[iAvec] * rsk[iSec-1];

							// Workers have rational expectations up to 2007
							// and static expectations from 2008
							if (iYear <= nYear - 3)
							{
								for (iReg = 0; iReg < nReg; iReg++)
								{
									PI_Coef[iReg] = PI_RE[iPIRE1];
								}
							}
							else
							{
								for (iReg = 0; iReg < nReg; iReg++)
								{
									PI_Coef[iReg] = PI_SE[iPI1];
								}
							}

							Emax[iSec-1] = Emax_Hat(PI_Coef,rsk_tomorrow,ExperTomorrow,lastWage);
							Emax0 = Emax_Hat(PI_Coef,rsk_tomorrow,gamma[0]*exper,lastWage);
						}
						else
						{
							Emax[iSec-1] = 0.0;
							Emax0 = 0.0;
						}
					}


					//-------------------------------------------------------------
					// Wage and and value function
					//-------------------------------------------------------------
					female = iFemale;
					educ = iEduc;
					elig = iElig;
					typ = iType + 1;

					for (iSec = 1; iSec < nSector + 1; iSec++)
					{
						if (iLagSec == iSec)
						{
							if (typ == 1)
							{
								w[iSec-1] = rsk[iSec-1] * exp( beta[0+(iSec-1)*nParamBeta]*female + beta[1+(iSec-1)*nParamBeta]*educ
												         + beta[2+(iSec-1)*nParamBeta]*exper  + beta[3+(iSec-1)*nParamBeta]*exper*exper );
							}
							else if (typ == 2)
							{
								w[iSec-1] = rsk[iSec-1] * exp( beta[0+(iSec-1)*nParamBeta]*female + beta[1+(iSec-1)*nParamBeta]*educ
												         + beta[2+(iSec-1)*nParamBeta]*exper  + beta[3+(iSec-1)*nParamBeta]*exper*exper + lambda[iSec-1] );
							}

							Cost[iSec-1] = 0.0;
						}
						else if (iLagSec == 0)
						{
							if (typ == 1)
							{
								w[iSec-1] = rsk[iSec-1] * exp( beta[0+(iSec-1)*nParamBeta]*female + beta[1+(iSec-1)*nParamBeta]*educ
												         + beta[2+(iSec-1)*nParamBeta]*exper  + beta[3+(iSec-1)*nParamBeta]*exper*exper );
							}
							else if (typ == 2)
							{
								w[iSec-1] = rsk[iSec-1] * exp( beta[0+(iSec-1)*nParamBeta]*female + beta[1+(iSec-1)*nParamBeta]*educ
												         + beta[2+(iSec-1)*nParamBeta]*exper  + beta[3+(iSec-1)*nParamBeta]*exper*exper + lambda[iSec-1] );
							}

							Cost[iSec-1] = 0.0;
						}
						else
						{
							experChg = gamma[iSec]*exper;

							if (typ == 1)
							{
								w[iSec-1] = rsk[iSec-1] * exp( beta[0+(iSec-1)*nParamBeta]*female + beta[1+(iSec-1)*nParamBeta]*educ
												         + beta[2+(iSec-1)*nParamBeta]*experChg  + beta[3+(iSec-1)*nParamBeta]*experChg*experChg );
							}
							else if (typ == 2)
							{
								w[iSec-1] = rsk[iSec-1] * exp( beta[0+(iSec-1)*nParamBeta]*female + beta[1+(iSec-1)*nParamBeta]*educ
												         + beta[2+(iSec-1)*nParamBeta]*experChg  + beta[3+(iSec-1)*nParamBeta]*experChg*experChg + lambda[iSec-1] );
							}

							Cost[iSec-1] = exp( xi[iLagSec-1 + 0 * nSector] + xi[iSec-1 + 1 * nSector] + kappa[0]*female +
											  kappa[1]*educ + kappa[2]*iAge + kappa[3]*iAge*iAge );
						}
					}

					if (elig == 0)
					{
						benefit = WA;
					}
					else
					{
						benefit = fmin(compRate * lastWage, UIbar);
					}



					// Compute value functions and integrate over eta and eps
					MeanEmax = 0.0;
					for (iDraw = 0; iDraw < nDraws; iDraw++)
					{
						VM = 0.0;
						for (iSec = 1; iSec < nSector + 1; iSec++)
						{
							V[iSec-1] = (1 - Delta[iDeltaRE])*(w[iSec-1]*exp(eps[iEpsEmax])  + rho*Emax[iSec-1] - Cost[iSec-1]) + Delta[iDeltaRE]*(benefit + rho*Emax0);

							if (V[iSec-1] < 0)
							{
								V[iSec-1] = 0.00001;
							}

							// Find max
							if (V[iSec-1] > VM)
							{
								VM = V[iSec-1];
							}

						}

						SumExpV = 0.0;
						for (iSec = 1; iSec < nSector + 1; iSec++)
						{
							SumExpV += exp( (V[iSec-1] - VM)/nu[0] );
						}

						MeanEmax += VM + nu[0]*log(SumExpV);
					}
					MeanEmax = MeanEmax / nDraws;


					// Fill regression matrix and vector
					X[iN + 0  * Intp] = 1.0;
					X[iN + 1  * Intp] = rsk[0];
					X[iN + 2  * Intp] = rsk[1];
					X[iN + 3  * Intp] = rsk[2];
					X[iN + 4  * Intp] = rsk[3];
					X[iN + 5  * Intp] = rsk[4];
					X[iN + 6  * Intp] = lastWage;
					X[iN + 7  * Intp] = exper;
					X[iN + 8  * Intp] = rsk[0]*rsk[0];
					X[iN + 9  * Intp] = rsk[1]*rsk[1];
					X[iN + 10 * Intp] = rsk[2]*rsk[2];
					X[iN + 11 * Intp] = rsk[3]*rsk[3];
					X[iN + 12 * Intp] = rsk[4]*rsk[4];
					X[iN + 13 * Intp] = lastWage*lastWage;
					X[iN + 14 * Intp] = exper*exper;
					X[iN + 15 * Intp] = rsk[0]*rsk[1];
					X[iN + 16 * Intp] = rsk[0]*rsk[2];
					X[iN + 17 * Intp] = rsk[0]*rsk[3];
					X[iN + 18 * Intp] = rsk[0]*rsk[4];
					X[iN + 19 * Intp] = rsk[0]*lastWage;
					X[iN + 20 * Intp] = rsk[0]*exper;
					X[iN + 21 * Intp] = rsk[1]*rsk[2];
					X[iN + 22 * Intp] = rsk[1]*rsk[3];
					X[iN + 23 * Intp] = rsk[1]*rsk[4];
					X[iN + 24 * Intp] = rsk[1]*lastWage;
					X[iN + 25 * Intp] = rsk[1]*exper;
					X[iN + 26 * Intp] = rsk[2]*rsk[3];
					X[iN + 27 * Intp] = rsk[2]*rsk[4];
					X[iN + 28 * Intp] = rsk[2]*lastWage;
					X[iN + 29 * Intp] = rsk[2]*exper;
					X[iN + 30 * Intp] = rsk[3]*rsk[4];
					X[iN + 31 * Intp] = rsk[3]*lastWage;
					X[iN + 32 * Intp] = rsk[3]*exper;
					X[iN + 33 * Intp] = rsk[4]*lastWage;
					X[iN + 34 * Intp] = rsk[4]*exper;
					X[iN + 35 * Intp] = lastWage*exper;

					Y[iN] = MeanEmax;

				} // End iN loop

				// OLS
				fcn_linreg(Y, X, Intp, nReg, coef, &SST, &SSE);

				// Store coefficients in PI_SE
				for (iReg = 0; iReg < nReg; iReg++)
				{
					PI_SE[iPIEmax] = coef[iReg];
				}

				// Store R-squared
				R_sq_RE[iRsqRE] = 1.0 - SSE / SST;

			} // end iLagSec
		} // end age


		// Free GSL memory
		gsl_rng_free(s);

	} // End year loop


	// Free allocated memory
	free(eps);
	free(rsk);
	free(beta);
	free(sigma);
	free(gamma);
	free(kappa);
	free(xi);
	free(nu);
	free(lambda);
	free(phi);
	free(Emax);
	free(PI_Coef);
	free(rsk_tomorrow);
	free(w);
	free(V);
	free(Cost);

	free(X);
	free(Y);
	free(coef);

} // End Emax_RE

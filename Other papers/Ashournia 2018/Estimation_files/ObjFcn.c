//------------------------------------------------------------------------------
// ObjFcn.c
//
// By Damoun Ashournia, University of Oxford
//
//------------------------------------------------------------------------------
// Functions to compute the SMD loss function using simulated data.
//------------------------------------------------------------------------------


//------------------------------------------------------------------------
// Module Loss Function
//------------------------------------------------------------------------
double LossFcn(double *param)
{
	//--------------------------------------------------------------------
	// Declarations and definitions
	//--------------------------------------------------------------------
	double *PI_RE, *R_sq_RE, *PI_SE, *R_sq_SE;
	double *Low_RSk, *Up_RSk, *avec, *rsk_sup, *rsk_dem, *rsk_eq, *sk_sup, *sk_sup_eq, *step, *beta, *sigma, *gamma, *delta, *kappa, *xi, *nu, *lambda, *phi;
	double *rsk_tomorrow, *PI_Coef, *w, *V, *Cost, *check, *Data, *z;
	double crit, eps_aux, eta_aux, Emax, VMAX, Low_Wage, Up_Wage, benefit, lastwage, Emax0;

	// Parallelization
	int *Female, *Educ, *Elig, *Type, iParallel, tid, nthreads;

	// Auxiliary regressions
	double *auxSim;

	int *emp, *emp_eq, *iter, *flag, *yearDummies;
	int iAge, iPar, iPar2, iCh, iSec, iYear, iCoh, iN, iCnt, iReg, iCount, iSize, iLagCh, iAv, RE_loop, it, iRsk, converged, *choice;
	int age, female, educ, elig, lagCh, typ, *type;
	double exper, ExperTomorrow, term2, prob;

	double *XWage, *XOther, *Y, *coef, SST, SSE, *choiceWgt;

	// Loss function
	double Loss;

	// Probabilities for drawing unemployment shocks.
	double Punemp;

	// File pointer
	FILE *fp;


	// Running time
	clock_t cstart = 0;
  	clock_t cend = 0;

	// Allocate memory
	PI_RE        = calloc(nReg*(nYear-1)*(nAge-1)*nChoices*CatGender*CatEduc*CatElig*CatTypes,sizeof(double));
	R_sq_RE      = calloc((nYear-1)*(nAge-1)*nChoices*CatGender*CatEduc*CatElig*CatTypes,sizeof(double));
	PI_SE        = calloc(nReg*(nAge-1)*nChoices*CatGender*CatEduc*CatElig*CatTypes,sizeof(double));
	R_sq_SE      = calloc((nAge-1)*nChoices*CatGender*CatEduc*CatElig*CatTypes,sizeof(double));
	Low_RSk      = calloc(nSector,sizeof(double));
	Up_RSk       = calloc(nSector,sizeof(double));
	avec         = calloc(nSector*nYear,sizeof(double));
	rsk_sup      = calloc(nSector*MaxIt,sizeof(double));
	rsk_dem      = calloc(nSector*MaxIt,sizeof(double));
	rsk_eq       = calloc(nSector*nYear,sizeof(double));
	sk_sup       = calloc(nSector,sizeof(double));
	sk_sup_eq    = calloc(nSector*nYear,sizeof(double));
	step         = calloc(nSector,sizeof(double));
	beta         = calloc(nSector*nParamBeta,sizeof(double));
	kappa        = calloc(4,sizeof(double));
	xi           = calloc(2*nChoices,sizeof(double));
	sigma        = calloc(nChoices,sizeof(double));
	gamma        = calloc(nChoices,sizeof(double));
    delta        = calloc(nChoices,sizeof(double));
	nu           = calloc(1,sizeof(double));
	lambda       = calloc(nChoices,sizeof(double));
    phi          = calloc(4,sizeof(double));
	rsk_tomorrow = calloc(nSector,sizeof(double));
	PI_Coef      = calloc(nReg,sizeof(double));
	w            = calloc(nChoices,sizeof(double));
	V            = calloc(nChoices,sizeof(double));
	Cost         = calloc(nChoices,sizeof(double));
	check        = calloc(nSector*nYear,sizeof(double));
	z            = calloc(nSector*nYear,sizeof(double));

	emp         = calloc(nChoices*MaxIt,sizeof(int));
	emp_eq      = calloc(nChoices*nYear,sizeof(int));
	iter        = calloc(nYear,sizeof(int));
	flag        = calloc(nSector,sizeof(int));
	yearDummies = calloc(nPop*nAge*nYear*nYear,sizeof(int));
	choice      = calloc(nChoices,sizeof(int));
	choiceWgt   = calloc(nChoices,sizeof(double));

	Female = calloc(CatGender*CatEduc*CatElig*CatTypes,sizeof(double));
	Educ   = calloc(CatGender*CatEduc*CatElig*CatTypes,sizeof(double));
	Elig   = calloc(CatGender*CatEduc*CatElig*CatTypes,sizeof(double));
	Type   = calloc(CatGender*CatEduc*CatElig*CatTypes,sizeof(double));
	type   = calloc(nPop*nCoh,sizeof(int));

	// Auxiliary regressions
	Data   = calloc(nPop*nAge*nYear*13,sizeof(double));
	auxSim = calloc(nParamAux,sizeof(double));

	XWage  = calloc(nPop*nAge*nYear*nRegWage,sizeof(double));
	XOther = calloc(nPop*nAge*nYear*nRegOther,sizeof(double));
	Y      = calloc(nPop*nAge*nYear,sizeof(double));
	coef   = calloc(nPop*nAge*nYear,sizeof(double));


	// GSL Random Number Generator declarations
	int unemp;
	const gsl_rng_type * T;
	gsl_rng * r;

	T = gsl_rng_mt19937; // Algorithm
	r = gsl_rng_alloc(T); // Allocate memory
	gsl_rng_set(r,5566778); // Set seed


	// Definitions
	#define iPopCohYearNext		iN + iCoh * nPop + (iYear+1) * nCoh * nPop
	#define size				nPop * nAge * nYear
	#define iAvecObj			iSec + iYear * nSector
	#define iPISE				iReg + (age-30)*nReg + lagCh*(nAge-1)*nReg + female*nChoices*(nAge-1)*nReg + educ*CatGender*nChoices*(nAge-1)*nReg + elig*CatEduc*CatGender*nChoices*(nAge-1)*nReg + (typ-1)*CatElig*CatEduc*CatGender*nChoices*(nAge-1)*nReg
	#define iPIRE				iReg + iYear*nReg + (age-30)*(nYear-1)*nReg + lagCh*(nAge-1)*(nYear-1)*nReg + female*nChoices*(nAge-1)*(nYear-1)*nReg + educ*CatGender*nChoices*(nAge-1)*(nYear-1)*nReg + elig*CatEduc*CatGender*nChoices*(nAge-1)*(nYear-1)*nReg + (typ-1)*CatElig*CatEduc*CatGender*nChoices*(nAge-1)*(nYear-1)*nReg
	#define iCohort				iCoh + iYear * nCoh


	//--------------------------------------------------------------------
	// Unpack parameters
	//--------------------------------------------------------------------
	unpackPars(param, beta, sigma, gamma, xi, kappa, nu, lambda, phi);

	//---------------------------------------------------------------------------------
	// Upper and lower bounds for draws of skill returns and wages
	//---------------------------------------------------------------------------------
	Low_RSk[0] = 50;
	Up_RSk[0]  = 350;

	Low_RSk[1] = 55;
	Up_RSk[1]  = 400;

	Low_RSk[2] = 60;
	Up_RSk[2]  = 400;

	Low_RSk[3] = 65;
	Up_RSk[3]  = 550;

	Low_RSk[4] = 60;
	Up_RSk[4]  = 600;

	Low_Wage = 65;
	Up_Wage  = 1200;


	// Initiallize avec to ones
	for (iAv = 0; iAv < nSector*nYear; iAv++)
	{
		avec[iAv] = 1.0;
	}

	// Types
	for (iCoh = 0; iCoh < nCoh; iCoh++)
	{
		for (iN = 0; iN < nPop; iN++)
		{
			female = FemaleData[iPopCoh];
			educ   = EducData[iPopCoh];
			exper  = ExperData[iN + iCoh * nPop + 0 * nCoh * nPop];

			term2 = exp(phi[0] + phi[1]*female + phi[2]*educ + phi[3]*exper);
			prob  = 1.0 / (1.0 + term2);
			if (unif_type[iTyp] < prob)
			{
				type[iTyp] = 1;
			}
			else
			{
				type[iTyp] = 2;
			}
		}
	}


	//--------------------------------------------------------------------
	// Female, Educ, Elig and Type for parallelization
	//--------------------------------------------------------------------
	Female[0]  = 0;
	Female[1]  = 0;
	Female[2]  = 0;
	Female[3]  = 0;
	Female[4]  = 0;
	Female[5]  = 0;
	Female[6]  = 0;
	Female[7]  = 0;
	Female[8]  = 1;
    Female[9]  = 1;
    Female[10] = 1;
    Female[11] = 1;
    Female[12] = 1;
    Female[13] = 1;
    Female[14] = 1;
    Female[15] = 1;

	Educ[0]  = 0;
	Educ[1]  = 0;
	Educ[2]  = 0;
	Educ[3]  = 0;
	Educ[4]  = 1;
	Educ[5]  = 1;
	Educ[6]  = 1;
	Educ[7]  = 1;
	Educ[8]  = 0;
	Educ[9]  = 0;
	Educ[10] = 0;
	Educ[11] = 0;
	Educ[12] = 1;
	Educ[13] = 1;
	Educ[14] = 1;
	Educ[15] = 1;

	Elig[0]  = 0;
	Elig[1]  = 0;
	Elig[2]  = 1;
	Elig[3]  = 1;
	Elig[4]  = 0;
	Elig[5]  = 0;
	Elig[6]  = 1;
	Elig[7]  = 1;
	Elig[8]  = 0;
	Elig[9]  = 0;
	Elig[10] = 1;
	Elig[11] = 1;
	Elig[12] = 0;
	Elig[13] = 0;
	Elig[14] = 1;
	Elig[15] = 1;

	Type[0]  = 0;
	Type[1]  = 1;
	Type[2]  = 0;
	Type[3]  = 1;
	Type[4]  = 0;
	Type[5]  = 1;
	Type[6]  = 0;
	Type[7]  = 1;
	Type[8]  = 0;
	Type[9]  = 1;
	Type[10] = 0;
	Type[11] = 1;
	Type[12] = 0;
	Type[13] = 1;
	Type[14] = 0;
	Type[15] = 1;

	//--------------------------------------------------------------------
	// Solve model with static expectations
	//--------------------------------------------------------------------
	cstart = clock();

	// Fork a team of threads
	#pragma omp parallel private(tid, iParallel) shared(PI_SE, R_sq_SE)
	{

		tid = omp_get_thread_num();
    	// Only master thread does this
    	if (tid == 0)
      	{
			nthreads = omp_get_num_threads();
			// printf("Number of threads = %d\n", nthreads);
      	}

      	#pragma omp for nowait
		for (iParallel = 0; iParallel < CatGender * CatEduc * CatElig * CatTypes; iParallel++)
		{
			Emax_SE(tid, param, Female[iParallel], Educ[iParallel], Elig[iParallel], Type[iParallel], Low_RSk, Up_RSk, Low_Wage, Up_Wage, PI_SE, R_sq_SE);
		}
	}

	cend = clock();
	// printf ("EMAX SE: %.3f cpu sec\n", ((double)cend - (double)cstart)* 1.0e-6);

	// Write Emax regression coefficients for static expectations to file
	fp = fopen("Output/Emax_Coef.txt", "w");
	fprintf(fp, "Age\tFemale\tEduc\tElig\tType\tLagCh\tReg\tPI_SE\n");
	for (iAge = 0; iAge < nAge-1; iAge++)
	{
		age = iAge + 30;
		for (female = 0; female < CatGender; female++)
		{
			for (educ = 0; educ < CatEduc; educ++)
			{
				for (elig = 0; elig < CatElig; elig++)
				{
					for (typ = 1; typ < CatTypes + 1; typ++)
					{
						for (lagCh = 0; lagCh < nChoices; lagCh++)
						{
							for (iReg = 0; iReg < nReg; iReg++)
							{
								fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\n", age, female, educ, elig, typ, lagCh, iReg, PI_SE[iPISE]);
							}
						}
					}
				}
			}
		}
	}
	fclose(fp);


	//--------------------------------------------------------------------
	// Solve model with forward looking expectations
	//--------------------------------------------------------------------
	RE_loop = 1;

	while (RE_loop == 1)//(RE_loop <= 2) for forward-looking expectations only
	{
		if (RE_loop >= 2)
		{
			cstart = clock();

			// Fork a team of threads
			#pragma omp parallel private(tid, iParallel) shared(PI_RE, R_sq_RE)
			{

				tid = omp_get_thread_num();
				// Only master thread does this
				if (tid == 0)
				{
					nthreads = omp_get_num_threads();
					//printf("Number of threads = %d\n", nthreads);
				}


				#pragma omp for nowait
				for (iParallel = 0; iParallel < CatGender * CatEduc * CatElig * CatTypes; iParallel++)
				{
					Emax_RE(param, Female[iParallel], Educ[iParallel], Elig[iParallel], Type[iParallel], avec, Low_RSk, Up_RSk, Low_Wage, Up_Wage, PI_SE, PI_RE, R_sq_RE);
				}
			}

			cend = clock();
			// printf ("EMAX RE: %.3f cpu sec\n", ((double)cend - (double)cstart)* 1.0e-6);
			// printf("RE_loop = %d; converged = %d\n", RE_loop, converged);
		}

		// Set first iteration skill return
		for (iSec = 0; iSec < nSector; iSec++)
		{
			rsk_sup[iSec + 0 * nSector] = rsk_global[iSec];
		}

		crit = 0.001; // Convergence criteria for rsk_sup and rsk_dem

		for (iYear = 0; iYear < nYear; iYear++)
		{
			it = 0; // iteration counter

			if (iYear > 0)
			{
				// Use last period result as initial point for rsk_sup
				for (iRsk = 0; iRsk < nSector*MaxIt; iRsk++)
				{
					rsk_sup[iRsk] = 0.0;
					rsk_dem[iRsk] = 0.0;
				}

				for (iSec = 0; iSec < nSector; iSec++)
				{
					rsk_sup[iSec + 0 * nSector] = fmin(fmax(rsk_eq[iSec + (iYear-1)*nSector],Low_RSk[iSec]),Up_RSk[iSec]);
				}
			}

			converged = 0;

			while (converged==0 && it < MaxIt)
			{
				for (iCh = 0; iCh < nChoices; iCh++)
				{
					if (iCh>0)
					{
						sk_sup[iCh-1] = 0.0;
					}
					emp[iCh + it * nChoices] = 0;
				}

				//--------------------------------------------------------
				// Compute aggregate skill supply
				//--------------------------------------------------------
				for (iCoh = 0; iCoh < nCoh; iCoh++)
				{
					for (iN = 0; iN < nPop; iN++)
					{
						for (iCh = 0; iCh < nChoices; iCh++)
						{
							age      = (iYear + 1996) - (iCoh + 1931);
							female   = FemaleData[iPopCoh];
							educ     = EducData[iPopCoh];
							lagCh    = LagChoiceData[iPopCohYear];
							exper    = ExperData[iPopCohYear];
							lastwage = LagWageData[iPopCohYear];
							eps_aux  = sigma[iCh] * epsSim[iEps];
							eta_aux  = -log(-log(etaSim[iEta]))*nu[0] - 0.577215665*nu[0];
							typ      = type[iTyp];

							// Eligibility for UI
							if (UIFundData[iPopCoh] == 1 && UnempYearsData[iPopCohYear] <= 4)
							{
								elig = 1;
								EligData[iPopCohYear] = 1;
							}
							else
							{
								elig = 0;
								EligData[iPopCohYear] = 0;
							}

							if (age <= LastAge && age >= FirstAge)
							{
								if (age < LastAge)
								{
									for (iSec = 0; iSec < nSector; iSec++)
									{
										rsk_tomorrow[iSec] = avec[iAvecObj] * rsk_sup[iSec + it * nSector];
									}

									if (lagCh == iCh)
									{
										ExperTomorrow = exper + 1.0;
									}
									else
									{
										ExperTomorrow = gamma[iCh]*exper;
									}

									if (iYear < nYear - 3 && RE_loop >= 2)
									{
										for (iReg = 0; iReg < nReg; iReg++)
										{
											// PI_Coef[iReg] = PI_RE[iPIRE];
											PI_Coef[iReg] = PI_SE[iPISE];
										}
									}
									else
									{
										for (iReg = 0; iReg < nReg; iReg++)
										{
											PI_Coef[iReg] = PI_SE[iPISE];
										}
									}

									Emax = Emax_Hat(PI_Coef,rsk_tomorrow,ExperTomorrow,lastwage);
									Emax0 = Emax_Hat(PI_Coef,rsk_tomorrow,gamma[0]*exper,lastwage);
								}
								else
								{
									Emax = 0.0;
									Emax0 = 0.0;
								}

								if (lagCh == iCh)
								{
									// exper = exper;
								}
								else
								{
									exper = gamma[iCh]*exper;
								}

								// Wage and value function
								if (typ == 1)
								{
									w[iCh] = rsk_sup[iCh + it * nSector] *
										exp(  beta[0+iCh*nParamBeta]*female + beta[1+iCh*nParamBeta]*educ
											+ beta[2+iCh*nParamBeta]*exper  + beta[3+iCh*nParamBeta]*exper*exper) * exp(eps_aux);
								}
								else if (typ == 2)
								{
									w[iCh] = rsk_sup[iCh + it * nSector] *
										exp(  beta[0+iCh*nParamBeta]*female + beta[1+iCh*nParamBeta]*educ
											+ beta[2+iCh*nParamBeta]*exper  + beta[3+iCh*nParamBeta]*exper*exper + lambda[iCh]) * exp(eps_aux);
								}

								if (lagCh == iCh || lagCh == 0)
								{
									Cost[iCh] = 0.0;
								}
								else
								{
									Cost[iCh] = exp( xi[lagCh-1 + 0 * nSector] + xi[iCh + 1 * nSector] + kappa[0]*female +
										  kappa[1]*educ + kappa[2]*(age-30) + kappa[3]*(age-30)*(age-30) );
								}


								if (elig == 0)
								{
									benefit = WA;
								}
								else
								{
									benefit = fmin(compRate * lastwage, UIbar);
								}

								V[iCh] = (1 - delta[iDelta])*(w[iCh] - Cost[iCh] + rho*Emax) + delta[iDelta]*(benefit + rho*Emax0) + eta_aux;
							}
						} // end iCh


						if (age <= LastAge && age >= FirstAge)
						{
							// Find choice (= maxV) and store next period states
							VMAX = Find_Max(V, nSector);

							// Unemployment shock
  							Punemp = gsl_rng_uniform_pos(r);

  							if (Punemp <= delta[iDelta])
  							{
  								unemp = 1;
  							}
  							else
  							{
  								unemp = 0;
  							}

  							if (unemp == 0)
  							{
  								for (iSec = 1; iSec < nSector + 1; iSec++)
								{
									if (VMAX == V[iSec-1])
									{
										ChoiceData[iPopCohYear]      = iSec;
										WageData[iPopCohYear]        = w[iSec-1];
										CostSim[iPopCohYear]         = Cost[iSec-1];
										VSim[iPopCohYear]            = VMAX;
										emp[iSec + it * nChoices] = emp[iSec + it * nChoices] + 1;
										sk_sup[iSec-1]               = sk_sup[iSec-1] + CohortWgt[iCohort] * (w[iSec-1] / rsk_sup[(iSec-1) + it * nSector]);
									}
								}
  							}
  							else
  							{
  								ChoiceData[iPopCohYear]   = 0;
  								WageData[iPopCohYear]     = benefit;
  								CostSim[iPopCohYear]      = 0;
  								VSim[iPopCohYear]         = VMAX;
  								emp[0 + it * nChoices] = emp[0 + it * nChoices] + 1;
  							}



							if (age < LastAge && iYear < nYear - 1)
							{
								if (ChoiceData[iPopCohYear] == lagCh && lagCh > 0)
								{
									ExperData[iPopCohYearNext]      = exper + 1.0;
									LagWageData[iPopCohYearNext]    = WageData[iPopCohYear];
									UnempYearsData[iPopCohYearNext] = 0;
								}
								else if (ChoiceData[iPopCohYear] == 0)
								{
									ExperData[iPopCohYearNext]      = gamma[0]*exper;
									LagWageData[iPopCohYearNext]    = LagWageData[iPopCohYear];
									UnempYearsData[iPopCohYearNext] = UnempYearsData[iPopCohYear] + 1;
								}
								else
								{
									ExperData[iPopCohYearNext]      = gamma[ChoiceData[iPopCohYear]]*exper;
									LagWageData[iPopCohYearNext]    = WageData[iPopCohYear];
									UnempYearsData[iPopCohYearNext] = 0;
								}

								LagChoiceData[iPopCohYearNext] = ChoiceData[iPopCohYear];

							}
						}

					} // end iN
				} // end iCoh



				// rsk is zero if there is no labor supply for some occupation, so it will be infinity.
				// If labor supply is zero, double rsk instead
				for (iSec = 1; iSec < nSector + 1; iSec++)
				{
					if (sk_sup[iSec-1] > 0)
					{
						rsk_dem[(iSec-1) + it * nSector] =
							fmin(IncomeShares[iYear + (iSec-1)*nYear + 0 * nSector * nYear] * OutputData[iYear + (iSec-1)*nYear] / sk_sup[iSec-1], 1000.0);
					}
					else
					{
						rsk_dem[(iSec-1) + it * nSector] =
							fmin(2*rsk_sup[(iSec-1) + it * nSector], 1000.0);
					}

					check[(iSec-1) + iYear * nSector] =
						fabs(log(rsk_sup[(iSec-1) + it * nSector]) - log(rsk_dem[(iSec-1) + it * nSector]));
				}

				// Check for convergence
				if (check[0 + iYear*nSector] <= crit && check[1 + iYear*nSector] <= crit && check[2 + iYear*nSector] <= crit && check[3 + iYear*nSector] <= crit && check[4 + iYear*nSector] <= crit)
				{
					converged = 1;
				}

				if (converged == 0 && it < MaxIt - 1)
				{
					// Adjust the step for updating the skill price
					if (it < MaxIt/2)
					{
						for (iSec = 1; iSec < nSector + 1; iSec++)
						{
							step[iSec-1] = 1.0/6.0;
						}
					}
					else
					{
						for (iSec = 1; iSec < nSector + 1; iSec++)
						{
							if (SameSign(1.0, rsk_dem[(iSec-1) + it*nSector]-rsk_dem[(iSec-1) + (it-1)*nSector]) !=
								SameSign(1.0, rsk_dem[(iSec-1) + (it-1)*nSector]-rsk_dem[(iSec-1) + (it-2)*nSector])
								&& flag[iSec-1] == 0)
							{
								flag[iSec-1] = it;
								step[iSec-1] = step[iSec-1] / 1.3;
							}

							if (SameSign(1.0, rsk_dem[(iSec-1) + it*nSector]-rsk_dem[(iSec-1) + (it-1)*nSector]) !=
								SameSign(1.0, rsk_dem[(iSec-1) + (it-1)*nSector]-rsk_dem[(iSec-1) + (it-2)*nSector])
								&& it >= flag[iSec-1] + 5)
							{
								flag[iSec-1] = it;
								step[iSec-1] = step[iSec-1] / 1.3;
							}
								step[iSec-1] = step[iSec-1] / 1.3;
						}
					}

					// Update the return to skill
					for (iSec = 1; iSec < nSector + 1; iSec++)
					{
						rsk_sup[(iSec-1) + (it+1)*nSector] = rsk_sup[(iSec-1) + it*nSector] -
							step[iSec-1] * (rsk_sup[(iSec-1) + it*nSector] - rsk_dem[(iSec-1) + it*nSector]);
					}
				}
				else
				{
					// Record equilibrium quantities
					for (iSec = 0; iSec < nSector + 1; iSec++)
					{
						if (iSec > 0)
						{
							rsk_eq[(iSec-1) + iYear*nSector]    = rsk_dem[(iSec-1) + it*nSector];
							sk_sup_eq[(iSec-1) + iYear*nSector] = sk_sup[(iSec-1)];
							emp_eq[iSec + iYear*nChoices]    = emp[iSec + it*nChoices];
							iter[iYear] = it;
						}
						else
						{
							emp_eq[iSec + iYear*nChoices]    = emp[iSec + it*nChoices];
							iter[iYear] = it;
						}

					}
				}

				it += 1;

			} // end convergence of skill return for iYear

			//printf("%d %f %f\n", iYear, rsk_eq[0 + iYear*nSector], rsk_eq[1 + iYear*nSector]);
			//printf("%d\n", emp_eq[0 + iYear*nChoices]);

		} // End iYear loop


		for (iSec = 1; iSec < nSector + 1; iSec++)
		{
			for (iYear = 0; iYear < nYear - 1; iYear++)
			{
				avec[(iSec-1) + iYear*nSector] =
					rsk_eq[(iSec-1) + (iYear+1)*nSector] / rsk_eq[(iSec-1) + iYear*nSector];
			}


			avec[(iSec-1) + (nYear-1)*nSector] = 1.0;
		}


		// Write skill returns to file
		if (RE_loop == 1)
		{
			fp = fopen("Output/RSk_loop1.txt", "w");
			fprintf(fp, "Year\tRSk_1\tRSk_2\tRSk_3\tRSk_4\tRSk_5\n");
			for (iYear = 0; iYear < nYear; iYear++)
			{
				fprintf(fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\n", iYear+1996, rsk_eq[0 + iYear * nSector], rsk_eq[1 + iYear * nSector], rsk_eq[2 + iYear * nSector], rsk_eq[3 + iYear * nSector], rsk_eq[4 + iYear * nSector]);
			}
			fclose(fp);
		}
		else if (RE_loop == 2)
		{
			fp = fopen("Output/RSk_loop2.txt", "w");
			fprintf(fp, "Year\tRSk_1\tRSk_2\tRSk_3\tRSk_4\tRSk_5\n");
			for (iYear = 0; iYear < nYear; iYear++)
			{
				fprintf(fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\n", iYear+1996, rsk_eq[0 + iYear * nSector], rsk_eq[1 + iYear * nSector], rsk_eq[2 + iYear * nSector], rsk_eq[3 + iYear * nSector], rsk_eq[4 + iYear * nSector]);
			}
			fclose(fp);
		}

		RE_loop += 1;

	} // End RE while loop

	for (iYear = 0; iYear < nYear; iYear++)
	{
		for (iSec = 1; iSec < nSector + 1; iSec++)
		{
			z[(iSec-1) + iYear * nSector] = OutputData[iYear + (iSec-1) * nYear] /
				( pow(sk_sup_eq[(iSec-1) + iYear*nSector],IncomeShares[iYear + (iSec-1)*nYear + 0 * nSector * nYear])
				* pow(CapitalData[(iSec-1) + iYear * nSector],1.0-IncomeShares[iYear + (iSec-1)*nYear + 0 * nSector * nYear]) );
		}
	}


	// Use equilibrium value for skill returns as initial guess in next iteration
	for (iSec = 1; iSec < nSector + 1; iSec++)
	{
		rsk_global[iSec-1] = fmin( fmax( rsk_eq[(iSec-1) + 0 *nSector] , Low_RSk[iSec-1] ) , Up_RSk[iSec-1] );
	}


	//--------------------------------------------------------------------
	// Construct dataset for use in auxiliary regressions
	//--------------------------------------------------------------------

	// 'choice' counts how many make each choice for use in allocating memory
	// to GSL matrices and vectors for regressions.
	choice[0] = 0;
	choice[1] = 0;
	choice[2] = 0;
	choice[3] = 0;
	choice[4] = 0;
	choice[5] = 0;

	choiceWgt[0] = 0;
	choiceWgt[1] = 0;
	choiceWgt[2] = 0;
	choiceWgt[3] = 0;
	choiceWgt[4] = 0;
	choiceWgt[5] = 0;

	iCount = 0;
	for (iYear = 0; iYear < nYear; iYear++)
	{
		for (iCoh = 0; iCoh < nCoh; iCoh++)
		{
			age = (iYear + 1996) - (iCoh + 1931);

			for (iN = 0; iN < nPop; iN++)
			{

				if (age <= LastAge && age >= FirstAge)
				{
					if (ChoiceData[iPopCohYear]==0)
					{
						choice[0] += 1;
						choiceWgt[0] += CohortWgt[iCoh + iYear * nCoh];
					}
					else if (ChoiceData[iPopCohYear]==1)
					{
						choice[1] += 1;
						choiceWgt[1] += CohortWgt[iCoh + iYear * nCoh];
					}
					else if (ChoiceData[iPopCohYear]==2)
					{
						choice[2] += 1;
						choiceWgt[2] += CohortWgt[iCoh + iYear * nCoh];
					}
					else if (ChoiceData[iPopCohYear]==3)
					{
						choice[3] += 1;
						choiceWgt[3] += CohortWgt[iCoh + iYear * nCoh];
					}
					else if (ChoiceData[iPopCohYear]==4)
					{
						choice[4] += 1;
						choiceWgt[4] += CohortWgt[iCoh + iYear * nCoh];
					}
					else if (ChoiceData[iPopCohYear]==5)
					{
						choice[5] += 1;
						choiceWgt[5] += CohortWgt[iCoh + iYear * nCoh];
					}

					Data[iCount + 0 * size]  = iN + 1;
					Data[iCount + 1 * size]  = iCoh + 1931;
					Data[iCount + 2 * size]  = iYear + 1996;
					Data[iCount + 3 * size]  = FemaleData[iPopCoh];
					Data[iCount + 4 * size]  = EducData[iPopCoh];
					Data[iCount + 5 * size]  = (iYear + 1996) - (iCoh + 1931);
					Data[iCount + 6 * size]  = ((iYear + 1996) - (iCoh + 1931)) * ((iYear + 1996) - (iCoh + 1931));
					Data[iCount + 7 * size]  = ExperData[iPopCohYear];
					Data[iCount + 8 * size]  = ChoiceData[iPopCohYear];
					Data[iCount + 9 * size]  = LagChoiceData[iPopCohYear];
					Data[iCount + 10 * size] = log( WageData[iPopCohYear] );
					Data[iCount + 11 * size] = sqrt( CohortWgt[iCoh + iYear * nCoh] );

					iCount += 1;
				} // end if

			}
		}
	}


	// Write generated data to file
	fp = fopen("Output/Data.txt", "w");
	fprintf(fp, "N\tCohort\tYear\tFemale\tEduc\tElig\tType\tSector\tLagSec\tUnempYrs\tExper\tWage\tLagWage\tCost\tValFun\tCohortWgt\n");
	for (iCoh = 0; iCoh < nCoh; iCoh++)
	{
		for (iN = 0; iN < nPop; iN++)
		{
			for (iYear = 0; iYear < nYear; iYear++)
			{
				age = (iYear + 1996) - (iCoh + 1931);

				if (age >= FirstAge && age <= LastAge)
				{
					fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
						iN, iCoh+1931, iYear+1996, FemaleData[iPopCoh], EducData[iPopCoh], EligData[iPopCohYear],
						type[iPopCoh], ChoiceData[iPopCohYear], LagChoiceData[iPopCohYear], UnempYearsData[iPopCohYear],
						ExperData[iPopCohYear], WageData[iPopCohYear],
						LagWageData[iPopCohYear], CostSim[iPopCohYear], VSim[iPopCohYear], CohortWgt[iCoh + iYear * nCoh]);
				}
			}
		}
	}
	fclose(fp);


	fp = fopen("Output/InitSimulation.txt", "w");
	fprintf(fp, "N\tCohort\tYear\tFemale\tEduc\tElig\tSector\tLagSec\tUnempYrs\tExper\tLagWage\tCohortWgt\n");
	for (iCoh = LastYr - LastAge - 1931; iCoh < LastYr - FirstAge - 1931 + 1; iCoh++)
	{
		iYear = nYear - 1;

		for (iN = 0; iN < nPop; iN++)
		{
			fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\n",
				iN, iCoh+1931, iYear+1996, FemaleData[iPopCoh], EducData[iPopCoh], EligData[iPopCohYear],
				ChoiceData[iPopCohYear], LagChoiceData[iPopCohYear], UnempYearsData[iPopCohYear],
				ExperData[iPopCohYear], LagWageData[iPopCohYear], CohortWgt[iCoh + iYear * nCoh]);
		}
	}
	fclose(fp);


	fp = fopen("Output/Outcomes.txt", "w");
	fprintf(fp, "Year\tZ1\tZ2\tZ3\tZ4\tZ5\tRSk_eq1\tRSk_eq2\tRSk_eq3\tRSk_eq4\tRSk_eq5\tSk_sup_eq1\tSk_sup_eq2\tSk_sup_eq3\tSk_sup_eq4\tSk_sup_eq5\n");
	for (iYear = 0; iYear < nYear; iYear++)
	{
		fprintf(fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
			iYear+1996, z[0 + iYear * nSector], z[1 + iYear * nSector], z[2 + iYear * nSector], z[3 + iYear * nSector], z[4 + iYear * nSector],
			rsk_eq[0 + iYear * nSector], rsk_eq[1 + iYear * nSector], rsk_eq[2 + iYear * nSector], rsk_eq[3 + iYear * nSector], rsk_eq[4 + iYear * nSector],
			sk_sup_eq[0 + iYear * nSector], sk_sup_eq[1 + iYear * nSector], sk_sup_eq[2 + iYear * nSector], sk_sup_eq[3 + iYear * nSector], sk_sup_eq[4 + iYear * nSector]);
	}
	fclose(fp);




	//--------------------------------------------------------------------
	// Auxiliary regressions on simulated data
	//--------------------------------------------------------------------

	// Year dummies
	for (iYear = 0; iYear < nYear; iYear++)
	{
		for (iCount = 0; iCount < size; iCount++)
		{
			yearDummies[iCount + iYear * size] = (Data[iCount + 2 * size] == iYear + 1996);
		}
	}

	//----------------------------------------------------------------
	// Wage regressions
	//----------------------------------------------------------------

	for (iSec = 1; iSec < nSector + 1; iSec++)
	{

		iCount = 0;
		for (iSize = 0; iSize < size; iSize++)
		{

			if (Data[iSize + 8 * size]==iSec)
			{
				// Fill regression matrix and vector
				XWage[iCount + 0  * choice[iSec]] = Data[iSize + 3  * size]        * Data[iSize + 11 * size];
				XWage[iCount + 1  * choice[iSec]] = Data[iSize + 4  * size]        * Data[iSize + 11 * size];
				XWage[iCount + 2  * choice[iSec]] = Data[iSize + 5  * size]        * Data[iSize + 11 * size];
				XWage[iCount + 3  * choice[iSec]] = Data[iSize + 6  * size]        * Data[iSize + 11 * size];
				XWage[iCount + 4  * choice[iSec]] = Data[iSize + 7  * size]        * Data[iSize + 11 * size];
				XWage[iCount + 5  * choice[iSec]] = yearDummies[iSize + 0  * size] * Data[iSize + 11 * size];
				XWage[iCount + 6  * choice[iSec]] = yearDummies[iSize + 1  * size] * Data[iSize + 11 * size];
				XWage[iCount + 7  * choice[iSec]] = yearDummies[iSize + 2  * size] * Data[iSize + 11 * size];
				XWage[iCount + 8  * choice[iSec]] = yearDummies[iSize + 3  * size] * Data[iSize + 11 * size];
				XWage[iCount + 9  * choice[iSec]] = yearDummies[iSize + 4  * size] * Data[iSize + 11 * size];
				XWage[iCount + 10 * choice[iSec]] = yearDummies[iSize + 5  * size] * Data[iSize + 11 * size];
				XWage[iCount + 11 * choice[iSec]] = yearDummies[iSize + 6  * size] * Data[iSize + 11 * size];
				XWage[iCount + 12 * choice[iSec]] = yearDummies[iSize + 7  * size] * Data[iSize + 11 * size];
				XWage[iCount + 13 * choice[iSec]] = yearDummies[iSize + 8  * size] * Data[iSize + 11 * size];
				XWage[iCount + 14 * choice[iSec]] = yearDummies[iSize + 9  * size] * Data[iSize + 11 * size];
				XWage[iCount + 15 * choice[iSec]] = yearDummies[iSize + 10 * size] * Data[iSize + 11 * size];
				XWage[iCount + 16 * choice[iSec]] = yearDummies[iSize + 11 * size] * Data[iSize + 11 * size];
				XWage[iCount + 17 * choice[iSec]] = yearDummies[iSize + 12 * size] * Data[iSize + 11 * size];

				Y[iCount] = Data[iSize + 10 * size] * Data[iSize + 11 * size];
				iCount += 1;
			}

		}

		// OLS
		fcn_linreg(Y, XWage, choice[iSec], nRegWage, coef, &SST, &SSE);

		// Store coefficients
		for (iReg = 0; iReg < nRegWage; iReg++)
		{
			auxSim[iReg + (iSec-1) * nRegWage] = coef[iReg];
			// printf("n: %d, coef: %9.6f\n", iReg, coef[iReg]);
		}

		// Variance of residuals (RMSE)
		auxSim[(iSec-1) + nRegWage*nSector] = sqrt( SSE / choiceWgt[iSec] );

	}


	//----------------------------------------------------------------
	// LPM for choice (X is also used in transition regressions)
	//----------------------------------------------------------------
	for (iSize = 0; iSize < size; iSize++)
	{
		XOther[iSize + 0  * size] = Data[iSize + 10 * size]        * Data[iSize + 11 * size];
		XOther[iSize + 1  * size] = Data[iSize + 3  * size]        * Data[iSize + 11 * size];
		XOther[iSize + 2  * size] = Data[iSize + 4  * size]        * Data[iSize + 11 * size];
		XOther[iSize + 3  * size] = Data[iSize + 5  * size]        * Data[iSize + 11 * size];
		XOther[iSize + 4  * size] = Data[iSize + 6  * size]        * Data[iSize + 11 * size];
		XOther[iSize + 5  * size] = Data[iSize + 7  * size]        * Data[iSize + 11 * size];
		XOther[iSize + 6  * size] = yearDummies[iSize + 0  * size] * Data[iSize + 11 * size];
		XOther[iSize + 7  * size] = yearDummies[iSize + 1  * size] * Data[iSize + 11 * size];
		XOther[iSize + 8  * size] = yearDummies[iSize + 2  * size] * Data[iSize + 11 * size];
		XOther[iSize + 9  * size] = yearDummies[iSize + 3  * size] * Data[iSize + 11 * size];
		XOther[iSize + 10 * size] = yearDummies[iSize + 4  * size] * Data[iSize + 11 * size];
		XOther[iSize + 11 * size] = yearDummies[iSize + 5  * size] * Data[iSize + 11 * size];
		XOther[iSize + 12 * size] = yearDummies[iSize + 6  * size] * Data[iSize + 11 * size];
		XOther[iSize + 13 * size] = yearDummies[iSize + 7  * size] * Data[iSize + 11 * size];
		XOther[iSize + 14 * size] = yearDummies[iSize + 8  * size] * Data[iSize + 11 * size];
		XOther[iSize + 15 * size] = yearDummies[iSize + 9  * size] * Data[iSize + 11 * size];
		XOther[iSize + 16 * size] = yearDummies[iSize + 10 * size] * Data[iSize + 11 * size];
		XOther[iSize + 17 * size] = yearDummies[iSize + 11 * size] * Data[iSize + 11 * size];
		XOther[iSize + 18 * size] = yearDummies[iSize + 12 * size] * Data[iSize + 11 * size];
	}

	for (iSec = 1; iSec < nSector + 1; iSec++)
	{
		// Dep. var.: dummy for choice
		for (iSize = 0; iSize < size; iSize++)
		{
			Y[iSize] = ( Data[iSize + 8 * size] == iSec ) * Data[iSize + 11 * size];
		}

		// OLS
		fcn_linreg(Y, XOther, size, nRegOther, coef, &SST, &SSE);

		// Store coefficients
		for (iReg = 0; iReg < nRegOther; iReg++)
		{
			auxSim[iReg + (iSec-1) * nRegOther + nRegWage*nSector + nSector] = coef[iReg];
			//printf("%5d %5d %9.4f\n", iSec, iReg, coef[iReg]);
		}
	}

	//----------------------------------------------------------------
	// LPM for transition
	//----------------------------------------------------------------
	iCount = 0;
	for (iLagCh = 1; iLagCh < nSector + 1; iLagCh++)
	{
		for (iSec = 1; iSec < nSector + 1; iSec++)
		{
			// Dep. var.: transition dummy
			for (iSize = 0; iSize < size; iSize++)
			{
				Y[iSize] = (Data[iSize + 9 * size]==iLagCh && Data[iSize + 8 * size]==iSec) * Data[iSize + 11 * size];
			}

			// OLS
			fcn_linreg(Y, XOther, size, nRegOther, coef, &SST, &SSE);

			// Store coefficients
			for (iReg = 0; iReg < nRegOther; iReg++)
			{
				auxSim[iReg + iCount * nRegOther + nRegWage*nSector + nSector + nRegOther*nSector] = coef[iReg];
			}

			iCount += 1;
		}
	}

	//----------------------------------------------------------------
	// Loss function
	//----------------------------------------------------------------
	Loss = 0.0;

	// Use inverse covariance matrix for weighting
	for (iPar = 0; iPar < nParamAux; iPar++)
	{
		for (iPar2 = 0; iPar2 < nParamAux; iPar2++)
		{
			Loss += ( auxSim[iPar] - auxData[iPar] ) * invCov[iPar + iPar2 * nParamAux] *
					( auxSim[iPar2] - auxData[iPar2] );
		}
	}

	// Restrict nu and sigma to be positive
	if (nu[0] <= 0)
	{
		Loss = 999999999999.9;
	}
	else if (Find_Min(sigma, nSector) <= 0)
	{
		Loss = 999999999999.9;
	}

	//----------------------------------------------------------------
	// Write coef to file
	//----------------------------------------------------------------
	fp = fopen("Output/CompareCoefficients.txt", "w");

	for (iPar = 0; iPar < nParamAux; iPar++)
	{
		fprintf(fp, "%d\t%lf\t%lf\t%lf\n", iPar+1, auxData[iPar], auxSim[iPar], sqrt( covData[iPar + iPar * nParamAux] ) );
	}

	fclose(fp);


	//----------------------------------------------------------------
	// Print some stats to screen
	//----------------------------------------------------------------
	printf("Min and Max R_sq SE and RE:\n");
	printf("%f\t%f\n", Find_Min(R_sq_SE,(nAge-1)*nChoices*CatGender*CatEduc*CatElig), Find_Max(R_sq_SE,(nAge-1)*nChoices*CatGender*CatEduc*CatElig));
	printf("%f\t%f\n\n", Find_Min(R_sq_RE,(nYear-1)*(nAge-1)*nChoices*CatGender*CatEduc*CatElig), Find_Max(R_sq_RE,(nYear-1)*(nAge-1)*nChoices*CatGender*CatEduc*CatElig));

	printf("Loss:\t\t%f\n", Loss);


	//----------------------------------------------------------------
	// Free allocated memory
	//----------------------------------------------------------------
	free(Low_RSk);
	free(Up_RSk);
	free(PI_RE);
	free(R_sq_RE);
	free(PI_SE);
	free(R_sq_SE);
	free(avec);
	free(rsk_sup);
	free(rsk_dem);
	free(rsk_eq);
	free(sk_sup);
	free(sk_sup_eq);
	free(emp);
	free(emp_eq);
	free(iter);
	free(step);
	free(flag);
	free(choice);
	free(choiceWgt);
	free(yearDummies);
	free(beta);
	free(sigma);
	free(kappa);
	free(xi);
	free(gamma);
    free(delta);
	free(nu);
	free(lambda);
	free(phi);
	free(rsk_tomorrow);
	free(PI_Coef);
	free(w);
	free(V);
	free(Cost);
	free(check);
	free(z);
	free(Data);
	free(Female);
	free(Educ);
	free(Elig);
	free(type);
	free(Type);

	free(auxSim);

	free(XWage);
	free(XOther);
	free(Y);
	free(coef);

	// Free allocated memory for random number generator
	gsl_rng_free (r);

	return Loss;
} // end LossFcn



//------------------------------------------------------------------------------
// Module calling the SMD loss function.
//------------------------------------------------------------------------------
double ObjFcn(unsigned n, const double *in_pars)
{
	//----------------------------------------------------------
	// Initializations
	//----------------------------------------------------------

	// Counters
	int iParam, k;

	// Arrays
	double *param, loss;

	// File pointers
	FILE *fpPars;


	//----------------------------------------------------------
	// Memory Allocation
	//----------------------------------------------------------
	param = calloc(nParam,sizeof(double));


	//----------------------------------------------------------
	// Increment function evaluations counter
	//----------------------------------------------------------
	++count;

	//----------------------------------------------------------
	// Scale and print parameters
	//----------------------------------------------------------

	// Get un-masked parameters (if masked get starting value of parameter)
	k = 0;
	for (iParam = 0; iParam < nParam; iParam++)
	{
		if (mask[iParam]==1)
		{
			param[iParam] = in_pars[k];
			k += 1;
		}
		else if (mask[iParam]==0)
		{
			param[iParam] = param_start[iParam];
		}
	}

	// Scale back parameters and print
	for (iParam = 0; iParam < nParam; iParam++)
	{
		param[iParam] = param[iParam] * scaling[iParam];
	}

	printf("----------------------------------------------------------------\n");
	printf("Calls to ObjFcn: %d\n", count);
	printf("\n");
	printf("   beta1     beta2    beta3    beta4 \n");
	printf("%9.6f %9.6f %9.6f %9.6f\n", param[0],  param[1],  param[2],  param[3]);
	printf("%9.6f %9.6f %9.6f %9.6f\n", param[4],  param[5],  param[6],  param[7]);
	printf("%9.6f %9.6f %9.6f %9.6f\n", param[8],  param[9],  param[10], param[11]);
	printf("%9.6f %9.6f %9.6f %9.6f\n", param[12], param[13], param[14], param[15]);
	printf("%9.6f %9.6f %9.6f %9.6f\n", param[16], param[17], param[18], param[19]);
	printf("   sigma1    sigma2    sigma3    sigma4    sigma5\n");
	printf("%9.6f %9.6f %9.6f %9.6f %9.6f\n", param[20], param[21], param[22], param[23], param[24]);
	printf("   gamma0    gamma1    gamma2    gamma3    gamma4    gamma5\n");
	printf("%9.6f %9.6f %9.6f %9.6f %9.6f %9.6f\n", param[25], param[26], param[27], param[28], param[29], param[30]);
	printf("   xiOut2    xiOut3    xiOut4    xiOut5\n");
	printf("%9.6f %9.6f %9.6f %9.6f\n", param[31], param[32], param[33], param[34]);
	printf("   xiIn1     xiIn2     xiIn3     xiIn4     xiIn5\n");
	printf("%9.6f %9.6f %9.6f %9.6f %9.6f\n", param[35], param[36], param[37], param[38], param[39]);
	printf("   kappa\n");
	printf("%9.6f %9.6f %9.6f %9.6f\n", param[40], param[41], param[42], param[43]);
	printf("   nu\n");
	printf("%9.6f\n", param[44]);
	printf("  lambda1   lambda2   lambda3   lambda4   lambda5\n");
	printf("%9.6f %9.6f %9.6f %9.6f %9.6f\n", param[45], param[46], param[47], param[48], param[49]);
	printf("    phi1      phi2      phi3      phi4 \n");
	printf("%9.6f %9.6f %9.6f %9.6f\n\n", param[50],  param[51],  param[52],  param[53]);


	// Write current parameters to file
	fpPars = fopen("Output/Current_Param.txt", "w");
	for (iParam = 0; iParam < nParam; iParam++)
	{
		fprintf(fpPars, "%lf\n", param[iParam]);
	}
	fclose(fpPars);


	//----------------------------------------------------------
	// Compute indirect inference loss function
	//----------------------------------------------------------
	loss = LossFcn(param);

	//----------------------------------------------------------
	// Free momeory
	//----------------------------------------------------------
	free(param);


	return loss;
}

//------------------------------------------------------------------------------
// Simulation5.c
//
// By Damoun Ashournia, University of Oxford
//
//------------------------------------------------------------------------------
// Simulation5():  Both shocks, half UI period.
//------------------------------------------------------------------------------
int Simulation5()
{
	//----------------------------------------------------------
	// Initializations
	//----------------------------------------------------------

	// Counters
	int iCnt, iN, iCoh, iYear, iAge, iSec, iEduc,
		iAv, iReg, iIter, iRsk, it, it2, converged, converged2, iSec2;

	// Reading InitSimulation and Outcomes
	int age, coh, female, educ, elig, typ, choice, lagSec, unempyrs, RE_loop, year;
	double cohwgt, z1, z2, z3, z4, z5, rsk_eq1, rsk_eq2, rsk_eq3, rsk_eq4, rsk_eq5, sk_sup_eq1, sk_sup_eq2, sk_sup_eq3, sk_sup_eq4, sk_sup_eq5, lastwage,
		eps_aux, *eta_aux, benefit, Emax, Emax0, Punemp, VMAX, *Pindex, OutSum, exper, ExperTomorrow;
	

	// Arrays
	int *EducSim, *EducInit, *FemaleSim, *FemaleInit, *ChoiceSim, *ChoiceInit,
		*EligSim, *EligInit, *LagSectorSim, *LagSectorInit, *UnempYrsSim,
		*UnempYrsInit;
	
	double *alpha_cons, *alpha_prod, *CohortWgt, *CapitalSim,
		*WageSim, *LagWageSim, *LagWageInit, *rsk_tomorrow, *WelfareSim,
		*Cap, *Output, *ExperSim, *ExperInit;

	double *z, *rsk_eq, *sk_sup, *sk_sup_eq, *price, *prod, *emp_eq,
			*rsk_sup, *rsk_dem, *emp, *Welfare, *Welfare_Educ, *Real_Wages,
			*Real_Wages_Educ, *PI_Coef, *w, *V, *Cost, *Benefits;

	double *Low_RSk, *Up_RSk, *avec, Low_Wage, Up_Wage;

	double *rK_eq, *rK_sup, *rK_dem, *p_sup, *p_dem;

	// Parameters
	double *beta, *sigma, *gamma, *kappa, *xi, *nu, *lambda, *phi;

	// Convergence criteria
	double crit1, crit2, crit3, check11, check12, check13, check14, check15, check21, check22, check23, check24, check25, check31, check32, check33;
	double *step, *step2, *step3;
	int *iter, *iter2, *flag, *flag2, *flag3;

	// Cohort expansion for simulation
	int SimExpansion, n2;

	// Parallelization
	int *Female, *Educ, *Elig, *Type, iParallel, tid, nthreads;

	// File pointer
	FILE *fp;

	// Emax
	double *PI_RE, *R_sq_RE, *PI_SE, *R_sq_SE;

	// GSL Random Number Generator declarations
	double *epsSim, *etaSim, *unif_type, term2, prob;
	int unemp;

	const gsl_rng_type *T;
	gsl_rng *r;

	T = gsl_rng_mt19937; // Algorithm
	r = gsl_rng_alloc(T); // Allocate memory
	gsl_rng_set(r,5566778); // Set seed
	
	const gsl_rng_type *Q;
	gsl_rng *s;

	Q = gsl_rng_mt19937; // Algorithm
	s = gsl_rng_alloc(Q); // Allocate memory
	gsl_rng_set(s,123456); // Set seed

	//----------------------------------------------------------
	// Memory Allocation
	//----------------------------------------------------------
	EducSim       = calloc(nPop_sim*nCoh_sim,sizeof(int));
	EducInit      = calloc(nPop*nCoh_init,sizeof(int));
	EligSim       = calloc(nPop_sim*nCoh_sim,sizeof(int));
	EligInit      = calloc(nPop*nCoh_init,sizeof(int));
	FemaleSim     = calloc(nPop_sim*nCoh_sim,sizeof(int));
	FemaleInit    = calloc(nPop*nCoh_init,sizeof(int));

	ChoiceSim     = calloc(nPop_sim*nCoh_sim*nAge,sizeof(int));
	ChoiceInit    = calloc(nPop*nCoh_init*nAge,sizeof(int));
	LagSectorSim  = calloc(nPop_sim*nCoh_sim*nAge,sizeof(int));
	LagSectorInit = calloc(nPop*nCoh_init*nAge,sizeof(int));
	UnempYrsSim   = calloc(nPop_sim*nCoh_sim*nAge,sizeof(int));
	UnempYrsInit  = calloc(nPop*nCoh_init*nAge,sizeof(int));
	ExperSim      = calloc(nPop_sim*nCoh_sim*nAge,sizeof(double));
	ExperInit     = calloc(nPop*nCoh_init*nAge,sizeof(double));
	LagWageSim    = calloc(nPop_sim*nCoh_sim*nAge,sizeof(double));
	LagWageInit   = calloc(nPop*nCoh_init*nAge,sizeof(double));
	WageSim       = calloc(nPop_sim*nCoh_sim*nAge,sizeof(double));
	WelfareSim    = calloc(nPop_sim*nCoh_sim*nAge,sizeof(double));
	CapitalSim    = calloc(nSector*(LastYr_sim-FirstYr+1),sizeof(double));

	CohortWgt = calloc(nCoh_sim,sizeof(double));

	z               = calloc(nSector*(LastYr_sim-FirstYr+1),sizeof(double));
	rsk_eq          = calloc(nSector*(LastYr_sim-FirstYr+1),sizeof(double));
	sk_sup          = calloc(nSector,sizeof(double));
	sk_sup_eq       = calloc(nSector*(LastYr_sim-FirstYr+1),sizeof(double));
	price           = calloc(nSector*(LastYr_sim-FirstYr+1),sizeof(double));
	prod            = calloc(nSector*(LastYr_sim-FirstYr+1),sizeof(double));
	emp_eq          = calloc((nSector+1)*(LastYr_sim-FirstYr+1),sizeof(double));
	emp             = calloc((nSector+1)*MaxIt,sizeof(double));
	Welfare         = calloc((LastYr_sim-FirstYr_sim+1),sizeof(double));
	Benefits        = calloc((LastYr_sim-FirstYr_sim+1),sizeof(double));
	Welfare_Educ    = calloc(CatEduc*(LastYr_sim-FirstYr_sim+1),sizeof(double));
	Real_Wages      = calloc((LastYr_sim-FirstYr_sim+1),sizeof(double));
	Real_Wages_Educ = calloc(CatEduc*(LastYr_sim-FirstYr_sim+1),sizeof(double));
	iter            = calloc((LastYr_sim-FirstYr+1),sizeof(int));
	iter2           = calloc((LastYr_sim-FirstYr+1),sizeof(int));
	rsk_tomorrow    = calloc(nSector,sizeof(double));

	rK_eq  = calloc(nSector*(LastYr_sim-FirstYr+1),sizeof(double));
	rK_sup = calloc(nSector*MaxIt,sizeof(double));
	rK_dem = calloc(nSector*MaxIt,sizeof(double));
	p_sup  = calloc(MaxIt*3,sizeof(double)); // 3: number of non traded sectors
	p_dem  = calloc(MaxIt*3,sizeof(double));

	alpha_cons = calloc(nSector,sizeof(double));
	alpha_prod = calloc(nSector*(LastYr_sim-FirstYr+1),sizeof(double));

	Low_RSk = calloc(nSector,sizeof(double));
	Up_RSk  = calloc(nSector,sizeof(double));
	avec    = calloc(nSector*(LastYr_sim-FirstYr+1),sizeof(double));
	rsk_sup = calloc(nSector*MaxIt,sizeof(double));
	rsk_dem = calloc(nSector*MaxIt,sizeof(double));

	beta   = calloc(nSector*nParamBeta,sizeof(double));
	kappa  = calloc(4,sizeof(double));
	xi     = calloc(2*nSector,sizeof(double));
	sigma  = calloc(nSector,sizeof(double));
	gamma  = calloc(nSector+1,sizeof(double));
	nu     = calloc(1,sizeof(double));
	lambda = calloc(nSector,sizeof(double));
    phi    = calloc(4,sizeof(double));

	Female = calloc(CatGender*CatEduc*CatElig*CatTypes,sizeof(double));
	Educ   = calloc(CatGender*CatEduc*CatElig*CatTypes,sizeof(double));
	Elig   = calloc(CatGender*CatEduc*CatElig*CatTypes,sizeof(double));
	Type   = calloc(CatGender*CatEduc*CatElig*CatTypes,sizeof(double));

	PI_RE   = calloc(nReg*(LastYr_sim-(Yr_shock+1)-1)*(nAge-1)*(nSector+1)*CatGender*CatEduc*CatElig*CatTypes,sizeof(double));
	R_sq_RE = calloc((LastYr_sim-(Yr_shock+1)-1)*(nAge-1)*(nSector+1)*CatGender*CatEduc*CatElig*CatTypes,sizeof(double));
	PI_SE   = calloc(nReg*(nAge-1)*(nSector+1)*CatGender*CatEduc*CatElig*CatTypes,sizeof(double));
	R_sq_SE = calloc((nAge-1)*(nSector+1)*CatGender*CatEduc*CatElig*CatTypes,sizeof(double));

	PI_Coef = calloc(nReg,sizeof(double));
	w       = calloc(nSector,sizeof(double));
	V       = calloc(nSector,sizeof(double));
	Cost    = calloc(nSector,sizeof(double));

	Pindex = calloc((LastYr_sim-FirstYr_sim+1),sizeof(double));
	Cap    = calloc(nSector,sizeof(double));
	Output = calloc(nSector*(LastYr_sim-FirstYr+1),sizeof(double));

	step  = calloc(nSector,sizeof(double));
	step2 = calloc(nSector,sizeof(double));
	step3 = calloc(3,sizeof(double));
	flag  = calloc(nSector,sizeof(int));
	flag2 = calloc(nSector,sizeof(int));
	flag3 = calloc(3,sizeof(int));

	// GSL Random Number Generator
	epsSim    = calloc(nPop_sim*nAge*nSector,sizeof(double));
	etaSim    = calloc(nPop_sim*nAge*nSector,sizeof(double));
	eta_aux   = calloc(nSector,sizeof(double));
	unif_type = calloc(nPop_sim*nCoh_sim,sizeof(double));

	// Initialize values
	benefit = 0.0;
	check11 = 999.9;
	check12 = 999.9;
	check13 = 999.9;
	check14 = 999.9;
	check15 = 999.9;
	check21 = 999.9;
	check22 = 999.9;
	check23 = 999.9;
	check24 = 999.9;
	check25 = 999.9;
	check31 = 999.9;
	check32 = 999.9;
	check33 = 999.9;


	//----------------------------------------------------------
	// Definitions
	//----------------------------------------------------------
	#define iInit		iN + iCoh * nPop
	#define iInitAge	iN + iCoh * nPop + iAge * nCoh_init * nPop

	#define iSimPopCoh		iN + iCoh * nPop_sim
	#define iSimPopCohAge	iN + iCoh * nPop_sim + iAge     * nCoh_sim * nPop_sim
	#define iPopCohAge		iN + iCoh * nPop_sim + (age-30) * nCoh_sim * nPop_sim
	#define iPopCohAgeNext	iN + iCoh * nPop_sim + (age-30+1) * nCoh_sim * nPop_sim

	#define iEps	iN + (age-30) * nPop_sim + (iSec-1) * nAge * nPop_sim
	#define iEta	iN + (age-30) * nPop_sim + (iSec-1) * nAge * nPop_sim

	#define iAvec2 	(iSec2-1) + iYear * nSector

	#define iPISE	iReg + (age-30)*nReg + lagSec*(nAge-1)*nReg + female*(nSector+1)*(nAge-1)*nReg + educ*CatGender*(nSector+1)*(nAge-1)*nReg + elig*CatEduc*CatGender*(nSector+1)*(nAge-1)*nReg + (typ-1)*CatElig*CatEduc*CatGender*(nSector+1)*(nAge-1)*nReg
	#define iPIRE	iReg + (iYear-163)*nReg + (age-30)*(LastYr_sim-(Yr_shock+1)-1)*nReg + lagSec*(nAge-1)*(LastYr_sim-(Yr_shock+1)-1)*nReg + female*(nSector+1)*(nAge-1)*(LastYr_sim-(Yr_shock+1)-1)*nReg + educ*CatGender*(nSector+1)*(nAge-1)*(LastYr_sim-(Yr_shock+1)-1)*nReg + elig*CatEduc*CatGender*(nSector+1)*(nAge-1)*(LastYr_sim-(Yr_shock+1)-1)*nReg + (typ-1)*CatElig*CatEduc*CatGender*(nSector+1)*(nAge-1)*(LastYr_sim-(Yr_shock+1)-1)*nReg

	#define iDelta  iYear + (age-30) * (LastYr_sim - FirstYr + 1) + female * nAge * (LastYr_sim - FirstYr + 1) + educ * CatGender * nAge * (LastYr_sim - FirstYr + 1) + lagSec * CatEduc * CatGender * nAge * (LastYr_sim - FirstYr + 1)

	//----------------------------------------------------------
	// Unpack paramters
	//----------------------------------------------------------
	unpackPars(param, beta, sigma, gamma, xi, kappa, nu, lambda, phi);

	//----------------------------------------------------------
	// Consumption Shares
	//----------------------------------------------------------
	alpha_cons[0] = 0.011447676;
	alpha_cons[1] = 0.175704772;
	alpha_cons[2] = 0.061006087;
	alpha_cons[3] = 0.274791681;
	alpha_cons[4] = 0.477049784;


	for (iSec = 0; iSec < nSector; iSec++)
	{
		for (iYear = 0; iYear < nYear; iYear++)
		{
			alpha_prod[iSec + iYear * nSector] = IncomeShares[iYear + iSec * nYear + 0 * nSector * nYear];
		}
	}


	//----------------------------------------------------------
	// Convergence criteria
	//----------------------------------------------------------
	crit1 = 0.00005;
	crit2 = 0.00005;
	crit3 = 0.00005;


	//----------------------------------------------------------
	// Read initial conditions generated by estimation
	//----------------------------------------------------------
	fp = fopen("Data/InitSimulation.txt", "r");

	if (fp == NULL)
	{
	    perror("Failed to open file \"InitSimulation.txt\"");
	    exit(EXIT_FAILURE);
	}

	// Read first line (Variable names)
	char buffer[98];
	fgets(buffer, 98, fp);

	for (iCnt = 0; iCnt < nAge * nPop; iCnt++)
	{
		fscanf(fp, "%d %d %d %d %d %d %d %d %d %lf %lf %lf",
	   		&iN, &coh, &iYear, &female, &educ, &elig, &choice, &lagSec, &unempyrs, &exper, &lastwage, &cohwgt );

		iCoh = coh - FirstCoh_sim;
		iAge = iYear - coh - 30;

	   	FemaleInit[iInit]       = female;
	   	EducInit[iInit]			= educ;
	   	EligInit[iInit]			= elig;
	   	ChoiceInit[iInitAge]	= choice;
	   	LagSectorInit[iInitAge]	= lagSec;
	   	UnempYrsInit[iInitAge]	= unempyrs;
	   	ExperInit[iInitAge]     = exper;
	   	LagWageInit[iInitAge]   = lastwage;
	   	CohortWgt[iCoh]			= cohwgt;
	}

	fclose(fp);

	
	SimExpansion = (nPop_sim / nPop);

	for (iCoh = 0; iCoh < nCoh_sim; iCoh++)
	{
		if (iCoh + 1943 <= LastCoh)
		{
			CohortWgt[iCoh] = CohortWgt[iCoh] / SimExpansion;
		}
		else
		{
			CohortWgt[iCoh] = CohortWgt[LastCoh-1943];
		}
	}


	//----------------------------------------------------------
	// From FirstYr_sim + 1, all cohorts look the same as the
	// 1978 cohort in terms of gender, education, and
	// eligibility, and experience
	//----------------------------------------------------------
	for (iCoh = 0; iCoh < nCoh_sim; iCoh++)
	{
		if (iCoh + 1943 <= LastCoh)
		{
			iAge = 2008 - (iCoh + 1943) - 30;
			for (iN = 0; iN < nPop_sim; iN++)
			{
				n2 = (iN / SimExpansion);
				FemaleSim[iSimPopCoh]       = FemaleInit[n2 + iCoh * nPop];
				EducSim[iSimPopCoh]         = EducInit[n2 + iCoh * nPop];
				EligSim[iSimPopCoh]         = EligInit[n2 + iCoh * nPop];
				ChoiceSim[iSimPopCohAge]    = ChoiceInit[n2 + iCoh * nPop + iAge * nCoh_init * nPop];
				LagSectorSim[iSimPopCohAge] = LagSectorInit[n2 + iCoh * nPop + iAge * nCoh_init * nPop];
				UnempYrsSim[iSimPopCohAge]  = UnempYrsInit[n2 + iCoh * nPop + iAge * nCoh_init * nPop];
				ExperSim[iSimPopCohAge]     = ExperInit[n2 + iCoh * nPop + iAge * nCoh_init * nPop];
				LagWageSim[iSimPopCohAge]   = LagWageInit[n2 + iCoh * nPop + iAge * nCoh_init * nPop];
			}
		}
		else
		{
			iAge = 0;
			for (iN = 0; iN < nPop_sim; iN++)
			{
			n2 = (iN / SimExpansion);
			FemaleSim[iSimPopCoh]       = FemaleInit[n2 + 35 * nPop];
			EducSim[iSimPopCoh]         = EducInit[n2 + 35 * nPop];
			EligSim[iSimPopCoh]         = EligInit[n2 + 35 * nPop];
			ChoiceSim[iSimPopCohAge]    = ChoiceInit[n2 + 35 * nPop + iAge * nCoh_init * nPop];
			LagSectorSim[iSimPopCohAge] = LagSectorInit[n2 + 35 * nPop + iAge * nCoh_init * nPop];
			UnempYrsSim[iSimPopCohAge]  = UnempYrsInit[n2 + 35 * nPop + iAge * nCoh_init * nPop];
			ExperSim[iSimPopCohAge]     = ExperInit[n2 + 35 * nPop + iAge * nCoh_init * nPop];
			LagWageSim[iSimPopCohAge]   = LagWageInit[n2 + 35 * nPop + iAge * nCoh_init * nPop];
			}
		}
	}


	//-----------------------------------------------------------
	// Read equilibrium returns to skill, productivity generated
	// in estimation code at optimal structural parameters
	//-----------------------------------------------------------
	fp = fopen("Data/Outcomes.txt", "r");

	if (fp == NULL)
	{
	    perror("Failed to open file \"Outcomes.txt\"");
	    exit(EXIT_FAILURE);
	}

	// Read first line (Variable names)
	char buffer_outcomes[115];
	fgets(buffer_outcomes, 115, fp);

	for (iCnt = 0; iCnt < nYear; iCnt++)
	{
		fscanf(fp, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	   		&iYear, &z1, &z2, &z3, &z4, &z5, &rsk_eq1, &rsk_eq2, &rsk_eq3, &rsk_eq4, &rsk_eq5, &sk_sup_eq1, &sk_sup_eq2, &sk_sup_eq3, &sk_sup_eq4, &sk_sup_eq5 );

		z[0 + (iYear - 1996) * nSector]      = z1;
		z[1 + (iYear - 1996) * nSector]      = z2;
		z[2 + (iYear - 1996) * nSector]      = z3;
		z[3 + (iYear - 1996) * nSector]      = z4;
		z[4 + (iYear - 1996) * nSector]      = z5;
		rsk_eq[0 + (iYear - 1996) * nSector] = rsk_eq1;
		rsk_eq[1 + (iYear - 1996) * nSector] = rsk_eq2;
		rsk_eq[2 + (iYear - 1996) * nSector] = rsk_eq3;
		rsk_eq[3 + (iYear - 1996) * nSector] = rsk_eq4;
		rsk_eq[4 + (iYear - 1996) * nSector] = rsk_eq5;
		sk_sup_eq[0 + (iYear - 1996) * nSector] = sk_sup_eq1;
		sk_sup_eq[1 + (iYear - 1996) * nSector] = sk_sup_eq2;
		sk_sup_eq[2 + (iYear - 1996) * nSector] = sk_sup_eq3;
		sk_sup_eq[3 + (iYear - 1996) * nSector] = sk_sup_eq4;
		sk_sup_eq[4 + (iYear - 1996) * nSector] = sk_sup_eq5;
	}

	fclose(fp);

	for (iYear = 0; iYear < LastYr_sim - FirstYr + 1; iYear++)
	{
		for (iSec = 0; iSec < nSector; iSec++)
		{
			if (iYear < nYear)
			{
				CapitalSim[iSec + iYear * nSector] = CapitalData[iSec + iYear * nSector];
				rK_eq[iSec + iYear * nSector] = rKData[iSec + iYear * nSector];
			}
			else
			{
				CapitalSim[iSec + iYear * nSector] = CapitalData[iSec + (nYear-1) * nSector];
				rK_eq[iSec + iYear * nSector] = 0.0;
			}
		}
	}


	//----------------------------------------------------------
	// Price of all sectors are fixed until the shock
	//----------------------------------------------------------
	for (iYear = 0; iYear < nYear; iYear++)
	{
		for (iSec = 0; iSec < nSector; iSec++)
		{
			price[iSec + iYear * nSector] = 1.0;
			prod[iSec + iYear * nSector]  = z[iSec + iYear * nSector];
		}
	}

	for (iYear = nYear; iYear < LastYr_sim - FirstYr + 1; iYear++)
	{
		for (iSec = 0; iSec < nSector; iSec++)
		{
			price[iSec + iYear * nSector]      = 1.0;
			prod[iSec + iYear * nSector]       = z[iSec + (nYear-1) * nSector];
			alpha_prod[iSec + iYear * nSector] = alpha_prod[iSec + (nYear-1) * nSector];
		}
	}


	//----------------------------------------------------------
	// Price of sector 2 gets a once and for all shock
	//----------------------------------------------------------
	for (iYear = Yr_shock - FirstYr; iYear < LastYr_sim - FirstYr + 1; iYear++)
	{
		price[1 + iYear * nSector] = 0.9;
	}


	//----------------------------------------------------------
	// Maximum UI period is halved from 4 to 2 years
	//----------------------------------------------------------
	UImax = 2;

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
	for (iAv = 0; iAv < nSector*(LastYr_sim-FirstYr+1); iAv++)
	{
		avec[iAv] = 1.0;
	}

	// Type probabilities
	for (iCoh = 0; iCoh < nCoh_sim; iCoh++)
	{
		for (iN = 0; iN < nPop_sim; iN++)
		{
			unif_type[iSimPopCoh] = gsl_rng_uniform_pos(r);
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

	printf("\n");
	printf("---------------------------------------\n");
	printf(" Computing Emax SE\n");
	printf("---------------------------------------\n");
	printf("\n");
	// Fork a team of threads
	#pragma omp parallel private(tid, iParallel) shared(PI_SE, R_sq_SE)
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
			Emax_SE(tid, param, Female[iParallel], Educ[iParallel], Elig[iParallel], Type[iParallel], Low_RSk, Up_RSk, Low_Wage, Up_Wage, PI_SE, R_sq_SE);
		}
	}
	printf("\n");
	printf("---------------------------------------\n");
	printf(" Emax SE computed\n");
	printf("---------------------------------------\n");
	printf("\n");


	//--------------------------------------------------------------------
	// Solve model with rational expectations
	//--------------------------------------------------------------------
	RE_loop = 1;

	while (RE_loop == 1)//(RE_loop <= 2)
	{
		if (RE_loop >= 2)
		{

			printf("\n");
			printf("---------------------------------------\n");
			printf(" Computing Emax RE\n");
			printf("---------------------------------------\n");
			printf("\n");
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
			printf("\n");
			printf("---------------------------------------\n");
			printf(" Emax RE computed\n");
			printf("---------------------------------------\n");
			printf("\n");
		}

		//------------------------------------------------------
		// Initializing
		//------------------------------------------------------
		for (iYear = 0; iYear < (LastYr_sim-FirstYr+1); iYear++)
		{
			for (iSec = 0; iSec < nSector+1; iSec++)
			{
				emp_eq[iSec + iYear * (nSector+1)] = 0.0;
			}
			iter[iYear] = 0;
		}

		for (iSec = 0; iSec < nSector; iSec++)
		{
			for (iIter = 0; iIter < MaxIt; iIter++)
			{
				rsk_sup[iSec + iIter * nSector] = 0.0;
			}
			rsk_sup[iSec + 0 * nSector] = rsk_eq[iSec + 12 * nSector];
		}

		for (iSec = 0; iSec < 3; iSec++) {
			p_sup[0 + iSec * MaxIt] = 1.0;
			p_dem[0 + iSec * MaxIt] = 0.0;
		}

		
		//------------------------------------------------------
		// Simulation
		//------------------------------------------------------
		for (iYear = FirstYr_sim - FirstYr; iYear < LastYr_sim - FirstYr + 1; iYear++)//(iYear = FirstYr_sim - FirstYr; iYear < LastYr_sim - FirstYr + 1; iYear++)
		{
			year = iYear + FirstYr;
			printf("Year: %d\n", year);

			//--------------------------------------------------
			// Random draws
			//--------------------------------------------------
			for (iSec = 1; iSec < nSector + 1; iSec++)
			{
				for (age = 30; age < nAge + 30; age++)
				{
					for (iN = 0; iN < nPop_sim; iN++)
					{
						epsSim[iEps] = gsl_ran_gaussian(r, 1);
						etaSim[iEta] = gsl_rng_uniform_pos(r);
					}
				}
			}

			if (year > FirstYr_sim)
			{
				// Use last period result as initial point for rsk_sup and rK_sup
				for (iRsk = 0; iRsk < nSector*MaxIt; iRsk++)
				{
					rsk_sup[iRsk] = 0.0;
					rsk_dem[iRsk] = 0.0;
					rK_sup[iRsk] = 0.0;
					rK_dem[iRsk] = 0.0;
				}

				for (iSec = 0; iSec < nSector; iSec++)
				{
					rsk_sup[iSec + 0 * nSector] = rsk_eq[iSec + (iYear-1)*nSector];//rsk_eq[iSec + 12*nSector];
					rK_sup[iSec + 0 * nSector] = rK_eq[iSec + (iYear-1)*nSector];
				}

				for (iIter = 0; iIter < MaxIt*3; iIter++)
				{
					p_sup[iIter] = 0.0;
					p_dem[iIter] = 0.0;
				}
				p_sup[0 + 0 * MaxIt] = price[2 + (iYear-1) * nSector];
				p_sup[0 + 1 * MaxIt] = price[3 + (iYear-1) * nSector];
				p_sup[0 + 2 * MaxIt] = price[4 + (iYear-1) * nSector];
			}

			it = 0; // iteration counter
			converged = 0; // convergence indicator


			//---------------------------------------------------
			// Convergence for returns to skill and capital
			//---------------------------------------------------
			while (converged==0 && it < MaxIt)
			{
				// converged = 1;
				// Initialize to zero
				for (iSec = 0; iSec < nSector + 1; iSec++)
				{
					if (iSec > 0) {
						sk_sup[iSec-1] = 0.0;
					}
					emp[iSec + it * (nSector+1)] = 0.0;
				}

				Welfare[iYear]    = 0.0;
				Real_Wages[iYear] = 0.0;

				for (iEduc = 0; iEduc < CatEduc; iEduc++)
				{
					Welfare_Educ[iEduc + iYear * CatEduc]    = 0.0;
					Real_Wages_Educ[iEduc + iYear * CatEduc] = 0.0;
				}

				// Compute aggregate skill supply
				for (iCoh = 0; iCoh < nCoh_sim; iCoh++)
				{
					age      = (iYear + 2008 - 12) - (iCoh + 1943);
					if (age <= LastAge && age >= FirstAge)
					{
						for (iN = 0; iN < nPop_sim; iN++)
						{
							for (iSec = 1; iSec < nSector + 1; iSec++)
							{
								//printf("age: %d\n", age-30);
								female          = FemaleSim[iSimPopCoh];
								educ            = EducSim[iSimPopCoh];
								lagSec          = LagSectorSim[iPopCohAge];
								exper           = ExperSim[iPopCohAge];
								lastwage        = LagWageSim[iPopCohAge];
								eps_aux         = sigma[iSec-1] * epsSim[iEps];
								eta_aux[iSec-1] = -log(-log(etaSim[iEta]))*nu[0] - 0.577215665*nu[0];

								term2 = exp(phi[0] + phi[1]*female + phi[2]*educ + phi[3]*exper);
								prob  = 1.0 / (1.0 + term2);
								if (unif_type[iSimPopCoh] < prob)
								{
									typ = 1;
								}
								else
								{
									typ = 2;
								}
								

								// Eligibility for UI
								if (EligSim[iSimPopCoh] == 1 && UnempYrsSim[iPopCohAge] <= UImax)
								{
									elig = 1;
								}
								else
								{
									elig = 0;
								}
							
								if (age < LastAge)
								{
									if (iYear + 1996 >= Yr_shock && iYear + 1996 <= LastYr_sim-2 && RE_loop >= 2)
									{
										for (iSec2 = 1; iSec2 < nSector + 1; iSec2++)
										{
											rsk_tomorrow[iSec2-1] = avec[iAvec2] * rsk_sup[(iSec2-1) + it * nSector];
										}
									}
									else
									{
										for (iSec2 = 1; iSec2 < nSector + 1; iSec2++)
										{
											rsk_tomorrow[iSec2-1] = rsk_sup[(iSec2-1) + it * nSector];
										}
									}
										
									if (lagSec == iSec)
									{
										ExperTomorrow = exper + 1.0;
									}
									else
									{
										ExperTomorrow = gamma[iSec]*exper;
									}
									
									if (year >= Yr_shock && year <= LastYr_sim-3 && RE_loop >= 2)
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

									Emax  = Emax_Hat(PI_Coef,rsk_tomorrow,ExperTomorrow,lastwage);
									Emax0 = Emax_Hat(PI_Coef,rsk_tomorrow,gamma[0]*exper,lastwage);
								}
								else
								{
									Emax  = 0.0;
									Emax0 = 0.0;
								}

								if (lagSec == iSec)
								{
									// exper = exper;
								}
								else
								{
									exper = gamma[iSec]*exper;
								}

								// Wage and value function
								if (typ == 1)
								{
									w[iSec-1] = rsk_sup[(iSec-1) + it * nSector] * 
										exp(  beta[0+(iSec-1)*nParamBeta]*female + beta[1+(iSec-1)*nParamBeta]*educ
											+ beta[2+(iSec-1)*nParamBeta]*exper  + beta[3+(iSec-1)*nParamBeta]*exper*exper) * exp(eps_aux);
								}
								else if (typ == 2)
								{
									w[iSec-1] = rsk_sup[(iSec-1) + it * nSector] * 
										exp(  beta[0+(iSec-1)*nParamBeta]*female + beta[1+(iSec-1)*nParamBeta]*educ
											+ beta[2+(iSec-1)*nParamBeta]*exper  + beta[3+(iSec-1)*nParamBeta]*exper*exper + lambda[iSec-1]) * exp(eps_aux);
								}


								if (lagSec == iSec || lagSec == 0)
								{
									Cost[iSec-1] = 0.0;
								}
								else
								{
									Cost[iSec-1] = exp( xi[lagSec-1 + 0 * nSector] + xi[iSec-1 + 1 * nSector] + kappa[0]*female + 
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
					
								V[iSec-1] = (1 - Delta[iDelta])*(w[iSec-1] - Cost[iSec-1] + rho*Emax) + Delta[iDelta]*(benefit + rho*Emax0) + eta_aux[iSec-1];
								
							} // end iSec


							// Find choice (= maxV) and store next period states
							VMAX = Find_Max(V, nSector);
							

							// Unemployment shock
  							Punemp = gsl_rng_uniform_pos(r);

  							if (Punemp <= Delta[iDelta])
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
										ChoiceSim[iPopCohAge]         = iSec;
										WageSim[iPopCohAge]           = w[iSec-1];
										Real_Wages[iYear-12]         += w[iSec-1] * CohortWgt[iCoh];
										Welfare[iYear-12]            += (w[iSec-1] - Cost[iSec-1] + eta_aux[iSec-1]) * CohortWgt[iCoh];
										WelfareSim[iPopCohAge]        = (w[iSec-1] - Cost[iSec-1] + eta_aux[iSec-1]);
										emp[iSec + it * (nSector+1)] += 1.0 + CohortWgt[iCoh];
										sk_sup[iSec-1]               += CohortWgt[iCoh] * (w[iSec-1] / rsk_sup[(iSec-1) + it * nSector]);
									}
								}
  							}
  							else
  							{
  								ChoiceSim[iPopCohAge]      = 0;
  								WageSim[iPopCohAge]        = benefit;
  								Benefits[iYear-12]        += benefit * CohortWgt[iCoh];
  								WelfareSim[iPopCohAge]     = benefit;
  								emp[0 + it * (nSector+1)] += 1.0 + CohortWgt[iCoh];
  							}


							if (age < LastAge && year < LastYr_sim - 1)
							{
								if (ChoiceSim[iPopCohAge] == lagSec && lagSec > 0)
								{
									ExperSim[iPopCohAgeNext]    = exper + 1.0;
									LagWageSim[iPopCohAgeNext]  = WageSim[iPopCohAge];
									UnempYrsSim[iPopCohAgeNext] = 0;
								}
								else if (ChoiceSim[iPopCohAge] == 0)
								{
									ExperSim[iPopCohAgeNext]    = gamma[0]*exper;
									LagWageSim[iPopCohAgeNext]  = LagWageSim[iPopCohAge];
									UnempYrsSim[iPopCohAgeNext] = UnempYrsSim[iPopCohAge] + 1;
								}
								else
								{
									ExperSim[iPopCohAgeNext]    = gamma[ChoiceSim[iPopCohAge]]*exper;
									LagWageSim[iPopCohAgeNext]  = WageSim[iPopCohAge];
									UnempYrsSim[iPopCohAgeNext] = 0;
								}
								LagSectorSim[iPopCohAgeNext] = ChoiceSim[iPopCohAge];
							}
						
						} // end iN
					} // if age
				} // end iCoh

				// Price index
				Pindex[iYear-12] = pow(price[0 + iYear * nSector],alpha_cons[0])
								 * pow(price[1 + iYear * nSector],alpha_cons[1])
								 * pow(p_sup[it + 0 * MaxIt],alpha_cons[2])
								 * pow(p_sup[it + 1 * MaxIt],alpha_cons[3])
								 * pow(p_sup[it + 2 * MaxIt],alpha_cons[4]);

				z[0 + iYear * nSector] = price[0 + iYear * nSector] * prod[0 + iYear * nSector] * Pindex[0] / Pindex[iYear-12];
				z[1 + iYear * nSector] = price[1 + iYear * nSector] * prod[1 + iYear * nSector] * Pindex[0] / Pindex[iYear-12];
				z[2 + iYear * nSector] = p_sup[it + 0 * MaxIt] * prod[2 + iYear * nSector] * Pindex[0] / Pindex[iYear-12];
				z[3 + iYear * nSector] = p_sup[it + 1 * MaxIt] * prod[3 + iYear * nSector] * Pindex[0] / Pindex[iYear-12];
				z[4 + iYear * nSector] = p_sup[it + 2 * MaxIt] * prod[4 + iYear * nSector] * Pindex[0] / Pindex[iYear-12];

				//----------------------------------------------------
				// Compute allocation of capital and rK
				//----------------------------------------------------

				// Allow capital to grow with skill supply until steady
				// state. To do that, rK is fixed until SS is reached
				if (year < HalfYr_sim - 10)
				{
					for (iSec = 1; iSec < nSector + 1; iSec++)
					{
						CapitalSim[iSec-1 + iYear * nSector] = sk_sup[iSec-1] *  pow( ( (1-alpha_prod[iSec-1 + (nYear-1) * nSector]) * z[iSec-1 + iYear * nSector] / rKData[iSec-1 + (nYear-1) * nSector]  ), 1.0 / alpha_prod[iSec-1 + (nYear-1) * nSector]);
					}
				}

				//----------------------------------------------------
				// rK and capital that results from rsk_sup and p_sup
				//----------------------------------------------------

				if (year == FirstYr_sim)
				{
					rK_sup[0 + 0 * nSector] = rKData[0 + (nYear-1) * nSector];
					rK_sup[1 + 0 * nSector] = rKData[1 + (nYear-1) * nSector];
					rK_sup[2 + 0 * nSector] = rKData[2 + (nYear-1) * nSector];
					rK_sup[3 + 0 * nSector] = rKData[3 + (nYear-1) * nSector];
					rK_sup[4 + 0 * nSector] = rKData[4 + (nYear-1) * nSector];
				}
				else
				{
					rK_sup[0 + 0 * nSector] = rK_eq[0 + (iYear-1) * nSector];
					rK_sup[1 + 0 * nSector] = rK_eq[1 + (iYear-1) * nSector];
					rK_sup[2 + 0 * nSector] = rK_eq[2 + (iYear-1) * nSector];
					rK_sup[3 + 0 * nSector] = rK_eq[3 + (iYear-1) * nSector];
					rK_sup[4 + 0 * nSector] = rK_eq[4 + (iYear-1) * nSector];
					if (it > 1)
					{
						rK_sup[0 + 0 * nSector] = rK_sup[0 + (it2-1) * nSector];
						rK_sup[1 + 0 * nSector] = rK_sup[1 + (it2-1) * nSector];
						rK_sup[2 + 0 * nSector] = rK_sup[2 + (it2-1) * nSector];
						rK_sup[3 + 0 * nSector] = rK_sup[3 + (it2-1) * nSector];
						rK_sup[4 + 0 * nSector] = rK_sup[4 + (it2-1) * nSector];	
					}
				}

				it2 = 0;
				converged2 = 0;
				for (iSec = 1; iSec < nSector + 1; iSec++)
				{
					flag2[iSec-1] = 0;
				}


				while (converged2 == 0 && it2 < MaxIt)
				{
					for (iSec = 1; iSec < nSector + 1; iSec++)
					{
						Cap[iSec-1] = sk_sup[iSec-1] *  pow( ( (1-alpha_prod[iSec-1 + (nYear-1) * nSector]) * z[iSec-1 + iYear * nSector] / rK_sup[iSec-1 + it2 * nSector]  ), 1.0 / alpha_prod[iSec-1 + (nYear-1) * nSector]);

						rK_dem[iSec-1 + it2 * nSector] = (1 - alpha_prod[iSec-1 + (nYear-1) * nSector]) * z[iSec-1 + iYear * nSector] * 
							pow(sk_sup[iSec-1],alpha_prod[iSec-1 + (nYear-1) * nSector]) *
							pow(Cap[iSec-1],-alpha_prod[iSec-1 + (nYear-1) * nSector]);
					}

					check21 = fabs( log(rK_sup[0 + it2 * nSector]) - log(rK_dem[0 + it2 * nSector]) );
					check22 = fabs( log(rK_sup[1 + it2 * nSector]) - log(rK_dem[1 + it2 * nSector]) );
					check23 = fabs( log(rK_sup[2 + it2 * nSector]) - log(rK_dem[2 + it2 * nSector]) );
					check24 = fabs( log(rK_sup[3 + it2 * nSector]) - log(rK_dem[3 + it2 * nSector]) );
					check25 = fabs( log(rK_sup[4 + it2 * nSector]) - log(rK_dem[4 + it2 * nSector]) );

					if (check21 <= crit2 && check22 <= crit2 && check23 <= crit2 && check24 <= crit2 && check25 <= crit2)
					{
						converged2 = 1;
					}

					if (converged2 == 0 && it2 < MaxIt - 1)
					{
						if (it2 < 30)
						{
							for (iSec = 1; iSec < nSector + 1; iSec++)
							{
								step2[iSec-1] = 1.0 / 10.0;
							}
						}
						else
						{
							for (iSec = 1; iSec < nSector + 1; iSec++)
							{
								if (SameSign(1.0, rK_dem[iSec-1 + it2 * nSector]     - rK_dem[iSec-1 + (it2-1) * nSector]) !=
									SameSign(1.0, rK_dem[iSec-1 + (it2-1) * nSector] - rK_dem[iSec-1 + (it2-2) * nSector])
									&& flag2[iSec-1] == 0 )
								{
									flag2[iSec-1] = it2;
									step2[iSec-1] = step2[iSec-1] / 1.3;
								}

								if (SameSign(1.0, rK_dem[iSec-1 + it2 * nSector]     - rK_dem[iSec-1 + (it2-1) * nSector]) !=
									SameSign(1.0, rK_dem[iSec-1 + (it2-1) * nSector] - rK_dem[iSec-1 + (it2-2) * nSector])
									&& it2 >= flag2[iSec-1] + 5 )
								{
									flag2[iSec-1] = it2;
									step2[iSec-1] = step2[iSec-1] / 1.3;
								}
							}
						}

						// Update the return to capital
						for (iSec = 1; iSec < nSector + 1; iSec++)
						{
							rK_sup[iSec-1 + (it2+1) * nSector] = rK_sup[iSec-1 + it2 * nSector]
								- step2[iSec-1] * (rK_sup[iSec-1 + it2 * nSector] - rK_dem[iSec-1 + it2 * nSector]);
						}
					}
					else
					{
						// Record equilibrium quantities
						for (iSec = 1; iSec < nSector + 1; iSec++)
						{
							rK_eq[iSec-1 + iYear * nSector] = rK_dem[iSec-1 + it2 * nSector];
							iter2[iYear] = it2;
							// printf("it: %d, rK: %15.11f\n", iter2[iYear], rK_eq[iSec-1 + iYear * nSector]);
						}
					}

					it2 += 1;
				} // end convergence of rK_sup and rK_dem

				for (iSec = 1; iSec < nSector + 1; iSec++)
				{
					Output[iSec-1 + iYear * nSector] = z[iSec-1 + iYear * nSector] *
						pow(sk_sup[iSec-1],alpha_prod[iSec-1 + (nYear-1) * nSector]) *
						pow(Cap[iSec-1],1.0-alpha_prod[iSec-1 + (nYear-1) * nSector]);
				}

				//------------------------------------------------------------
				// Given rsk_sup, rK_eq(rsk_sup) and Cap(rsk_sup) compute eq.
				// NT sector prices
				//------------------------------------------------------------
				OutSum = 0.0;
				for (iSec = 1; iSec < nSector + 1; iSec++) {
					OutSum += Output[iSec-1 + iYear * nSector];
				}

				for (iSec = 0; iSec < 3; iSec++) {
					p_dem[it + iSec * MaxIt] = ( (alpha_cons[iSec+2] / (1-alpha_cons[iSec+2])) * (OutSum - Output[iSec+2 + iYear * nSector]) ) /
						( prod[iSec+2 + iYear * nSector] * pow(sk_sup[iSec+2],alpha_prod[iSec+2 + (nYear-1) * nSector]) *
					  	pow(Cap[iSec+2],1-alpha_prod[iSec+2 + (nYear-1) * nSector]) );
					p_dem[it + iSec * MaxIt] = p_dem[it + iSec * MaxIt] * (Pindex[iYear-12]/Pindex[0]);	
				}

				for (iSec = 1; iSec < nSector + 1; iSec++)
				{
					if (sk_sup[iSec-1] > 0)
					{
						rsk_dem[(iSec-1) + it * nSector] = 
							fmin(alpha_prod[iSec-1 + (nYear-1) * nSector] * Output[iSec-1 + iYear * nSector] / (sk_sup[iSec-1]), 1000.0);
					}
					else
					{
						rsk_dem[(iSec-1) + it * nSector] =
							fmin(2*rsk_sup[(iSec-1) + it * nSector], 1000.0);
					}
				}

				check11 = fabs(log(rsk_sup[0 + it * nSector]) - log(rsk_dem[0 + it * nSector]));
				check12 = fabs(log(rsk_sup[1 + it * nSector]) - log(rsk_dem[1 + it * nSector]));
				check13 = fabs(log(rsk_sup[2 + it * nSector]) - log(rsk_dem[2 + it * nSector]));
				check14 = fabs(log(rsk_sup[3 + it * nSector]) - log(rsk_dem[3 + it * nSector]));
				check15 = fabs(log(rsk_sup[4 + it * nSector]) - log(rsk_dem[4 + it * nSector]));
				
				check31  = fabs(log(p_sup[it + 0 * MaxIt]) - log(p_dem[it + 0 * MaxIt]));
				check32  = fabs(log(p_sup[it + 1 * MaxIt]) - log(p_dem[it + 1 * MaxIt]));
				check33  = fabs(log(p_sup[it + 2 * MaxIt]) - log(p_dem[it + 2 * MaxIt]));

				if (check11 <= crit1 && check12 <= crit1 && check13 <= crit1 && check14 <= crit1 && check15 <= crit1 && 
					check31 <= crit3 && check32 <= crit3 && check33 <= crit3)
				{
					converged = 1;
				}

				// Step size updating
				if (converged == 0 && it < MaxIt - 1)
				{
					if (it < 30)
					{
						for (iSec = 1; iSec < nSector + 1; iSec++)
						{
							step[iSec-1] = 1.0/6.0;
						}
						for (iSec = 0; iSec < 3; iSec++)
						{
							step3[iSec] = 1.0/6.0;
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
						}

						for (iSec = 0; iSec < 3; iSec++) {
							if (SameSign(1.0, p_dem[it + iSec * MaxIt] - p_dem[it-1 + iSec * MaxIt]) !=
								SameSign(1.0, p_dem[it-1 + iSec * MaxIt] - p_dem[it-2 + iSec * MaxIt])
								&& flag3[iSec] == 0)
							{
								flag3[iSec] = it;
								step3[iSec] = step3[iSec] / 1.3;
							}

							if (SameSign(1.0, p_dem[it + iSec * MaxIt] - p_dem[it-1 + iSec * MaxIt]) !=
								SameSign(1.0, p_dem[it-1 + iSec * MaxIt] - p_dem[it-2 + iSec * MaxIt])
								&& it >= flag3[iSec] + 5)
							{
								flag3[iSec] = it;
								step3[iSec] = step3[iSec] / 1.3;
							}
						}
					}


					// Update skill returns and NT prices
					for (iSec = 1; iSec < nSector + 1; iSec++)
					{
						rsk_sup[(iSec-1) + (it+1)*nSector] = rsk_sup[(iSec-1) + it*nSector] - 
							step[iSec-1] * (rsk_sup[(iSec-1) + it*nSector] - rsk_dem[(iSec-1) + it*nSector]);
					}
					p_sup[it+1 + 0 * MaxIt] = p_sup[it + 0 * MaxIt] - step3[0] * (p_sup[it + 0 * MaxIt] - p_dem[it + 0 * MaxIt]);
					p_sup[it+1 + 1 * MaxIt] = p_sup[it + 1 * MaxIt] - step3[1] * (p_sup[it + 1 * MaxIt] - p_dem[it + 1 * MaxIt]);
					p_sup[it+1 + 2 * MaxIt] = p_sup[it + 2 * MaxIt] - step3[2] * (p_sup[it + 2 * MaxIt] - p_dem[it + 2 * MaxIt]);
				} // end step updating
				else
				{
					// Record equilibrium quantities
					for (iSec = 0; iSec < nSector + 1; iSec++)
					{
						if (iSec > 0)
						{
							rsk_eq[(iSec-1) + iYear*nSector]    = rsk_dem[(iSec-1) + it*nSector];
							sk_sup_eq[(iSec-1) + iYear*nSector] = sk_sup[(iSec-1)];
							emp_eq[iSec + iYear*(nSector+1)]    = emp[iSec + it*(nSector+1)];
							CapitalSim[iSec-1 + iYear*nSector]  = Cap[iSec-1];
							price[2 + iYear*nSector]            = p_dem[it + 0 * MaxIt];
							price[3 + iYear*nSector]            = p_dem[it + 1 * MaxIt];
							price[4 + iYear*nSector]            = p_dem[it + 2 * MaxIt];
							Welfare[iYear-12]                   = Welfare[iYear-12] + rK_eq[iSec-1 + iYear * nSector] * CapitalSim[iSec-1 + iYear*nSector];
							iter[iYear] = it;
						}
						else
						{
							emp_eq[iSec + iYear*(nSector+1)] = emp[iSec + it*(nSector+1)];
							iter[iYear] = it;
						}
					}
					
				}

				it += 1;
			} // end convergence for returns to skill

			printf("\n");
			printf("RE Iteration: %d\n", RE_loop);
			printf("\n");
			printf("rK1: %9.3f, rK2: %9.3f, rK3: %9.3f\n", rK_eq[0 + iYear * nSector], rK_eq[1 + iYear * nSector], rK_eq[2 + iYear * nSector]);
			printf("rK4: %9.3f, rK5: %9.3f\n", rK_eq[3 + iYear * nSector], rK_eq[4 + iYear * nSector]);
			printf("\n");
			printf("rsk1: %9.3f, rsk2: %9.3f, rsk3: %9.3f\n", rsk_eq[0 + iYear*nSector], rsk_eq[1 + iYear*nSector], rsk_eq[2 + iYear*nSector]);
			printf("rsk4: %9.3f, rsk5: %9.3f\n", rsk_eq[3 + iYear*nSector], rsk_eq[4 + iYear*nSector]);
			printf("\n");
			printf("sk1: %9.3f, sk2: %9.3f, sk3: %9.3f\n", sk_sup_eq[0 + iYear*nSector], sk_sup_eq[1 + iYear*nSector], sk_sup_eq[2 + iYear*nSector]);
			printf("sk5: %9.3f, sk4: %9.3f\n", sk_sup_eq[3 + iYear*nSector], sk_sup_eq[4 + iYear*nSector]);
			printf("\n");
			printf("Unemp: %9.3f, Emp1: %9.3f, Emp2: %9.3f\n", emp_eq[0 + iYear*(nSector+1)], emp_eq[1 + iYear*(nSector+1)], emp_eq[2 + iYear*(nSector+1)]);
			printf("Emp3: %9.3f, Emp4: %9.3f, Emp5: %9.3f\n", emp_eq[3 + iYear*(nSector+1)], emp_eq[4 + iYear*(nSector+1)], emp_eq[5 + iYear*(nSector+1)]);
			printf("\n");
			printf("Fraction unemployed: %9.3f\n", emp_eq[0 + iYear*(nSector+1)] / (emp_eq[0 + iYear*(nSector+1)]+emp_eq[1 + iYear*(nSector+1)]+emp_eq[2 + iYear*(nSector+1)]+emp_eq[3 + iYear*(nSector+1)]+emp_eq[4 + iYear*(nSector+1)]+emp_eq[5 + iYear*(nSector+1)]));
			printf("\n");
			printf("Capital1: %9.3f, Capital2: %9.3f, Capital3: %9.3f\n", CapitalSim[0 + iYear*nSector], CapitalSim[1 + iYear*nSector], CapitalSim[2 + iYear*nSector]);
			printf("Capital4: %9.3f, Capital5: %9.3f\n", CapitalSim[3 + iYear*nSector], CapitalSim[4 + iYear*nSector]);
			printf("\n");
			printf("Welfare: %9.3f\n", Welfare[iYear-12]);
			printf("\n");
			printf("Benefits: %9.3f\n", Benefits[iYear-12]);
			printf("\n");
			printf("Price1: %9.3f, Price2: %9.3f, Price3: %9.3f\n", price[0 + iYear*nSector], price[1 + iYear*nSector], price[2 + iYear*nSector]);
			printf("Price4: %9.3f, Price5: %9.3f\n", price[3 + iYear*nSector], price[4 + iYear*nSector]);
			printf("\n");
			printf("Iter: %d\n", iter[iYear]);
			printf("\n");
			printf("Iter2: %d\n", iter2[iYear]);
			printf("\n");
			printf("Check1: %9.3f, %9.3f,  %9.3f,  %9.3f,  %9.3f\n", check11, check12, check13, check14, check15);
			printf("\n");
			printf("Check2: %9.3f, %9.3f,  %9.3f,  %9.3f,  %9.3f\n", check21, check22, check23, check24, check25);
			printf("\n");
			printf("Check3: %9.3f, %9.3f, %9.3f\n", check31, check32, check33);
			printf("\n");
			printf("---------------------------------------------------------------------------\n");


		} // End iYear loop

		for (iSec = 1; iSec < nSector + 1; iSec++)
		{
			for (iYear = nYear-1; iYear < (LastYr_sim-FirstYr+1) - 1; iYear++)
			{
				avec[(iSec-1) + iYear*nSector] = 
					rsk_eq[(iSec-1) + (iYear+1)*nSector] / rsk_eq[(iSec-1) + iYear*nSector];
			}
			avec[(iSec-1) + (LastYr_sim-FirstYr)*nSector] = 1.0;
		}
		
		// Write skill returns to file
		if (RE_loop == 1)
		{
			fp = fopen("Output/RSk_loop1.txt", "w");
			fprintf(fp, "Year\tRSk_1\tRSk_2\tRSk_3\tRSk_4\tRSk_5\tAvec1\tAvec2\tAvec3\tAvec4\tAvec5\n");
			for (iYear = FirstYr_sim - FirstYr; iYear < LastYr_sim - FirstYr; iYear++)//(iYear = FirstYr_sim - FirstYr; iYear < LastYr_sim - FirstYr + 1; iYear++)
			{
				fprintf(fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", iYear+1996, 
					rsk_eq[0 + iYear * nSector], rsk_eq[1 + iYear * nSector], rsk_eq[2 + iYear * nSector], rsk_eq[3 + iYear * nSector], rsk_eq[4 + iYear * nSector], 
					avec[0 + iYear*nSector], avec[1 + iYear*nSector], avec[2 + iYear*nSector], avec[3 + iYear*nSector], avec[4 + iYear*nSector]);
			}
			fclose(fp);
		}
		else if (RE_loop == 2)
		{
			fp = fopen("Output/RSk_loop2.txt", "w");
			fprintf(fp, "Year\tRSk_1\tRSk_2\tRSk_3\tRSk_4\tRSk_5\tAvec1\tAvec2\tAvec3\tAvec4\tAvec5\n");
			for (iYear = FirstYr_sim - FirstYr; iYear < LastYr_sim - FirstYr; iYear++)//(iYear = FirstYr_sim - FirstYr; iYear < LastYr_sim - FirstYr + 1; iYear++)
			{
				fprintf(fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", iYear+1996, 
					rsk_eq[0 + iYear * nSector], rsk_eq[1 + iYear * nSector], rsk_eq[2 + iYear * nSector], rsk_eq[3 + iYear * nSector], rsk_eq[4 + iYear * nSector], 
					avec[0 + iYear*nSector], avec[1 + iYear*nSector], avec[2 + iYear*nSector], avec[3 + iYear*nSector], avec[4 + iYear*nSector]);
			}
			fclose(fp);
		}
		else if (RE_loop == 3)
		{
			fp = fopen("Output/RSk_loop3.txt", "w");
			fprintf(fp, "Year\tRSk_1\tRSk_2\tRSk_3\tRSk_4\tRSk_5\tAvec1\tAvec2\tAvec3\tAvec4\tAvec5\n");
			for (iYear = FirstYr_sim - FirstYr; iYear < LastYr_sim - FirstYr; iYear++)//(iYear = FirstYr_sim - FirstYr; iYear < LastYr_sim - FirstYr + 1; iYear++)
			{
				fprintf(fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", iYear+1996, 
					rsk_eq[0 + iYear * nSector], rsk_eq[1 + iYear * nSector], rsk_eq[2 + iYear * nSector], rsk_eq[3 + iYear * nSector], rsk_eq[4 + iYear * nSector], 
					avec[0 + iYear*nSector], avec[1 + iYear*nSector], avec[2 + iYear*nSector], avec[3 + iYear*nSector], avec[4 + iYear*nSector]);
			}
			fclose(fp);
		}


		RE_loop += 1;
	}

	fp = fopen("Output/Simulation5.txt", "w");
	fprintf(fp, "Year\tPrice1\tPrice2\tPrice3\tPrice4\tPrice5\tPriceIndex\tZ1\tZ2\tZ3\tZ4\tZ5\tRSk1\tRSk2\tRSk3\tRSk4\tRSk5\tSk1\tSk2\tSk3\tSk4\tSk5\tRK1\tRK2\tRK3\tRK4\tRK5\tCapital1\tCapital2\tCapital3\tCapital4\tCapital5\tOutput1\tOutput2\tOutput3\tOutput4\tOutput5\tUnemp\tEmp1\tEmp2\tEmp3\tEmp4\tEmp5\tWelfare\tBenefits\tRWage\tIter\tIter2\n");
	for (iYear = FirstYr_sim - FirstYr; iYear < LastYr_sim - FirstYr; iYear++)//(iYear = FirstYr_sim - FirstYr; iYear < LastYr_sim - FirstYr + 1; iYear++)
	{
		fprintf(fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\n",
					 iYear+1996, price[0 + iYear*nSector], price[1 + iYear*nSector], price[2 + iYear*nSector], price[3 + iYear*nSector], price[4 + iYear*nSector], Pindex[iYear-12], 
					 z[0 + iYear * nSector], z[1 + iYear * nSector], z[2 + iYear * nSector], z[3 + iYear * nSector], z[4 + iYear * nSector], 
					 rsk_eq[0 + iYear*nSector], rsk_eq[1 + iYear*nSector], rsk_eq[2 + iYear*nSector], rsk_eq[3 + iYear*nSector], rsk_eq[4 + iYear*nSector], 
					 sk_sup_eq[0 + iYear*nSector], sk_sup_eq[1 + iYear*nSector], sk_sup_eq[2 + iYear*nSector], sk_sup_eq[3 + iYear*nSector], sk_sup_eq[4 + iYear*nSector],
					 rK_eq[0 + iYear * nSector], rK_eq[1 + iYear * nSector], rK_eq[2 + iYear * nSector], rK_eq[3 + iYear * nSector], rK_eq[4 + iYear * nSector],
					 CapitalSim[0 + iYear*nSector], CapitalSim[1 + iYear*nSector], CapitalSim[2 + iYear*nSector], CapitalSim[3 + iYear*nSector], CapitalSim[4 + iYear*nSector],
					 Output[0 + iYear * nSector], Output[1 + iYear * nSector], Output[2 + iYear * nSector], Output[3 + iYear * nSector], Output[4 + iYear * nSector],
					 emp_eq[0 + iYear*(nSector+1)], emp_eq[1 + iYear*(nSector+1)], emp_eq[2 + iYear*(nSector+1)], emp_eq[3 + iYear*(nSector+1)], emp_eq[4 + iYear*(nSector+1)], emp_eq[5 + iYear*(nSector+1)],
					 Welfare[iYear-12], Benefits[iYear-12], Real_Wages[iYear-12], iter[iYear], iter2[iYear]);
	}


	//----------------------------------------------------------
	// Free memory
	//----------------------------------------------------------
	free(EducSim);
	free(EducInit);
	free(EligSim);
	free(EligInit);
	free(FemaleSim);
	free(FemaleInit);

	free(ChoiceSim);
	free(ChoiceInit);
	free(LagSectorSim);
	free(LagSectorInit);
	free(UnempYrsSim);
	free(UnempYrsInit);
	free(ExperSim);
	free(ExperInit);
	free(LagWageSim);
	free(LagWageInit);
	free(WageSim);
	free(WelfareSim);
	free(CapitalSim);

	free(CohortWgt);
	free(z);
	free(rsk_eq);
	free(sk_sup);
	free(sk_sup_eq);
	free(price);
	free(prod);
	free(emp_eq);
	free(emp);
	free(Welfare);
	free(Benefits);
	free(Welfare_Educ);
	free(Real_Wages);
	free(Real_Wages_Educ);
	free(iter);
	free(iter2);
	free(rsk_tomorrow);

	free(rK_eq);
	free(rK_sup);
	free(rK_dem);

	free(alpha_cons);
	free(alpha_prod);

	free(Up_RSk);
	free(Low_RSk);
	free(avec);
	free(rsk_sup);
	free(rsk_dem);
	free(p_sup);
	free(p_dem);

	free(beta);
	free(sigma);
	free(kappa);
	free(xi);
	free(gamma);
	free(nu);
	free(lambda);
	free(phi);

	free(Female);
	free(Educ);
	free(Elig);
	free(Type);
	
	free(PI_RE);
	free(R_sq_RE);
	free(PI_SE);
	free(R_sq_SE);

	free(PI_Coef);
	free(w);
	free(V);
	free(Cost);

	free(Pindex);
	free(Cap);
	free(Output);

	free(step);
	free(step2);
	free(step3);
	free(flag);
	free(flag2);
	free(flag3);

	free(epsSim);
	free(etaSim);
	free(eta_aux);
	free(unif_type);

	gsl_rng_free(r);
	gsl_rng_free(s);


	return 0;
}
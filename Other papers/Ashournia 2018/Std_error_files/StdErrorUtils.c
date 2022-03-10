//------------------------------------------------------------------------------
//  StdErrorUtils.c
//
//  By Damoun Ashournia, University of Oxford
//
//------------------------------------------------------------------------------
// Functions to compute numerical derivative and standard errors.
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
    double *Low_RSk, *Up_RSk, *avec, *rsk_sup, *rsk_dem, *rsk_eq, *sk_sup, *sk_sup_eq, *step, *beta, *sigma, *gamma, *kappa, *xi, *nu, *lambda, *phi;
    double *rsk_tomorrow, *PI_Coef, *w, *V, *Cost, *check, *Data, *z;
    double crit, eps_aux, eta_aux, Emax, VMAX, Low_Wage, Up_Wage, benefit, lastwage, Emax0;

    // Parallelization
    int *Female, *Educ, *Elig, *Type, iParallel, tid, nthreads;
    
    // Auxiliary regressions
    double *auxSim;
    
    int *emp, *emp_eq, *iter, *flag, *yearDummies;
    int iAge, iPar, iPar2, iSec, iSec2, iYear, iCoh, iN, iCnt, iReg, iCount, iSize, iLagSec, iAv, RE_loop, it, iRsk, converged, *choice;
    int age, female, educ, elig, lagSec, typ, *type;
    double exper, ExperTomorrow, term2, prob;

    double *XWage, *XOther, *Y, *coef, SST, SSE, *choiceWgt;

    // Loss function
    double Loss;

    // Probabilities for drawing unemployment shocks.
    double Punemp;

    // File pointers
    FILE *fpEmaxSE, *fpRE1, *fpRE2, *fpData;


    // Running time
    clock_t cstart = 0;
    clock_t cend = 0;

    // Allocate memory
    PI_RE        = calloc(nReg*(nYear-1)*(nAge-1)*(nSector+1)*CatGender*CatEduc*CatElig*CatTypes,sizeof(double));
    R_sq_RE      = calloc((nYear-1)*(nAge-1)*(nSector+1)*CatGender*CatEduc*CatElig*CatTypes,sizeof(double));
    PI_SE        = calloc(nReg*(nAge-1)*(nSector+1)*CatGender*CatEduc*CatElig*CatTypes,sizeof(double));
    R_sq_SE      = calloc((nAge-1)*(nSector+1)*CatGender*CatEduc*CatElig*CatTypes,sizeof(double));
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
    xi           = calloc(2*nSector,sizeof(double));
    sigma        = calloc(nSector,sizeof(double));
    gamma        = calloc(nSector+1,sizeof(double));
    nu           = calloc(1,sizeof(double));
    lambda       = calloc(nSector,sizeof(double));
    phi          = calloc(4,sizeof(double));
    rsk_tomorrow = calloc(nSector,sizeof(double));
    PI_Coef      = calloc(nReg,sizeof(double));
    w            = calloc(nSector,sizeof(double));
    V            = calloc(nSector,sizeof(double));
    Cost         = calloc(nSector,sizeof(double));
    check        = calloc(nSector*nYear,sizeof(double));
    z            = calloc(nSector*nYear,sizeof(double));

    emp         = calloc((nSector+1)*MaxIt,sizeof(int));
    emp_eq      = calloc((nSector+1)*nYear,sizeof(int));
    iter        = calloc(nYear,sizeof(int));
    flag        = calloc(nSector,sizeof(int));
    yearDummies = calloc(nPop*nAge*nYear*nYear,sizeof(int));
    choice      = calloc(nSector+1,sizeof(int));
    choiceWgt   = calloc(nSector+1,sizeof(double));

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
    #define iPopCohYearNext     iN + iCoh * nPop + (iYear+1) * nCoh * nPop
    #define size                nPop * nAge * nYear
    #define iAvec2              (iSec2-1) + iYear * nSector
    #define iPISE               iReg + (age-30)*nReg + lagSec*(nAge-1)*nReg + female*(nSector+1)*(nAge-1)*nReg + educ*CatGender*(nSector+1)*(nAge-1)*nReg + elig*CatEduc*CatGender*(nSector+1)*(nAge-1)*nReg + (typ-1)*CatElig*CatEduc*CatGender*(nSector+1)*(nAge-1)*nReg
    #define iPIRE               iReg + iYear*nReg + (age-30)*(nYear-1)*nReg + lagSec*(nAge-1)*(nYear-1)*nReg + female*(nSector+1)*(nAge-1)*(nYear-1)*nReg + educ*CatGender*(nSector+1)*(nAge-1)*(nYear-1)*nReg + elig*CatEduc*CatGender*(nSector+1)*(nAge-1)*(nYear-1)*nReg + (typ-1)*CatElig*CatEduc*CatGender*(nSector+1)*(nAge-1)*(nYear-1)*nReg
    #define iCohort             iCoh + iYear * nCoh
    
    
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

    //--------------------------------------------------------------------
    // Solve model with rational expectations
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
        for (iSec = 1; iSec < nSector + 1; iSec++)
        {
            rsk_sup[(iSec-1) + 0 * nSector] = rsk_global[iSec-1];
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

                for (iSec = 1; iSec < nSector + 1; iSec++)
                {
                    rsk_sup[(iSec-1) + 0 * nSector] = fmin(fmax(rsk_eq[(iSec-1) + (iYear-1)*nSector],Low_RSk[(iSec-1)]),Up_RSk[(iSec-1)]);
                }
            }

            converged = 0;
            
            while (converged==0 && it < MaxIt)
            {
                for (iSec = 0; iSec < nSector + 1; iSec++)
                {
                    if (iSec>0)
                    {
                        sk_sup[iSec-1] = 0.0;   
                    }
                    emp[iSec + it * (nSector+1)] = 0;
                }
                
                //--------------------------------------------------------
                // Compute aggregate skill supply
                //--------------------------------------------------------
                for (iCoh = 0; iCoh < nCoh; iCoh++)
                {
                    for (iN = 0; iN < nPop; iN++)
                    {
                        for (iSec = 1; iSec < nSector + 1; iSec++)
                        {
                            age      = (iYear + 1996) - (iCoh + 1931);
                            female   = FemaleData[iPopCoh];
                            educ     = EducData[iPopCoh];
                            lagSec   = LagSectorData[iPopCohYear];
                            exper    = ExperData[iPopCohYear];
                            lastwage = LagWageData[iPopCohYear];
                            eps_aux  = sigma[iSec-1] * epsSim[iEps];
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
                                    for (iSec2 = 1; iSec2 < nSector + 1; iSec2++)
                                    {
                                        rsk_tomorrow[iSec2-1] = avec[iAvec2] * rsk_sup[(iSec2-1) + it * nSector];
                                    }
                                
                                    if (lagSec == iSec)
                                    {
                                        ExperTomorrow = exper + 1.0;
                                    }
                                    else
                                    {
                                        ExperTomorrow = gamma[iSec]*exper;
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
                    
                                V[iSec-1] = (1 - delta[iDelta])*(w[iSec-1] - Cost[iSec-1] + rho*Emax) + delta[iDelta]*(benefit + rho*Emax0) + eta_aux;
                            }
                        } // end iSec


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
                                        emp[iSec + it * (nSector+1)] = emp[iSec + it * (nSector+1)] + 1;
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
                                emp[0 + it * (nSector+1)] = emp[0 + it * (nSector+1)] + 1;
                            }
                                


                            if (age < LastAge && iYear < nYear - 1)
                            {
                                if (ChoiceData[iPopCohYear] == lagSec && lagSec > 0)
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

                                LagSectorData[iPopCohYearNext] = ChoiceData[iPopCohYear];

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
                            emp_eq[iSec + iYear*(nSector+1)]    = emp[iSec + it*(nSector+1)];
                            iter[iYear] = it;
                        }
                        else
                        {
                            emp_eq[iSec + iYear*(nSector+1)]    = emp[iSec + it*(nSector+1)];
                            iter[iYear] = it;   
                        }
                        
                    }
                }
                
                it += 1;

            } // end convergence of skill return for iYear

            //printf("%d %f %f\n", iYear, rsk_eq[0 + iYear*nSector], rsk_eq[1 + iYear*nSector]);
            //printf("%d\n", emp_eq[0 + iYear*(nSector+1)]);
            
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
                    Data[iCount + 9 * size]  = LagSectorData[iPopCohYear];
                    Data[iCount + 10 * size] = log( WageData[iPopCohYear] );
                    Data[iCount + 11 * size] = sqrt( CohortWgt[iCoh + iYear * nCoh] );

                    iCount += 1;
                } // end if

            }
        }
    }


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
    for (iLagSec = 1; iLagSec < nSector + 1; iLagSec++)
    {
        for (iSec = 1; iSec < nSector + 1; iSec++)
        {
            // Dep. var.: transition dummy
            for (iSize = 0; iSize < size; iSize++)
            {
                Y[iSize] = (Data[iSize + 9 * size]==iLagSec && Data[iSize + 8 * size]==iSec) * Data[iSize + 11 * size];
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
    // Write data for numerical derivatives
    //----------------------------------------------------------------
    fprintf(fpSEParams, "%d\t%12.10f\t", PAR_NUM, PAR_VAL);
    for (iPar = 0; iPar < nParamAux; iPar++) {
        fprintf(fpSEParams, "%12.10f\t", auxSim[iPar]);
    }
    fprintf(fpSEParams, "%lf\n", Loss);

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


//--------------------------------------------------------------
// Function to compute numerical derivatives of
// objective function at some parameter vector.
//--------------------------------------------------------------
int NumDerivatives(double *param_in) {
    
    //------------------------------------------------------
    // Declarations
    //------------------------------------------------------
    
    // Counters
    int iParam, iParam2, iGrid, iParAux;
    
    // File pointer
    FILE *fp;
    
    // Gradient vector
    double *G0;
    
    // Arrays
    double *param_low, *param_high, *parameters, *Y, *XR, *X, *auxDelta, *coef;
    double SST, SSE;
    
    // Grid points to interpolate over for each parameter
    int nGrid;
    
    
    //------------------------------------------------------
    // Allocate memory
    //------------------------------------------------------
    param_low  = calloc(nParEst, sizeof(double));
    param_high = calloc(nParEst, sizeof(double));
    parameters = calloc(nParEst, sizeof(double));
    
    // Grid points to interpolate over for each parameter
    nGrid = 20;
    
    Y        = calloc(nGrid, sizeof(double));
    XR       = calloc(nGrid*3, sizeof(double));
    X        = calloc(nParEst*nGrid, sizeof(double));
    auxDelta = calloc(nParEst*nGrid*nParamAux, sizeof(double));
    
    coef = calloc(3, sizeof(double));
    
    G0 = calloc(nParamAux*nParEst, nParam);
    
    //------------------------------------------------------
    // Read range over which interpolation is performed
    //------------------------------------------------------
    fp = fopen("Data/StdErrors_Range.txt", "r");
    
    if (fp == NULL)
    {
        perror("Failed to open file \"StdErrors_Range.txt\"");
        exit(EXIT_FAILURE);
    }
    
    for (iParam = 0; iParam < nParEst; iParam++)
    {
        fscanf(fp, "%*d %lf %lf", &param_low[iParam], &param_high[iParam]);
    }
    
    fclose(fp);

    
    // Check if ranges are well-defined
    for (iParam = 0; iParam < nParEst; iParam++)
    {
        if (param_low[iParam] > param_in[iParam] || param_high[iParam] < param_in[iParam]) {
            printf("Error: Parameter range not well-defined for parameter %d.\n", iParam+1);
            exit(EXIT_FAILURE);
        }
    }
    
    //------------------------------------------------------
    // Compute loss function at every grid-point for every
    // parameter, leaving other parameters to input values
    //------------------------------------------------------
    fpSEParams = fopen("Output/StdErrors/params.txt", "w");

    for (iParam = 0; iParam < nParEst; iParam++) {
        PAR_NUM = iParam + 1;
        for (iGrid = 0; iGrid < nGrid; iGrid++) {
            for (iParam2 = 0; iParam2 < nParEst; iParam2++) {
                if (iParam == iParam2) {
                    parameters[iParam2] = param_low[iParam] + ( (param_high[iParam] - param_low[iParam]) / nGrid ) * iGrid;
                    PAR_VAL = parameters[iParam2];
                }
                else {
                    parameters[iParam2] = param_in[iParam2];
                }
            }
            // Compute loss function
            printf("Par %2d of %d. Grid %2d of %d. L = %9.6f\n", iParam+1, nParEst, iGrid+1, nGrid, LossFcn(parameters));
        }
    }
    
    fclose(fpSEParams);
    
    
    //------------------------------------------------------
    // Compute numerical derivative
    //------------------------------------------------------
    fpSEParams = fopen("Output/StdErrors/params.txt", "r");
    
    if (fpSEParams == NULL)
    {
        perror("Failed to open file \"StdErros/params.txt\"");
        exit(EXIT_FAILURE);
    }
    
    for (iParam = 0; iParam < nParEst; iParam++) {
        for (iGrid = 0; iGrid < nGrid; iGrid++) {
            fscanf(fpSEParams, "%*d %lf", &X[iParam + iGrid * nParEst]);
            for (iParAux = 0; iParAux < nParamAux; iParAux++) {
                fscanf(fpSEParams, "%lf", &auxDelta[iParam + iGrid * nParEst + iParAux * nGrid * nParEst]);
            }
            fscanf(fpSEParams, "%*f");
        }
    }
    
    fclose(fpSEParams);
    
    
    fp = fopen("Output/StdErrors/G0.txt", "w");
    
    for (iParam = 0; iParam < nParEst; iParam++) {
        for (iParAux = 0; iParAux < nParamAux; iParAux++) {
            for (iGrid = 0; iGrid < nGrid; iGrid++) {
                Y[iGrid]              = auxDelta[iParam + iGrid * nParEst + iParAux * nGrid * nParEst];
                XR[iGrid + 0 * nGrid] = 1.0;
                XR[iGrid + 1 * nGrid] = X[iParam + iGrid * nParEst];
                XR[iGrid + 2 * nGrid] = X[iParam + iGrid * nParEst] * X[iParam + iGrid * nParEst];
            }
            fcn_linreg(Y, XR, nGrid, 3, coef, &SST, &SSE);
            G0[iParAux + iParam * nParamAux] = coef[2] + coef[3] * param_in[iParam];
            fprintf(fp, "%d\t%d\t%22.20f\n", iParam+1, iParAux+1, G0[iParAux + iParam * nParamAux]);
        }
    }
    
    fclose(fp);
    
    
    //------------------------------------------------------
    // Free memory
    //------------------------------------------------------
    free(param_low);
    free(param_high);
    free(parameters);
    free(Y);
    free(XR);
    free(X);
    free(auxDelta);
    free(coef);
    free(G0);
    
    
    return 0;
}

//--------------------------------------------------------------
// Function that uses numerical derivatives to
// compute standard errors.
//--------------------------------------------------------------
int ComputeStdErrors() {
    //------------------------------------------------------
    // Declarations
    //------------------------------------------------------
    
    // File pointer
    FILE *fp;
    
    // Counters
    int iPar, iAux, iPar2, S;
    
    // Arrays
    double *auxDelta, *auxDeltaTranspose, *GTmp, *G, *Cov;
    
    
    //------------------------------------------------------
    // Allocate memory
    //------------------------------------------------------
    auxDelta          = calloc(nParamAux*nParEst, sizeof(double));
    auxDeltaTranspose = calloc(nParamAux*nParEst, sizeof(double));
    
    GTmp = calloc(nParEst*nParamAux, sizeof(double));
    G    = calloc(nParEst*nParEst, sizeof(double));
    Cov  = calloc(nParEst*nParEst, sizeof(double));
    
    //------------------------------------------------------
    // Read data
    //------------------------------------------------------
    fp = fopen("Output/StdErrors/G0.txt", "r");
    
    if (fp == NULL)
    {
        perror("Failed to open file \"StdErros/G0.txt\"");
        exit(EXIT_FAILURE);
    }
    
    for (iPar = 0; iPar < nParEst; iPar++) {
        for (iAux = 0; iAux < nParamAux; iAux++) {
            fscanf(fp, "%*d %*d %lf", &auxDelta[iAux + iPar * nParamAux]);
            auxDeltaTranspose[iPar + iAux * nParEst] = auxDelta[iAux + iPar * nParamAux];
        }
    }
    
    fclose(fp);
    
    
    //------------------------------------------------------
    // Compute covariance matrix
    //------------------------------------------------------
    fcn_matmul(auxDeltaTranspose, invCov, nParEst, nParamAux, nParamAux, GTmp);
    fcn_matmul(GTmp, auxDelta, nParEst, nParamAux, nParEst, G);
    
    fcn_invX(G, nParEst, Cov);
    
    S = 1; // number of simulations
    
    fp = fopen("Output/StdErrors/VarCov.txt", "w");
    
    for (iPar = 0; iPar < nParEst; iPar++) {
        for (iPar2 = 0; iPar2 < nParEst; iPar2++) {
            Cov[iPar + iPar2 * nParEst] = ( (S + 1) / S ) * Cov[iPar + iPar2 * nParEst];
            if (iPar2 == nParEst - 1) {
                fprintf(fp, "%25.23f\n", Cov[iPar + iPar2 * nParEst]);
            }
            else {
                fprintf(fp, "%25.23f\t", Cov[iPar + iPar2 * nParEst]);
            }
        }
    }
    
    fclose(fp);
    
    fp = fopen("Output/StdErrors/StdErrors.txt", "w");
    
    for (iPar = 0; iPar < nParEst; iPar++) {
        fprintf(fp, "%25.23f\n", sqrt( Cov[iPar + iPar * nParEst] ) );
    }
    
    fclose(fp);
    
    
    //------------------------------------------------------
    // Free memory
    //------------------------------------------------------
    free(auxDelta);
    free(auxDeltaTranspose);
    free(GTmp);
    free(G);
    free(Cov);
    
    return 0;
}
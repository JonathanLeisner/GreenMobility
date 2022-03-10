//------------------------------------------------------------------------------
// Estimate.c
//
// By Damoun Ashournia, University of Oxford
//
//------------------------------------------------------------------------------
// Main estimation file. Minimizes the objective function using the Nelder-Mead
// simplex method.
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Use the snippet below to allow gcc4.9 to work with the Accelerate framework
//------------------------------------------------------------------------------
// #ifndef __has_extension
// #define __has_extension(x) 0
// #endif
//
// #define vImage_Utilities_h
// #define vImage_CVUtilities_h
//------------------------------------------------------------------------------

#include <stdlib.h>                    // Standard C library
#include <stdio.h>                     // Standard C input output library
#include <math.h>                      // C math library
#include <time.h>                      // In order to compute running time
#include <omp.h>                       // OpenMP
#include <nlopt.h>                     // NLopt for numerical optimization

// GSL
#include <gsl/gsl_rng.h>               // Random number generator
#include <gsl/gsl_randist.h>           // distributions to be drawn from

// Mac OS X
#include <Accelerate/Accelerate.h>     // XCode Accelerate framework - holds lapack and blas

//------------------------------------------------------------------------------
// Declare globals
//------------------------------------------------------------------------------

// Counter for objective function evaluations
int count;

// File pointers
FILE *fpParOut;

// Initial conditions
int *CohortBirthData, *FirstYearData, *AgeData, *EducData, *FemaleData, *UIFundData;
int *ChoiceData, *LagChoiceData, *UnempYearsData, *EligData;
double *WageData, *LagWageData, *CostSim, *VSim, *ExperData;

// Starting values
double *param_start;
int *mask;

// Cohort weights
double *CohortWgt, *CohortSizeData;;

// Aggregate data
double *OutputData, *IncomeShares, *CapitalData, *rsk_global, *rsk_year;

// Auxiliary regressions
double *auxData, *covData, *invCov;;

// GSL Random Number Generator declarations
double *epsSim, *etaSim, *unif_type;

// Spa parameters
int FirstYr = 1996;
int LastYr  = 2008;
int nYear   = 2008 - 1996 + 1;

int FirstAge = 30;
int LastAge  = 65;
int nAge     = 65 - 30 + 1;

int FirstCoh = 1996 - 65;
int LastCoh  = 2008 - 30;
int nCoh     = (2008 - 30) - (1996 - 65) + 1;

int nPop = 2000; // # of simulated people
int Intp = 1500; // # of people used in Emax regressions
int nReg = 36;   // # of regressors in Emax regressions

int nParam = 65; // # of parameters to be estimated in model
int nParamBeta = 4; // # of paramters in wage offer function
int nChoices    = 6; // # of choices
int nSector     = 5; // # of sectors
int nFactors    = 2; // # of factors in production function

int nDraws = 300; // # of points used for Monte Carlo integration

int nCohSizes = 2808; // # of observations in cohost size data

int CatGender = 2; // # of gender categories
int CatEduc   = 2; // # of education categories
int CatElig   = 2; // # of eligibility categories
int CatTypes  = 2; // # of types

int MaxIt = 100; // Max # of iteration for skill-price convergence algorithm

int nRegWage  = 18; // # of regressors in auxiliary wage regressions
int nRegOther = 19; // # of regressors in auxiliary employment regressions
int nParamAux = 665; // # of parameters in auxiliary regressions

double rho = 0.95; // Discount factor

double compRate = 0.9; // Degree of compensation in UI system
double UIbar    = 203090 * 0.846575342 / 1702; // Max UI benefit in 2008 times CPI(2000=1) divided by work-hours per year
double WA       = 13732 * 12 * 0.776152482 / 1702; // Welfare assistance (kontanthjaelp) in 2012 times months times CPI(2000=1) divided by work-hours per year


//------------------------------------------------------------------------------
// Definitions
//------------------------------------------------------------------------------
#define NM_MAXITER     100000
#define NM_TOLER       1e-4

#define iPopCohYear    iN + iCoh * nPop + iYear * nCoh * nPop
#define iPopCoh        iN + iCoh * nPop

#define iEps           iN + iCoh * nPop + iYear * nCoh * nPop + iCh * nYear * nCoh * nPop
#define iEta           iN + iCoh * nPop + iYear * nCoh * nPop + iCh * nYear * nCoh * nPop
#define iTyp           iN + iCoh * nPop


//------------------------------------------------------------------------------
// Load modules
//------------------------------------------------------------------------------
#include "CUtils.c"     // C utilities
// #include "Emax.c"       // Model solution
// #include "ObjFcn.c"     // Objective Function
#include "Setup.c"


//------------------------------------------------------------------------------
// Main function
//------------------------------------------------------------------------------
int main()
{
   //-----------------------------------------------
   // Initialize
   //-----------------------------------------------

   // Counters
   int iParam, iN, iCoh, iYear, iCh, iCnt, k;

   // Starting values
   double *param_in, *param_out, *param, *beta, *sigma, *gamma, *delta, *kappa, *xi, *nu, *lambda, *phi;

   // GSL Random Number Generator declarations
   const gsl_rng_type *T;
   gsl_rng *r;

   //-----------------------------------------------
   // Allocate memory
   //-----------------------------------------------

   // GSL Random Number Generator
   epsSim    = calloc(nPop*nCoh*nYear*nChoices,sizeof(double));
   etaSim    = calloc(nPop*nCoh*nYear*nChoices,sizeof(double));
   unif_type = calloc(nPop*nCoh,sizeof(double));

   // Starting values
   param_start = calloc(nParam,sizeof(double));
   param_in    = calloc(nParam,sizeof(double));
   param_out   = calloc(nParam,sizeof(double));
   param       = calloc(nParam,sizeof(double));
   mask        = calloc(nParam,sizeof(int));

   // Paramters
   beta        = calloc(nParamBeta*nSector,sizeof(double));
   sigma       = calloc(nChoices,sizeof(double));
   gamma       = calloc(nChoices,sizeof(double));
   delta       = calloc(nChoices,sizeof(double));
   kappa       = calloc(4,sizeof(double));
   xi          = calloc(2*nChoices,sizeof(double));
   nu          = calloc(1,sizeof(double));
   lambda      = calloc(nChoices,sizeof(double));
   phi         = calloc(4,sizeof(double));

   // Initial conditions
   CohortBirthData = calloc(nPop*nCoh*nYear,sizeof(int));
   ExperData       = calloc(nPop*nCoh*nYear,sizeof(double));
   ChoiceData      = calloc(nPop*nCoh*nYear,sizeof(int));
   LagChoiceData   = calloc(nPop*nCoh*nYear,sizeof(int));
   UnempYearsData  = calloc(nPop*nCoh*nYear,sizeof(int));
   EligData        = calloc(nPop*nCoh*nYear,sizeof(int));
   WageData        = calloc(nPop*nCoh*nYear,sizeof(double));
   LagWageData     = calloc(nPop*nCoh*nYear,sizeof(double));
   CostSim         = calloc(nPop*nCoh*nYear,sizeof(double));
   VSim            = calloc(nPop*nCoh*nYear,sizeof(double));
   FirstYearData   = calloc(nPop*nCoh,sizeof(int));
   AgeData         = calloc(nPop*nCoh,sizeof(int));
   EducData        = calloc(nPop*nCoh,sizeof(int));
   FemaleData      = calloc(nPop*nCoh,sizeof(int));
   UIFundData      = calloc(nPop*nCoh,sizeof(int));

   // Cohort sizes
   CohortSizeData = calloc(nCoh*nYear,sizeof(double));
   CohortWgt      = calloc(nCoh*nYear,sizeof(double));

   // Aggregate data
   OutputData   = calloc(nYear*nSector,sizeof(double));
   IncomeShares = calloc(nYear*nSector*nFactors,sizeof(double));
   CapitalData  = calloc(nYear*nSector,sizeof(double));
   rsk_global   = calloc(nSector,sizeof(double));
   rsk_year     = calloc(nYear*nSector,sizeof(double));

   // Auxiliary regressions
   auxData = calloc(nParamAux,sizeof(double));
   covData = calloc(nParamAux*nParamAux,sizeof(double));
   invCov  = calloc(nParamAux*nParamAux,sizeof(double));


   //-----------------------------------------------
   // Prepare for estimation
   //-----------------------------------------------
   Setup();

   // Unpack paramters
   unpackPars(param_start, beta, sigma, gamma, delta, xi, kappa, nu, lambda, phi);


   //-----------------------------------------------
   // Draw random variables for simulation in
   // ObjFcn.c once and for all.
   //-----------------------------------------------

   // Epsilon shock (normal)
   T = gsl_rng_mt19937;       // Algorithm
   r = gsl_rng_alloc(T);      // Allocate memory
   gsl_rng_set(r,5544332211); // Set seed

   for (iCh = 0; iCh < nChoices; iCh++)
   {
      for (iYear = 0; iYear < nYear; iYear++)
      {
         for (iCoh = 0; iCoh < nCoh; iCoh++)
         {
            for (iN = 0; iN < nPop; iN++)
            {
               epsSim[iEps] = gsl_ran_gaussian(r, 1);
            }
         }
      }
   }

   // Preference shock (EV - see eta_aux in LossFcn.c)
   gsl_rng_set(r,223344556); // Set seed
   for (iCh = 0; iCh < nChoices; iCh++)
   {
      for (iYear = 0; iYear < nYear; iYear++)
      {
         for (iCoh = 0; iCoh < nCoh; iCoh++)
         {
            for (iN = 0; iN < nPop; iN++)
            {
               etaSim[iEta] = gsl_rng_uniform_pos(r);
            }
         }
      }
   }

   // Type probability
   gsl_rng_set(r,443377118); // Set seed
   for (iCoh = 0; iCoh < nCoh; iCoh++)
   {
      for (iN = 0; iN < nPop; iN++)
      {
         unif_type[iTyp] = gsl_rng_uniform_pos(r);
      }
   }

   gsl_rng_free (r);


   //-----------------------------------------------
   // Set mask = 0 for parameters not to be
   // optimized over, then assign input parameters.
   //-----------------------------------------------
   iCnt = 0;

   // Beta
   for (iParam = 0; iParam < nParamBeta*nSector; iParam++)
   {
      mask[iCnt] = 0;
      iCnt += 1;
   }

   // Sigma
   for (iParam = 0; iParam < nChoices; iParam++)
   {
      mask[iCnt] = 0;
      iCnt += 1;
   }

   // Gamma
   for (iParam = 0; iParam < nChoices; iParam++)
   {
       mask[iCnt] = 0;
       iCnt += 1;
   }

   // Delta
   for (iParam = 0; iParam < nChoices; iParam++)
   {
       mask[iCnt] = 0;
       iCnt += 1;
   }

   // Xi
   for (iParam = 0; iParam < 2*nChoices; iParam++)
   {
       mask[iCnt] = 0;
       iCnt += 1;
   }

   // Kappa
   for (iParam = 0; iParam < 4; iParam++)
   {
       mask[iCnt] = 0;
       iCnt += 1;
   }

   // Nu
   mask[iCnt] = 1;
   iCnt += 1;

   // Lambda
   for (iParam = 0; iParam < nChoices; iParam++)
   {
       mask[iCnt] = 0;
       iCnt += 1;
   }

   // Phi
   for (iParam = 0; iParam < 4; iParam++)
   {
       mask[iCnt] = 0;
       iCnt += 1;
   }

   // Load param_in vector
   k = 0;
   for (iParam = 0; iParam < nParam; iParam++)
   {
      if (mask[iParam]==1)
      {
         param_in[k] = param_start[iParam];
         k += 1;
      }
   }

   // Number of parameters to estimate
   printf("Estimating %d parameter(s)\n", k);
   // OBS! k = number of paramters to be optimized over in Nelder-Mead algorithm

   //-----------------------------------------------
   // Call Nelder-Mead minimizer
   //-----------------------------------------------
   double *in_pars;
   in_pars = calloc(k,sizeof(double));

   for (iParam = 0; iParam < k; iParam++)
   {
      in_pars[iParam] = param_in[iParam];
   }
/*
   nlopt_opt opt;
   opt = nlopt_create(NLOPT_LN_NELDERMEAD, k);
   nlopt_set_min_objective(opt, ObjFcn, NULL);
   nlopt_set_xtol_rel(opt, 1e-4);

   double minf; // the minimum objective value, upon return

   if (nlopt_optimize(opt, in_pars, &minf) < 0) {
       printf("nlopt failed!\n");
   }
   else {
      printf("found minimum after %d evaluations\n", count);
      printf("found minimum at f() = %15.12f\n", minf);
   }

   nlopt_destroy(opt);
*/

   // Unmask parameters
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

   // Write output parameters to file
   fpParOut = fopen("Output/EstimationResults.txt", "w");

   for (iParam = 0; iParam < nParam; iParam++)
   {
      fprintf(fpParOut, "%lf\n", param_out[iParam]);
   }

   fclose(fpParOut);


   //-----------------------------------------------
   // Free allocated memory
   //-----------------------------------------------

   // GSL Random Number Generator
   free(epsSim);
   free(etaSim);
   free(unif_type);

   // Starting values
   free(param_start);
   free(param_in);
   free(param_out);
   free(param);
   free(mask);
   free(beta);
   free(sigma);
   free(gamma);
   free(delta);
   free(kappa);
   free(xi);
   free(nu);
   free(lambda);
   free(phi);

   free(in_pars);

   // Initial conditions
   free(CohortBirthData);
   free(ExperData);
   free(ChoiceData);
   free(LagChoiceData);
   free(UnempYearsData);
   free(EligData);
   free(WageData);
   free(LagWageData);
   free(CostSim);
   free(VSim);
   free(FirstYearData);
   free(AgeData);
   free(EducData);
   free(FemaleData);
   free(UIFundData);

   // Cohort sizes
   free(CohortSizeData);
   free(CohortWgt);

   // Aggregate data
   free(OutputData);
   free(IncomeShares);
   free(CapitalData);
   free(rsk_global);
   free(rsk_year);

   // Auxiliary regressions
   free(auxData);
   free(covData);
   free(invCov);

   return 0;
}

//------------------------------------------------------------------------------
// Simulate.c
//
// By Damoun Ashournia, University of Oxford
//
//------------------------------------------------------------------------------
// Main simulation file. Calls the simulation functions:
//   Simulation0():  No shocks.
//   Simulation1():  Manufacturing price shock
//   Simulation2():  Unemployment shock to manufacturing workers.
//   Simulation3():  Both the above shocks.  
//   Simulation4():  Both the shocks, no utility cost of switching.
//   Simulation5():  Both the shocks, half UI period.
//   Simulation6():  Both the shocks, 50% UI coverage.  
//   Simulation7():  Both shocks, lower UIbar.
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Use the snippet below to allow gcc4.9 to work with the Accelerate framework
//------------------------------------------------------------------------------
#ifndef __has_extension
#define __has_extension(x) 0
#endif

#define vImage_Utilities_h
#define vImage_CVUtilities_h
//------------------------------------------------------------------------------

#include <stdlib.h>                    // Standard C library
#include <stdio.h>                     // Standard C input output library
#include <math.h>		                  // C math library
#include <time.h>                      // In order to compute running time
#include <omp.h>                       // OpenMP

// GSL
#include <gsl/gsl_multimin.h>          // GSL multidimensional minimizer library
#include <gsl/gsl_rng.h>				   // Random number generator
#include <gsl/gsl_randist.h>			   // distributions to be drawn from

// Mac OS X
#include <Accelerate/Accelerate.h>     // XCode Accelerate framework - holds lapack and blas


//------------------------------------------------------------------------------
// Declare globals
//------------------------------------------------------------------------------

// Initial conditions
int *CohortBirthData, *FirstYearData, *AgeData, *EducData, *FemaleData, *UIFundData;
int *ChoiceData, *LagSectorData, *UnempYearsData;
double *WageData, *LagWageData, *ExperData;

// Cohort weights
double *CohortSizeData;

// Aggregate data
double *OutputData, *IncomeShares, *CapitalData, *rKData;

// Parameters
double *param, *delta, *Delta;


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

int nParam = 54; // # of parameters to be estimated in model
int nParamBeta = 4; // # of paramters in wage offer function
int nSector     = 5; // # of sectors
int nFactors    = 2; // # of factors in production function

int nDraws = 300; // # of points used for Monte Carlo integration

int nCohSizes = 2808; // # of observations in cohost size data

int CatGender = 2; // # of gender categories
int CatEduc   = 2; // # of education categories
int CatElig   = 2; // # of eligibility categories
int CatTypes  = 2; // # of types

int MaxIt = 80; // Max # of iteration for skill-price convergence algorithm

int nRegWage = 18; // # of regressors in auxiliary wage regressions
int nRegEmp  = 19; // # of regressors in auxiliary employment regressions
int nRegTr   = 19; // # of regressors in auxiliary transition regressions

double rho = 0.95; // Discount factor

double compRate = 0.9; // Degree of compensation in UI system
double UIbar    = 203090 * 0.846575342 / 1702; // Max UI benefit in 2008 times CPI(2000=1) divided by work-hours per year
double WA       = 13732 * 12 * 0.776152482 / 1702; // Welfare assistance (kontanthjaelp) in 2012 times months times CPI(2000=1) divided by work-hours per year
int UImax = 4;

// Variables used in counterfactual simulation
int nPop_sim     = 2000;
int sim_t        = 300;
int FirstYr_sim  = 2008;
int LastYr_sim   = 2008 + 300;
int FirstCoh_sim = 2008 - 65;
int LastCoh_sim  = 2008 + 300 - 30;
int HalfYr_sim   = 2158;
int Yr_shock     = 2158 + 1;

int nCoh_sim  = (2008 + 300 - 30) - (2008 - 65) + 1;
int nCoh_init = (2008 - 30) - (2008 - 65) + 1;



//------------------------------------------------------------------------------
// Load modules
//------------------------------------------------------------------------------
#include "CUtils.c"     // C utilities
#include "Emax.c"	      // Model solution
#include "Setup.c"      // Load data and parameter values
#include "Simulation0.c"
#include "Simulation1.c"
#include "Simulation2.c"
#include "Simulation3.c"
#include "Simulation4.c"
#include "Simulation5.c"
#include "Simulation6.c"
#include "Simulation7.c"


//------------------------------------------------------------------------------
// Main function
//------------------------------------------------------------------------------
int main()
{
   //-----------------------------------------------
   // Initialize
   //-----------------------------------------------
   
   // Parameters
   double *beta, *sigma, *gamma, *kappa, *xi, *nu, *lambda, *phi;


   //-----------------------------------------------
   // Allocate memory
   //-----------------------------------------------

   // Paramters
   param = calloc(nParam,sizeof(double));
   beta  = calloc(nParamBeta*nSector,sizeof(double));
   sigma = calloc(nSector,sizeof(double));
   gamma = calloc(nSector+1,sizeof(double));
   kappa = calloc(4,sizeof(double));
   xi    = calloc(2*nSector,sizeof(double));
   nu    = calloc(1,sizeof(double));
   lambda      = calloc(nSector,sizeof(double));
   phi         = calloc(4,sizeof(double));

   // Unemployment transitions
   delta = calloc(CatGender*nAge*CatEduc*(nSector+1),sizeof(double));
   Delta = calloc((LastYr_sim-FirstYr+1)*CatGender*nAge*CatEduc*(nSector+1),sizeof(double)); // for shock in Simulation2.c

   // Initial conditions
   CohortBirthData = calloc(nPop*nCoh*nYear,sizeof(int));
   ExperData       = calloc(nPop*nCoh*nYear,sizeof(double));
   ChoiceData      = calloc(nPop*nCoh*nYear,sizeof(int));
   LagSectorData   = calloc(nPop*nCoh*nYear,sizeof(int));
   UnempYearsData  = calloc(nPop*nCoh*nYear,sizeof(int));
   WageData        = calloc(nPop*nCoh*nYear,sizeof(double));
   LagWageData     = calloc(nPop*nCoh*nYear,sizeof(double));
   FirstYearData   = calloc(nPop*nCoh,sizeof(int));
   AgeData         = calloc(nPop*nCoh,sizeof(int));
   EducData        = calloc(nPop*nCoh,sizeof(int));
   FemaleData      = calloc(nPop*nCoh,sizeof(int));
   UIFundData      = calloc(nPop*nCoh,sizeof(int));

   // Cohort sizes
   CohortSizeData = calloc(nCoh*nYear,sizeof(double));

   // Aggregate data
   OutputData   = calloc(nYear*nSector,sizeof(double));
   IncomeShares = calloc(nYear*nSector*nFactors,sizeof(double));
   CapitalData  = calloc(nYear*nSector,sizeof(double));
   rKData       = calloc(nYear*nSector,sizeof(double));


   //-----------------------------------------------
   // Prepare for simulation
   //-----------------------------------------------
   Setup();

   // Unpack paramters
   unpackPars(param, beta, sigma, gamma, xi, kappa, nu, lambda, phi);

   //-----------------------------------------------
   // Call simulation functions
   //-----------------------------------------------
   Simulation0(); // No shocks
   Simulation1(); // Price shock to manufacturing
   Simulation2(); // Unemployment shock to manufacturing
   Simulation3(); // Both shocks
   Simulation4(); // Both shocks, no utility cost of switching.
   Simulation5(); // Both shocks, half UI period.
   Simulation6(); // Both shocks, 50% UI coverage.
   Simulation7(); // Both shocks, lower UIbar.

   //-----------------------------------------------
   // Free allocated memory
   //-----------------------------------------------

   // Starting values
   free(param);
   free(beta);
   free(sigma);
   free(gamma);
   free(kappa);
   free(xi);
   free(nu);
   free(lambda);
   free(phi);
   free(delta);
   free(Delta);

   // Initial conditions
   free(CohortBirthData);
   free(ExperData);
   free(ChoiceData);
   free(LagSectorData);
   free(UnempYearsData);
   free(WageData);
   free(LagWageData);
   free(FirstYearData);
   free(AgeData);
   free(EducData);
   free(FemaleData);
   free(UIFundData);

   // Cohort sizes
   free(CohortSizeData);

   // Aggregate data
   free(OutputData);
   free(IncomeShares);
   free(CapitalData);
   free(rKData);

   
   return 0;
}


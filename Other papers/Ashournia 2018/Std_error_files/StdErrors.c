//------------------------------------------------------------------------------
//  StdErrors.c
//
//  By Damoun Ashournia, University of Oxford
//
//------------------------------------------------------------------------------
// Main module to compute standard errors.
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
#include <math.h>                      // C math library
#include <time.h>                      // In order to compute running time
#include <omp.h>                       // OpenMP

// GSL
#include <gsl/gsl_rng.h>               // Random number generator
#include <gsl/gsl_randist.h>           // distributions to be drawn from

// Mac OS X
#include <Accelerate/Accelerate.h>     // XCode Accelerate framework - holds lapack and blas


//------------------------------------------------------------------------------
// Declare globals
//------------------------------------------------------------------------------

// Global file pointer for StdErrors/params.txt
FILE *fpSEParams;

// Random Numbers
double *epsSim, *etaSim, *unif_type;

// Counter for objective function evaluations
int count;

// File pointers
FILE *fpParOut;

// Initial conditions
int *CohortBirthData, *FirstYearData, *AgeData, *EducData, *FemaleData, *UIFundData;
int *ChoiceData, *LagSectorData, *UnempYearsData, *EligData;
double *WageData, *LagWageData, *CostSim, *VSim, *ExperData;

// Unemployment transitions
double *delta;

// Starting values
double *param_start, *scaling;
int *mask, nParEst;

// Cohort weights
double *CohortWgt, *CohortSizeData;;

// Aggregate data
double *OutputData, *IncomeShares, *CapitalData, *rsk_global, *rsk_year;

// Auxiliary regressions
double *auxData, *covData, *invCov;

// For writing grid over params to file in loss function
int PAR_NUM;
double PAR_VAL;

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

// Used in Setup.c
#define iPopCohYear    iN + iCoh * nPop + iYear * nCoh * nPop
#define iPopCoh        iN + iCoh * nPop

// Used in LossFcn.c
#define iEps           iN + iCoh * nPop + iYear * nCoh * nPop + (iSec-1) * nYear * nCoh * nPop
#define iEta           iN + iCoh * nPop + iYear * nCoh * nPop + (iSec-1) * nYear * nCoh * nPop
#define iTyp           iN + iCoh * nPop
#define iDelta         (age-30) + female * nAge + educ * CatGender * nAge + lagSec * CatEduc * CatGender * nAge


//------------------------------------------------------------------------------
// Load user functions
//------------------------------------------------------------------------------
#include "CUtils.c"
#include "Setup.c"
#include "Emax.c"
#include "StdErrorUtils.c"

int main(int argc, const char * argv[])
{
    //-----------------------------------------------
    // Declarations
    //-----------------------------------------------
    
    // Counters
    int iParam, iN, iCoh, iYear, iSec, iCnt, k;
    
    // Input parameter vector
    double *param_in;
    
    // GSL Random Number Generator declarations
    const gsl_rng_type *T;
    gsl_rng *r;
    
    //-----------------------------------------------
    // Allocate memory
    //-----------------------------------------------
    
    // Random Numbers
    epsSim = calloc(nPop*nCoh*nYear*nSector,sizeof(double));
    etaSim = calloc(nPop*nCoh*nYear*nSector,sizeof(double));
    unif_type = calloc(nPop*nCoh,sizeof(double));

    // Starting values
    param_start = calloc(nParam,sizeof(double));
    param_in    = calloc(nParam,sizeof(double));
    scaling     = calloc(nParam,sizeof(double));
    mask        = calloc(nParam,sizeof(int));

    // Unemployment transitions
    delta       = calloc(CatGender*nAge*CatEduc*(nSector+1),sizeof(double));
    
    // Initial conditions
    CohortBirthData = calloc(nPop*nCoh*nYear,sizeof(int));
    ExperData       = calloc(nPop*nCoh*nYear,sizeof(double));
    ChoiceData      = calloc(nPop*nCoh*nYear,sizeof(int));
    LagSectorData   = calloc(nPop*nCoh*nYear,sizeof(int));
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
    // Read in data and parameters
    //-----------------------------------------------
    Setup();
    

    //-----------------------------------------------
    // Draw random numbers for simulations
    //-----------------------------------------------
    
    // Skill shock (normal)
    T = gsl_rng_mt19937;       // Algorithm
    r = gsl_rng_alloc(T);      // Allocate memory
    gsl_rng_set(r,5544332211); // Set seed

    for (iSec = 1; iSec < nSector + 1; iSec++)
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
    for (iSec = 1; iSec < nSector + 1; iSec++)
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
    
    // Beta
    for (iParam = 0; iParam < nParamBeta*nSector; iParam++)
    {
       mask[iParam] = 1;
    }

    // Sigma
    for (iSec = 0; iSec < nSector; iSec++)
    {
       mask[iSec + nParamBeta*nSector] = 1;
    }

    // Gamma
    for (iSec = 0; iSec < nSector + 1; iSec++)
    {
       mask[iSec + (nParamBeta+1)*nSector] = 1;
    }

    // Xi
    for (iSec = 0; iSec < 2*nSector - 1; iSec++)
    {
       mask[iSec + (nParamBeta+2)*nSector + 1] = 1;
    }

    // Kappa
    for (iParam = 0; iParam < 4; iParam++)
    {
       mask[iParam + (nParamBeta+4)*nSector] = 1;
    }

    // Nu
    mask[(nParamBeta+5)*nSector - 1] = 1;

    // Lambda
    for (iParam = 0; iParam < nSector; iParam++)
    {
       mask[iParam + (nParamBeta+5)*nSector] = 1;
    }

    // Phi
    for (iParam = 0; iParam < 4; iParam++)
    {
       mask[iParam + (nParamBeta+6)*nSector] = 1;
    }

    for (iParam = 0; iParam < nParam; iParam++) {
        nParEst += mask[iParam];
    }

    k = 0;
    for (iParam = 0; iParam < nParam; iParam++)
    {
        if (mask[iParam]==1)
        {
            param_in[k] = param_start[iParam];
            k += 1;
        }
    }

    // Number of parameters to compute standard errors for
    printf("Computing standard errors for %d parameter(s)\n", k);
    
    
    //-----------------------------------------------
    // Compute standard errors
    //-----------------------------------------------
    NumDerivatives(param_in);
    
    ComputeStdErrors();
    
    //-----------------------------------------------
    // Free allocated memory
    //-----------------------------------------------

    // Random numbers
    free(epsSim);
    free(etaSim);
    free(unif_type);
    
    // Parameter values
    free(param_start);
    free(param_in);
    free(scaling);
    free(mask);

    // Unemployment transitions
    free(delta);
    
    // Initial conditions
    free(CohortBirthData);
    free(ExperData);
    free(ChoiceData);
    free(LagSectorData);
    free(UnempYearsData);
    free(EligData);
    free(WageData);
    free(LagWageData);
    free(CostSim);
    free(FirstYearData);
    free(AgeData);
    free(EducData);
    free(FemaleData);
    free(UIFundData);
    free(VSim);
    
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
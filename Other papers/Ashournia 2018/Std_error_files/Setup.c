//------------------------------------------------------------------------------
// Setup.c
//
// By Damoun Ashournia, University of Oxford
//
//------------------------------------------------------------------------------
// Module for setting up the model for computing standard errors.
//------------------------------------------------------------------------------


void Setup(void)
{
//------------------------------------------------------------------
// Initilizations
//------------------------------------------------------------------

    // Counters
    int iSec, iCoh, iN, iYear, iParam, iParam2, iAge, iFemale, iEduc;
	
    // File pointers
    FILE *fpAux, *fpCov, *fpIni, *fpAgg, *fpWS, *fpPar, *fpCoh, *fpTr;

   // Temporary storage for data
   int n, birth, year, age, educ, female, uifund, exper, choice, lagchoice, unempyears, *CohSize, *CohCoh, *CohYear;
   double wage, lagwage, untr;

   // Allocate memory
   CohSize = calloc(nCohSizes,sizeof(int));
   CohCoh  = calloc(nCohSizes,sizeof(int));
   CohYear = calloc(nCohSizes,sizeof(int));


   //------------------------------------------------------------------
   // Read in starting values for parameter vector
   //------------------------------------------------------------------
   fpPar = fopen("Data/EstimationResults.txt", "r");

   if (fpPar == NULL)
   {
       perror("Failed to open file \"EstimationResults.txt\"");
       exit(EXIT_FAILURE);
   }

   for (iParam = 0; iParam < nParam; iParam++)
   {
      fscanf(fpPar, "%lf", &param_start[iParam]);
      // printf("%d\t%9.5f\n", iParam+1, param_start[iParam]);
   }

   fclose(fpPar);


   //------------------------------------------------------------------
   // Read in data
   //------------------------------------------------------------------
   
    // Initial conditions
    fpIni = fopen("Data/initial_conditions.txt", "r");
    
    if (fpIni == NULL)
    {
        perror("Failed to open file \"initial_conditions.txt\"");
        exit(EXIT_FAILURE);
    }
    
    for (iCoh = 0; iCoh < nCoh; iCoh++)
    {
        for (iN = 0; iN < nPop; iN++)
        {
            // %lf is 'long float' = double in fscanf
            fscanf(fpIni, "%d %d %d %d %d %d %d %d %d %d %lf %lf %d",
                   &n, &birth, &year, &age, &educ, &female, &uifund, &exper, &choice, &lagchoice, &wage, &lagwage, &unempyears);
            
            iYear = year - FirstYr;
            
            CohortBirthData[iPopCohYear] = birth;
            FirstYearData[iPopCoh]       = year;
            AgeData[iPopCoh]             = age;
            EducData[iPopCoh]            = educ;
            FemaleData[iPopCoh]          = female;
            UIFundData[iPopCoh]          = uifund;
            ExperData[iPopCohYear]       = exper;
            ChoiceData[iPopCohYear]      = choice;
            LagSectorData[iPopCohYear]   = lagchoice;
            WageData[iPopCohYear]        = wage;
            LagWageData[iPopCohYear]     = lagwage;
            UnempYearsData[iPopCohYear]  = unempyears;
        }
    }
    
    fclose(fpIni);

   // Cohort sizes
   fpCoh = fopen("Data/cohort_sizes.txt", "r");

   if (fpCoh == NULL)
   {
       perror("Failed to open file \"cohort_sizes.txt\"");
       exit(EXIT_FAILURE);
   }

   for (iN = 0; iN < nCohSizes; iN++)
   {
      fscanf(fpCoh, "%d %d %*d %d", &CohCoh[iN], &CohYear[iN], &CohSize[iN]);
   }

   for (iN = 0; iN < nCohSizes; iN+=6)
   {
      CohortSizeData[(CohCoh[iN] - FirstCoh) + (CohYear[iN] - FirstYr) * nCoh ] = CohSize[iN] + CohSize[iN + 1] + CohSize[iN + 2] + CohSize[iN + 3] + CohSize[iN + 4] + CohSize[iN + 5];
      // printf("%d\t%9.5f\n", iN+1, CohortSizeData[(CohCoh[iN] - FirstCoh) + (CohYear[iN] - FirstYr) * nCoh ]);
   }

   for (iYear = 0; iYear < nYear; iYear++)
   {
      for (iCoh = 0; iCoh < nCoh; iCoh++)
      {
         age = (iYear + 1996) - (iCoh + 1931);
         if (age <= LastAge && age >= FirstAge)
         {
            CohortWgt[iCoh + iYear * nCoh] = CohortSizeData[iCoh + iYear * nCoh] / CohortSizeData[0];
            // printf("%d\t%d\t%9.5f\n", iYear+1996, iCoh+1931, CohortWgt[iCoh + iYear * nCoh]);
         }
      }
   }

   fclose(fpCoh);


   // Output data
   fpAgg = fopen("Data/Real Value Added.txt", "r");

   if (fpAgg == NULL)
   {
       perror("Failed to open file \"Real Value Added.txt\"");
       exit(EXIT_FAILURE);
   }

   // Read first line (Variable names)
   char buffer1[20];
   fgets(buffer1, 20, fpAgg);

  // Read data
  for (iYear = 0; iYear < nYear; iYear++)
  {
    // %*lf means that fscanf does not load entry to memory
    fscanf(fpAgg, "%*d %lf %lf %lf %lf %lf", &OutputData[iYear + 0 * nYear], &OutputData[iYear + 1 * nYear], &OutputData[iYear + 2 * nYear], &OutputData[iYear + 3 * nYear], &OutputData[iYear + 4 * nYear]);
    for (iSec = 0; iSec < nSector; iSec++)
    {
      OutputData[iYear + iSec * nYear] = OutputData[iYear + iSec * nYear] * 1000000 / 1702; // 1702 = work-hours per year
      // printf("%d\t%d\t%9.5f\n", iSec+1, iYear+1996, OutputData[iYear + iSec * nYear]);
    }
  }

  fclose(fpAgg);


   // Capital data
   fpAgg = fopen("Data/Capital.txt", "r");

   if (fpAgg == NULL)
   {
       perror("Failed to open file \"Capital.txt\"");
       exit(EXIT_FAILURE);
   }

   // Read first line (Variable names)
   char buffer3[105];
   fgets(buffer3, 105, fpAgg);

   // Read data
   for (iYear = 0; iYear < nYear; iYear++)
   {
      // %*lf means that fscanf does not load entry to memory
      fscanf(fpAgg, "%*d %lf %lf %lf %lf %lf %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f", &CapitalData[0 + iYear * nSector], &CapitalData[1 + iYear * nSector], &CapitalData[2 + iYear * nSector], &CapitalData[3 + iYear * nSector], &CapitalData[4 + iYear * nSector]);
      for (iSec = 0; iSec < nSector; iSec++)
      {
        CapitalData[iSec + iYear * nSector] = CapitalData[iSec + iYear * nSector] * 1000000 / 1702; // 1702 = work-hours per year
        // printf("%d\t%d\t%9.5f\n", iSec+1, iYear+1996, CapitalData[iSec + iYear * nSector]);
      }
   }

   fclose(fpAgg);


   // Income shares
   fpWS = fopen("Data/Income Shares.txt", "r");

   // Read first line (Variable names)
   char buffer2[45];
   fgets(buffer2, 45, fpWS);

   if (fpWS == NULL)
   {
       perror("Failed to open file \"Income Shares.txt\"");
       exit(EXIT_FAILURE);
   }

  for (iYear = 0; iYear < nYear; iYear++)
  {
    fscanf(fpWS, "%*d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
                            &IncomeShares[iYear + 0 * nYear + 0 * nSector * nYear], &IncomeShares[iYear + 0 * nYear + 1 * nSector * nYear],
                            &IncomeShares[iYear + 1 * nYear + 0 * nSector * nYear], &IncomeShares[iYear + 1 * nYear + 1 * nSector * nYear],
                            &IncomeShares[iYear + 2 * nYear + 0 * nSector * nYear], &IncomeShares[iYear + 2 * nYear + 1 * nSector * nYear],
                            &IncomeShares[iYear + 3 * nYear + 0 * nSector * nYear], &IncomeShares[iYear + 3 * nYear + 1 * nSector * nYear],
                            &IncomeShares[iYear + 4 * nYear + 0 * nSector * nYear], &IncomeShares[iYear + 4 * nYear + 1 * nSector * nYear]);
    // for (iSec = 0; iSec < nSector; iSec++)
    // {
    //   printf("%d\t%d\t%9.5f\t%9.5f\n", iSec+1, iYear+1996, IncomeShares[iYear + iSec * nYear + 0 * nSector * nYear], IncomeShares[iYear + iSec * nYear + 1 * nSector * nYear]);
    // }
  }

  fclose(fpWS);


  // Unemployment transitions
  fpTr = fopen("Data/UnempTransitions.txt", "r");

  if (fpTr == NULL)
  {
      perror("Failed to open file \"UnempTransitions.txt\"");
      exit(EXIT_FAILURE);
  }

  for (iN = 0; iN < CatGender*nAge*CatEduc*(nSector+1); iN++)
  {
    fscanf(fpTr, "%d %d %d %d %lf", &female, &age, &educ, &lagchoice, &untr);
    iFemale = female;
    iAge = age-30;
    iEduc = educ;
    iSec = lagchoice;
    delta[iAge + iFemale * nAge + iEduc * CatGender * nAge + iSec * CatEduc * CatGender * nAge] = untr;
  }

  // printf("%9.5f\n", delta[22 + 1 * nAge + 1 * CatGender * nAge + 2 * CatEduc * CatGender * nAge]);

  fclose(fpTr);


	//------------------------------------------------------------------
	// Read coefficients from auxiliary regressions
	//------------------------------------------------------------------

	// Parameters
	fpAux = fopen("Data/AuxiliaryParams.txt", "r");

   if (fpAux == NULL)
   {
       perror("Failed to open file \"AuxiliaryParams.txt\"");
       exit(EXIT_FAILURE);
   }

   for (iParam = 0; iParam < nParamAux; iParam++)
   {
   	fscanf(fpAux, "%lf", &auxData[iParam]);
    // printf("%d\t%9.5f\n", iParam+1, auxData[iParam]);
   }

   fclose(fpAux);


   // Covariance matrix
   fpCov = fopen("Data/CovBoot.txt", "r");

   if (fpCov == NULL)
   {
       perror("Failed to open file \"CovBoot.txt\"");
       exit(EXIT_FAILURE);
   }

   for (iParam = 0; iParam < nParamAux; iParam++)
   {
    for (iParam2 = 0; iParam2 < nParamAux; iParam2++)
    {
      fscanf(fpCov, "%lf", &covData[iParam + iParam2 * nParamAux]);
    }
   }

   fclose(fpCov);

   // Inverse covariance matrix
   fcn_invX(covData, nParamAux, invCov);



   //------------------------------------------------------------------
   // Set rsk_global to average from log wage regressions
   //------------------------------------------------------------------
    for (iSec = 0; iSec < nSector; iSec++)
    {
        rsk_global[iSec] = exp( (auxData[8 + iSec * nRegWage]  + auxData[9 + iSec * nRegWage] +
                                 auxData[10 + iSec * nRegWage] + auxData[11 + iSec * nRegWage] +
                                 auxData[12 + iSec * nRegWage] + auxData[13 + iSec * nRegWage] +
                                 auxData[14 + iSec * nRegWage] + auxData[15 + iSec * nRegWage] +
                                 auxData[16 + iSec * nRegWage] + auxData[17 + iSec * nRegWage] +
                                 auxData[18 + iSec * nRegWage] + auxData[19 + iSec * nRegWage] +
                                 auxData[20 + iSec * nRegWage]) / nYear );

      for (iYear = 0; iYear < nYear; iYear++)
      {
        rsk_year[iYear + iSec * nYear] = auxData[iYear+8 + iSec * nYear];
      }
   }

   //------------------------------------------------------------------
   // Free allocated memory
   //------------------------------------------------------------------
   free(CohYear);
   free(CohCoh);
   free(CohSize);

}
//------------------------------------------------------------------------------
// Setup.c
//
// By Damoun Ashournia, University of Oxford
//
//------------------------------------------------------------------------------
// Module for setting up the model for estimation.
//------------------------------------------------------------------------------

void Setup(void)
{
    //------------------------------------------------------------------
    // Initilizations
    //------------------------------------------------------------------

    // Counters
    int iSec, iCoh, iN, iYear, iParam, iParam2;

    // File pointers
    FILE *fp;

    // Temporary storage for data
    int n, birth, year, age, educ, female, uifund, choice, lagchoice, unempyears, *CohSize, *CohCoh, *CohYear;
    double exper, wage, lagwage;

    // Allocate memory
    CohSize = calloc(nCohSizes,sizeof(int));
    CohCoh  = calloc(nCohSizes,sizeof(int));
    CohYear = calloc(nCohSizes,sizeof(int));


    //------------------------------------------------------------------
    // Read in starting values for parameter vector
    //------------------------------------------------------------------
    fp = fopen("Data/Starting_Values.txt", "r");

    if (fp == NULL)
    {
        perror("Failed to open file \"Starting_Values.txt\"");
        exit(EXIT_FAILURE);
    }

    for (iParam = 0; iParam < nParam; iParam++)
    {
        fscanf(fp, "%*s %lf", &param_start[iParam]);
    }

    fclose(fp);


    //------------------------------------------------------------------
    // Read in data
    //------------------------------------------------------------------

    // Initial conditions
    fp = fopen("Data/initial_conditions.txt", "r");

    if (fp == NULL)
    {
        perror("Failed to open file \"initial_conditions.txt\"");
        exit(EXIT_FAILURE);
    }

    for (iCoh = 0; iCoh < nCoh; iCoh++)
    {
        for (iN = 0; iN < nPop; iN++)
        {
            fscanf(fp, "%d %d %d %d %d %d %d %lf %d %d %lf %lf %d",
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
            LagChoiceData[iPopCohYear]   = lagchoice;
            WageData[iPopCohYear]        = wage;
            LagWageData[iPopCohYear]     = lagwage;
            UnempYearsData[iPopCohYear]  = unempyears;
        }
    }

    fclose(fp);


    // Cohort sizes
    fp = fopen("Data/cohort_sizes.txt", "r");

    if (fp == NULL)
    {
        perror("Failed to open file \"cohort_sizes.txt\"");
        exit(EXIT_FAILURE);
    }

    for (iN = 0; iN < nCohSizes; iN++)
    {
        fscanf(fp, "%d %d %*d %d", &CohCoh[iN], &CohYear[iN], &CohSize[iN]);
    }

    for (iN = 0; iN < nCohSizes; iN+=6)
    {
        CohortSizeData[(CohCoh[iN] - FirstCoh) + (CohYear[iN] - FirstYr) * nCoh ] = CohSize[iN] + CohSize[iN + 1] + CohSize[iN + 2] + CohSize[iN + 3] + CohSize[iN + 4] + CohSize[iN + 5];
    }

    for (iYear = 0; iYear < nYear; iYear++)
    {
        for (iCoh = 0; iCoh < nCoh; iCoh++)
        {
            age = (iYear + 1996) - (iCoh + 1931);
            if (age <= LastAge && age >= FirstAge)
            {
                CohortWgt[iCoh + iYear * nCoh] = CohortSizeData[iCoh + iYear * nCoh] / CohortSizeData[0];
            }
        }
    }

    fclose(fp);


    // Output data
    fp = fopen("Data/Real Value Added.txt", "r");

    if (fp == NULL)
    {
        perror("Failed to open file \"Real Value Added.txt\"");
        exit(EXIT_FAILURE);
    }

    // Read first line (Variable names)
    char buffer1[20];
    fgets(buffer1, 20, fp);

    // Read data
    for (iYear = 0; iYear < nYear; iYear++)
    {
        fscanf(fp, "%*d %lf %lf %lf %lf %lf", &OutputData[iYear + 0 * nYear], &OutputData[iYear + 1 * nYear], &OutputData[iYear + 2 * nYear], &OutputData[iYear + 3 * nYear], &OutputData[iYear + 4 * nYear]);
        for (iSec = 0; iSec < nSector; iSec++)
        {
            OutputData[iYear + iSec * nYear] = OutputData[iYear + iSec * nYear] * 1000000 / 1702; // 1702 = work-hours per year
        }
    }

    fclose(fp);


    // Capital data
    fp = fopen("Data/Capital.txt", "r");

    if (fp == NULL)
    {
        perror("Failed to open file \"Capital.txt\"");
        exit(EXIT_FAILURE);
    }

    // Read first line (Variable names)
    char buffer3[105];
    fgets(buffer3, 105, fp);

    // Read data
    for (iYear = 0; iYear < nYear; iYear++)
    {
        fscanf(fp, "%*d %lf %lf %lf %lf %lf %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f", &CapitalData[0 + iYear * nSector], &CapitalData[1 + iYear * nSector], &CapitalData[2 + iYear * nSector], &CapitalData[3 + iYear * nSector], &CapitalData[4 + iYear * nSector]);
        for (iSec = 0; iSec < nSector; iSec++)
        {
            CapitalData[iSec + iYear * nSector] = CapitalData[iSec + iYear * nSector] * 1000000 / 1702; // 1702 = work-hours per year
        }
    }

    fclose(fp);


    // Income shares
    fp = fopen("Data/Income Shares.txt", "r");

    // Read first line (Variable names)
    char buffer2[45];
    fgets(buffer2, 45, fp);

    if (fp == NULL)
    {
        perror("Failed to open file \"Income Shares.txt\"");
        exit(EXIT_FAILURE);
    }

    for (iYear = 0; iYear < nYear; iYear++)
    {
        fscanf(fp, "%*d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &IncomeShares[iYear + 0 * nYear + 0 * nSector * nYear], &IncomeShares[iYear + 0 * nYear + 1 * nSector * nYear],
                &IncomeShares[iYear + 1 * nYear + 0 * nSector * nYear], &IncomeShares[iYear + 1 * nYear + 1 * nSector * nYear],
                &IncomeShares[iYear + 2 * nYear + 0 * nSector * nYear], &IncomeShares[iYear + 2 * nYear + 1 * nSector * nYear],
                &IncomeShares[iYear + 3 * nYear + 0 * nSector * nYear], &IncomeShares[iYear + 3 * nYear + 1 * nSector * nYear],
                &IncomeShares[iYear + 4 * nYear + 0 * nSector * nYear], &IncomeShares[iYear + 4 * nYear + 1 * nSector * nYear]);
    }

    fclose(fp);

    //------------------------------------------------------------------
    // Read coefficients from auxiliary regressions
    //------------------------------------------------------------------

    // Parameters
    fp = fopen("Data/AuxiliaryParams.txt", "r");

    if (fp == NULL)
    {
        perror("Failed to open file \"AuxiliaryParams.txt\"");
        exit(EXIT_FAILURE);
    }

    for (iParam = 0; iParam < nParamAux; iParam++)
    {
        fscanf(fp, "%lf", &auxData[iParam]);
    }

    fclose(fp);


    // Covariance matrix
    fp = fopen("Data/CovBoot.txt", "r");

    if (fp == NULL)
    {
        perror("Failed to open file \"CovBoot.txt\"");
        exit(EXIT_FAILURE);
    }

    for (iParam = 0; iParam < nParamAux; iParam++)
    {
        for (iParam2 = 0; iParam2 < nParamAux; iParam2++)
        {
            fscanf(fp, "%lf", &covData[iParam + iParam2 * nParamAux]);
        }
    }

    fclose(fp);

    // Inverse covariance matrix
    fcn_invX(covData, nParamAux, invCov);



    //------------------------------------------------------------------
    // Set rsk_global to average from log wage regressions
    //------------------------------------------------------------------
    for (iSec = 0; iSec < nSector; iSec++)
    {
        rsk_global[iSec] = exp( (auxData[5  + iSec * nRegWage] + auxData[6  + iSec * nRegWage] +
                                auxData[7  + iSec * nRegWage] + auxData[8  + iSec * nRegWage] +
                                auxData[9  + iSec * nRegWage] + auxData[10 + iSec * nRegWage] +
                                auxData[11 + iSec * nRegWage] + auxData[12 + iSec * nRegWage] +
                                auxData[13 + iSec * nRegWage] + auxData[14 + iSec * nRegWage] +
                                auxData[15 + iSec * nRegWage] + auxData[16 + iSec * nRegWage] +
                                auxData[17 + iSec * nRegWage]) / nYear );
    }

    //------------------------------------------------------------------
    // Free allocated memory
    //------------------------------------------------------------------
    free(CohYear);
    free(CohCoh);
    free(CohSize);

}

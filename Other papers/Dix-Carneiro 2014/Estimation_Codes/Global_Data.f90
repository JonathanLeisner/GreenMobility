Module Global_Data

! ********************************************
! Global variables used throughout the program
! ********************************************

implicit none

save

! Double Precision
integer, parameter :: DOUBLE      = SELECTED_REAL_KIND(p=10)

integer, parameter :: NSECTORS    = 7

! Number of parameters of the model
integer, parameter :: NPARAM      = 223

! Number of types in the model
integer, parameter :: NTYPES      = 3

! Number of Regressors used in the auxiliary models
integer, parameter :: NREG_Wage    = 24
integer, parameter :: NREG_Emp     = 24
integer, parameter :: NREG_Return  = 24
integer, parameter :: NREG_Tr      = 24
integer, parameter :: NREG_WageDif = 11
integer, parameter :: NREG_Pers    = 21
integer, parameter :: NREG_Freq    = 21


! **********************************
! Numerical Approximation Parameters
! **********************************

integer, parameter :: NREG       = 134

! # of Monte Carlo Integration Draws
integer, parameter :: NDRAWS     = 300
integer, parameter :: WARMUP     = 100

! # of Points Where the Value Function is Calculated 
! via Monte Carlo Integration
integer, parameter :: INTP       = 600

! # of simulated people per cohort
integer, parameter :: NPOP       = 2000

! Maximum number of iterations for the computation
! of equilibrium returns to skill
integer, parameter :: MAXIT      = 40



integer, parameter :: FIRST_AGE  = 25
integer, parameter :: LAST_AGE   = 60
integer, parameter :: NPER       = LAST_AGE - FIRST_AGE + 1


integer, parameter :: FIRST_YR      = 1995
integer, parameter :: LAST_YR       = 2005
integer, parameter :: FIRST_YR_EMAX = 1996
integer, parameter :: LAST_YR_EMAX  = 2040
integer, parameter :: FIRST_COH     = FIRST_YR - LAST_AGE
integer, parameter :: LAST_COH      = LAST_YR - FIRST_AGE
integer, parameter :: NGEN          = LAST_COH - FIRST_COH + 1

integer, parameter :: CAT_EDUC   = 4
integer, parameter :: CAT_GENDER = 2

integer, parameter :: COH_PROC = 3
integer, parameter :: NCOH_BATCH = 12

! Variables used in counterfactual simulation
integer, parameter :: NPOP_SIM      = 20000
integer, parameter :: SIM_T         = 300
integer, parameter :: FIRST_YR_SIM  = LAST_YR
integer, parameter :: LAST_YR_SIM   = LAST_YR + SIM_T
integer, parameter :: FIRST_COH_SIM = LAST_YR - LAST_AGE 
integer, parameter :: LAST_COH_SIM  = LAST_YR + SIM_T - FIRST_AGE

real(KIND=DOUBLE) param_global(NPARAM)
integer mask_global(NPARAM)

! Discount factor
real, parameter :: rho = 0.95

! Range of returns to skill I sample from
! In order to get the points where I compute
! the value function exactly
real(KIND=DOUBLE) LOW_RSK(NSECTORS,0:1), &
                  UP_RSK(NSECTORS,0:1)

real(KIND=DOUBLE) eps_global(NPOP,FIRST_COH:LAST_COH,FIRST_YR:LAST_YR,0:NSECTORS), &
                  eta_global(NPOP,FIRST_COH:LAST_COH,FIRST_YR:LAST_YR,0:NSECTORS), &
                  unif_global(NPOP,FIRST_COH:LAST_COH,1), &
                  rsk_global(NSECTORS,0:1), &
                  ExpansionFactorGlobal, &
                  scaling_global(NPARAM)
         
integer           FUNCTION_ITER, &
                  ALGORITHM, &
                  BATCH

integer           StdErrors_FLAG, &
                  PARAM_NUMBER
                  
real(KIND=DOUBLE) PARAM_VALUE


! Input Data

integer FirstYearData       (NPOP,FIRST_COH:LAST_COH),                                      &
        AgeData             (NPOP,FIRST_COH:LAST_COH),                                      &
        EducData            (NPOP,FIRST_COH:LAST_COH),                                      &
        GenderData          (NPOP,FIRST_COH:LAST_COH),                                      &       
        ChoiceData          (NPOP,FIRST_COH:LAST_COH,FIRST_YR-9:LAST_YR),                   &                       
        CohortSizeData      (FIRST_COH:LAST_COH,0:1),                                       &
        EducData_ExpFactor  (NPOP,FIRST_COH:LAST_COH,FIRST_YR:LAST_YR)    
        
               
real(KIND=DOUBLE) WageData            (NPOP,FIRST_COH:LAST_COH,FIRST_YR:LAST_YR),           &        
                  OutputData          (NSECTORS,FIRST_YR:LAST_YR),                          &
                  CapitalData         (FIRST_YR:LAST_YR),                                   &
                  rKData              (FIRST_YR:LAST_YR),                                   &
                  TotalHrWageData     (1:NSECTORS,FIRST_YR:LAST_YR,1:CAT_EDUC),             &
                  TotalWageData       (1:NSECTORS,FIRST_YR:LAST_YR,1:CAT_EDUC),             &
                  LaborSharesData     (1:NSECTORS,0:1,FIRST_YR:LAST_YR),                    &
                  CapitalSharesData   (1:NSECTORS,FIRST_YR:LAST_YR),                        &
                  WageData_ExpFactor  (NPOP,FIRST_COH:LAST_COH,FIRST_YR:LAST_YR) 
        
real(KIND=DOUBLE) betaData(NSECTORS,NREG_Wage),                                             &
                  CPWageData(NSECTORS,NREG_Wage,NREG_Wage),                                 & 
                  COVbetaData(NSECTORS,NREG_Wage,NREG_Wage),                                &
                  invCOVbetaData(NSECTORS,NREG_Wage,NREG_Wage),                             &
                  sigmaWageData(NSECTORS),                                                  &
                  gammaData(0:NSECTORS,NREG_Emp),                                           &
                  CPEmpData(0:NSECTORS,NREG_Emp,NREG_Emp),                                  &
                  COVgammaData(0:NSECTORS,NREG_Emp,NREG_Emp),                               &
                  invCOVgammaData(0:NSECTORS,NREG_Emp,NREG_Emp),                            &
                  sigmaEmpData(0:NSECTORS),                                                 &
                  phiData(0:NSECTORS,0:NSECTORS,NREG_Tr),                                   &
                  CPTrData(0:NSECTORS,0:NSECTORS,NREG_Tr,NREG_Tr),                          &
                  COVphiData(0:NSECTORS,0:NSECTORS,NREG_Tr,NREG_Tr),                        &
                  invCOVphiData(0:NSECTORS,0:NSECTORS,NREG_Tr,NREG_Tr),                     &
                  sigmaTrData(0:NSECTORS,0:NSECTORS),                                       &
                  VarSigma2Data(NSECTORS),                                                  &
                  sigmaWageDifData(NSECTORS),                                               &
                  VarSigma2WageDifData(NSECTORS),                                           &
                  xsi1998Data(0:NSECTORS,NREG_Pers),                                        &
                  CPPers1998Data(0:NSECTORS,NREG_Pers,NREG_Pers),                           &
                  COVxsi1998Data(0:NSECTORS,NREG_Pers,NREG_Pers),                           &
                  invCOVxsi1998Data(0:NSECTORS,NREG_Pers,NREG_Pers),                        &
                  sigmaPers1998Data(0:NSECTORS),                                            &
                  rhoData(0:NSECTORS,NREG_Return),                                          &
                  CPReturnData(0:NSECTORS,NREG_Return,NREG_Return),                         &
                  COVrhoData(0:NSECTORS,NREG_Return,NREG_Return),                           &
                  invCOVrhoData(0:NSECTORS,NREG_Emp,NREG_Emp),                              &
                  sigmaReturnData(0:NSECTORS),                                              &
                                  
                  xsi2000Data(0:NSECTORS,NREG_Pers),                                        &
                  CPPers2000Data(0:NSECTORS,NREG_Pers,NREG_Pers),                           &
                  COVxsi2000Data(0:NSECTORS,NREG_Pers,NREG_Pers),                           &
                  invCOVxsi2000Data(0:NSECTORS,NREG_Pers,NREG_Pers),                        &
                  sigmaPers2000Data(0:NSECTORS),                                            &
                  
                  xsi2005Data(0:NSECTORS,NREG_Pers),                                        &
                  CPPers2005Data(0:NSECTORS,NREG_Pers,NREG_Pers),                           &
                  COVxsi2005Data(0:NSECTORS,NREG_Pers,NREG_Pers),                           &
                  invCOVxsi2005Data(0:NSECTORS,NREG_Pers,NREG_Pers),                        &
                  sigmaPers2005Data(0:NSECTORS),                                            &
                  
                  etaData(0:NSECTORS,NREG_Freq),                                            &
                  CPFreqData(0:NSECTORS,NREG_Freq,NREG_Freq),                               &
                  COVetaData(0:NSECTORS,NREG_Freq,NREG_Freq),                               &
                  invCOVetaData(0:NSECTORS,NREG_Freq,NREG_Freq),                            &
                  sigmaFreqData(0:NSECTORS),                                                &                  
                  
                  NobsWageData(NSECTORS),                                                   &
                  NobsEmpData(0:NSECTORS),                                                  &
                  NobsTrData(0:NSECTORS,0:NSECTORS),                                        &         
                  NobsWageDifData(NSECTORS),                                                &
                  NobsReturnData(0:NSECTORS),                                               &

                  NobsPers1998Data(0:NSECTORS),                                             &
                  NobsPers2000Data(0:NSECTORS),                                             &
                  NobsPers2005Data(0:NSECTORS),                                             &
                  
                  NobsFreqData(0:NSECTORS)                                                
                  

end module Global_Data

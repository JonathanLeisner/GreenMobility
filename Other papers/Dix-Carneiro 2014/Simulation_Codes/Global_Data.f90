Module Global_Data

! ********************************************
! Global variables used throughout the program
! ********************************************

implicit none

save

! Double Precision
integer, parameter :: DOUBLE      = SELECTED_REAL_KIND(p=10)

integer, parameter :: NSECTORS    = 7

integer, parameter :: NTYPES    = 3

! Number of parameters of the model
integer, parameter :: NPARAM      = 223

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
integer, parameter :: MAXIT      = 80


integer, parameter :: FIRST_AGE  = 25
integer, parameter :: LAST_AGE   = 60
integer, parameter :: NPER       = LAST_AGE - FIRST_AGE + 1


integer, parameter :: FIRST_YR   = 1995
integer, parameter :: LAST_YR    = 2005
integer, parameter :: FIRST_COH  = FIRST_YR - LAST_AGE
integer, parameter :: LAST_COH   = LAST_YR - FIRST_AGE
integer, parameter :: NGEN       = LAST_COH - FIRST_COH + 1

integer, parameter :: CAT_EDUC   = 4
integer, parameter :: CAT_GENDER = 2

integer, parameter :: COH_PROC = 3
integer, parameter :: NCOH_BATCH = 12


! Variables used in counterfactual simulation
integer, parameter :: NPOP_SIM      = 70000
integer, parameter :: SIM_T         = 300
integer, parameter :: FIRST_YR_SIM  = LAST_YR
integer, parameter :: LAST_YR_SIM   = LAST_YR + SIM_T
integer, parameter :: FIRST_COH_SIM = LAST_YR - LAST_AGE 
integer, parameter :: LAST_COH_SIM  = LAST_YR + SIM_T - FIRST_AGE
integer, parameter :: HALF_YR_SIM   = floor((FIRST_YR_SIM+LAST_YR_SIM)/2.0)
integer, parameter :: YR_SHOCK      = HALF_YR_SIM + 1
integer, parameter :: YR_ANNMNT     = YR_SHOCK

real(KIND=DOUBLE) param_global(NPARAM)
integer mask_global(NPARAM)

! Discount factor
real, parameter :: rho = 0.95

real(KIND=DOUBLE) COST_MULT

! Range of returns to skill I sample from
! In order to get the points where I compute
! the value function exactly
real(KIND=DOUBLE) LOW_RSK(NSECTORS,0:1), &
                  UP_RSK(NSECTORS,0:1)

real(KIND=DOUBLE) eps_global(NPOP,FIRST_COH:LAST_COH,FIRST_YR:LAST_YR,0:NSECTORS),    &
                  rsk_global(NSECTORS,0:1), &
                  ExpansionFactorGlobal, &
                  scaling_global(NPARAM)
         
integer           FUNCTION_ITER, &
                  ALGORITHM, &
                  BATCH, &
                  CAPMOBILITY

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
                  aaGlobal            (1:NSECTORS,0:1,FIRST_YR:LAST_YR),                    &
                  WageData_ExpFactor  (NPOP,FIRST_COH:LAST_COH,FIRST_YR:LAST_YR) 
        
end module Global_Data

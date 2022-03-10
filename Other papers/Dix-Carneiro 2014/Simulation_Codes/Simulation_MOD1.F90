! ************************
! Perfect Capital Mobility
! ************************

Module Simulation_MOD1



Contains


Subroutine Simulation1(param)

USE Global_data
USE Emax_MOD
USE WelfareComp
USE Newton
USE ParallelCohorts_MOD
USE MPI

implicit none

integer, parameter :: dp = SELECTED_REAL_KIND(12,60)
! Affected sector: High-Tech = 3
integer, parameter :: AFFSECTOR = 3

real (dp) , intent(in)  :: param(NPARAM)

integer i, j, n, n2, t, s, it, it2, coh, gen, year, age, gender, & 
        educ, SimExpansion, status2, ISEED, &
        educ2, educ3, educ4, ed, RE_loop, tp, &
        s1, s2, ii, ind, coh_batch
        

! Model parameters
real(KIND=DOUBLE) theta(7), beta(NSECTORS,12), & 
                  sigma(0:NSECTORS), kappa(26), tau(0:NSECTORS), SigmaPref, &
                  omega(2:NTYPES,0:NSECTORS), lambda(NTYPES), gamma(2:NTYPES,0:NSECTORS+4), &
                  aa(NSECTORS,0:1,0:1), alpha_prod(NSECTORS,0:1,FIRST_YR:LAST_YR_SIM), sigma_prod(NSECTORS), &
                  rsk_init(NSECTORS,0:1), alpha_c(1:NSECTORS)     

real(KIND=DOUBLE) price(NSECTORS,FIRST_YR:LAST_YR_SIM), &
                  prod(NSECTORS,FIRST_YR:LAST_YR_SIM), &
                  Pindex(FIRST_YR:LAST_YR_SIM), &
                  rvec(0:NSECTORS)
         
! Iteration Loop for the Equilibrium 
real(KIND=DOUBLE) rsk_sup(NSECTORS,0:1,MAXIT), &
                  rsk_dem(NSECTORS,0:1,MAXIT), &
                  rsk_eq(NSECTORS,0:1,FIRST_YR:LAST_YR_SIM), &
                  avec(NSECTORS,0:1,FIRST_YR_SIM:LAST_YR_SIM), &
                  rsk_tomorrow(NSECTORS), &
                  sk_sup(NSECTORS,0:1), &
                  sk_sup_eq(NSECTORS,0:1,FIRST_YR:LAST_YR_SIM), &
                  cost(0:NSECTORS,0:NSECTORS), & 
                  cost25(0:NSECTORS,0:NSECTORS), & 
                  Emax(0:NSECTORS), & 
                  PI_COEF(NREG), &
                  cut(NSECTORS,0:1), &
                  w(0:NSECTORS), &
                  V(0:NSECTORS), &
                  VMAX, &
                  check1(NSECTORS,0:1), &
                  check2, &
                  step(NSECTORS,0:1), &
                  step2, & 
                  CohortWgt(FIRST_COH:LAST_COH_SIM,0:1), &
                  weight, &
                  eps_aux(0:NSECTORS), &
                  eta_aux(0:NSECTORS), &
                  CapitalStock, &
                  rK_sup(MAXIT), &
                  rK_dem(MAXIT), &
                  rK_eq(FIRST_YR:LAST_YR_SIM), &
                  output(NSECTORS,FIRST_YR:LAST_YR_SIM), &
                  crit1, &
                  crit2, &
                  emp(0:NSECTORS,MAXIT), & 
                  emp_eq(0:NSECTORS,FIRST_YR:LAST_YR_SIM), &
                  emp_by_skill(0:NSECTORS,0:1,FIRST_YR:LAST_YR_SIM), &
                  Welfare(FIRST_YR_SIM:LAST_YR_SIM), &
                  Welfare_Educ(0:1,FIRST_YR_SIM:LAST_YR_SIM), &
                  Real_Wages(FIRST_YR_SIM:LAST_YR_SIM), &
                  Real_Wages_Educ(0:1,FIRST_YR_SIM:LAST_YR_SIM), &
                  term2, &
                  term3, &
                  SumW(CAT_GENDER,CAT_EDUC), &
                  AllocEmp(NSECTORS,NTYPES,CAT_GENDER,CAT_EDUC), &
                  Const1, &
                  Const2, &
                  Const3, &
                  init, &
                  solution
                  
real(KIND=DOUBLE) A(FIRST_AGE+1:LAST_AGE,0:NSECTORS,NREG), &
                  B(FIRST_AGE+1:LAST_AGE,0:NSECTORS), &
                  A_RE((YR_ANNMNT+1):(LAST_YR_SIM-1),FIRST_AGE+1:LAST_AGE,0:NSECTORS,NREG), &
                  B_RE((YR_ANNMNT+1):(LAST_YR_SIM-1),FIRST_AGE+1:LAST_AGE,0:NSECTORS)


integer           iter(FIRST_YR:LAST_YR_SIM), &
                  iter2(FIRST_YR:LAST_YR_SIM), &
                  flag(NSECTORS,0:1), &
                  flag2, &
                  lag, &
                  exper(NSECTORS), &
                  ExperTomorrow(NSECTORS)
                  
logical           converged, &
                  converged2                 

real(KIND=DOUBLE) Cap(NSECTORS), &
                  CapitalSim(NSECTORS,FIRST_YR:LAST_YR_SIM), &
                  z(NSECTORS,FIRST_YR:LAST_YR_SIM)
                  

real(KIND=DOUBLE) output_inf1, output_inf2, output_A, Gains1, Gains2, &
                  AdjCost, perc_complete(FIRST_YR_SIM:LAST_YR_SIM), &
                  tot_output(FIRST_YR_SIM:LAST_YR_SIM), &
                  perc_complete2(FIRST_YR_SIM:LAST_YR_SIM), &
                  emp1, emp2, &
                  Welfare_inf1, Welfare_inf2, Welfare_A, WGains1, WGains2, &
                  WAdjCost, &
                  Wage_inf1, Wage_inf2, Wage_A, WageGains1, WageGains2, &
                  WageAdjCost, &
                  Welfare_inf_Educ1, Welfare_inf_Educ2, Welfare_EducA, &
                  WGains_Educ1, WGains_Educ2, WAdjCost_Educ, &
                  Welfare_inf_NonEduc1, Welfare_inf_NonEduc2, Welfare_NonEducA, &
                  WGains_NonEduc1, WGains_NonEduc2, WAdjCost_NonEduc, &
                  Wage_inf_Educ1, Wage_inf_Educ2, Wage_EducA, &
                  WageGains_Educ1, WageGains_Educ2, WageAdjCost_Educ, &
                  Wage_inf_NonEduc1, Wage_inf_NonEduc2, Wage_NonEducA, &
                  WageGains_NonEduc1, WageGains_NonEduc2, WageAdjCost_NonEduc
                  
real(KIND=DOUBLE) p_sup(4:7,MAXIT), p_dem(4:7,MAXIT), check3(4:7), step3(4:7), crit3
integer flag3                   
                  
integer Size2
integer, allocatable, dimension(:) :: n_index, coh_index
integer, allocatable, dimension(:,:,:) :: ExperSim
real(KIND=DOUBLE), allocatable, dimension(:,:,:) :: WelfareSim, WageSim
                  
integer          , allocatable, dimension(:,:)   :: EducSim, GenderSim, typSim, EducInit, GenderInit, typInit
integer          , allocatable, dimension(:,:,:) :: ChoiceSim, ChoiceInit

! Coefficients for the calculation of Emax_hat
real(KIND=DOUBLE), allocatable, dimension(:,:,:,:,:,:)   :: PI_Myopic
real(KIND=DOUBLE), allocatable, dimension(:,:,:,:,:)     :: R_sq_Myopic
real(KIND=DOUBLE), allocatable, dimension(:,:,:,:,:,:,:) :: PI_RE
real(KIND=DOUBLE), allocatable, dimension(:,:,:,:,:,:)   :: R_sq_RE
real(KIND=DOUBLE), allocatable, dimension(:,:,:)         :: eps, eta    


real(KIND=DOUBLE) epsAux(NPOP_SIM,COH_PROC,0:NSECTORS), &
                  etaAux(NPOP_SIM,COH_PROC,0:NSECTORS), &
                  WelfareSimOut(NPOP_SIM,COH_PROC), &
                  WageSimOut(NPOP_SIM,COH_PROC), &
                  CohortWgtAux(COH_PROC,0:1), &
                  sk_supOut(NSECTORS,0:1), &
                  WelfareOut, &
                  Welfare_EducOut(0:1), &
                  Real_WagesOut, &
                  Real_Wages_EducOut(0:1), &
                  empOut(0:NSECTORS), &
                  emp_by_skillOut(0:NSECTORS,0:1)

integer older_coh, &
        newer_coh, &
        older_gen, &
        newer_gen, &
        EducSimAux(NPOP_SIM,COH_PROC), &
        GenderSimAux(NPOP_SIM,COH_PROC), &
        ChoiceSimAux(NPOP_SIM,COH_PROC,-9:-1), &
        ChoiceSimOut(NPOP_SIM,COH_PROC), &
        ExperSimOut(NPOP_SIM,COH_PROC)
           
     
! Random Number Generation
TYPE (VSL_STREAM_STATE) :: stream_gauss, stream_gumbel, stream_unif
integer(kind=4)   errcode
integer           brng, method_gauss, method_gumbel, method_unif, seed
real(KIND=DOUBLE) rvec_eps(0:NSECTORS), rvec_eta(0:NSECTORS), unif(1), scale, loc

! MPI Variables
integer ierr, numproc, rank, tag, source, count1, count2, proc, ierror
integer status(MPI_STATUS_SIZE)


! ***************************************
! MPI functions for processor information
! ***************************************

Call MPI_COMM_SIZE (MPI_COMM_WORLD,numproc,ierr)
Call MPI_COMM_RANK (MPI_COMM_WORLD,rank,ierr)

allocate(PI_Myopic(FIRST_AGE+1:LAST_AGE,0:NSECTORS,1:CAT_GENDER,1:CAT_EDUC,1:NTYPES,NREG))
allocate(R_sq_Myopic(FIRST_AGE+1:LAST_AGE,0:NSECTORS,1:CAT_GENDER,1:CAT_EDUC,1:NTYPES))   

allocate(PI_RE(YR_ANNMNT+1:LAST_YR_SIM-1,FIRST_AGE+1:LAST_AGE,0:NSECTORS,1:CAT_GENDER,1:CAT_EDUC,1:NTYPES,NREG))
allocate(R_sq_RE(YR_ANNMNT+1:LAST_YR_SIM-1,FIRST_AGE+1:LAST_AGE,0:NSECTORS,1:CAT_GENDER,1:CAT_EDUC,1:NTYPES))

allocate(typSim(NPOP_SIM,FIRST_COH:LAST_COH_SIM))

CAPMOBILITY = 1

! *****************************
! Use only the Master processor
! *****************************

if (rank == 0) then    
         
    allocate(ChoiceInit(NPOP,(LAST_YR-LAST_AGE):(LAST_YR-FIRST_AGE),(FIRST_AGE-9):LAST_AGE))
    allocate(ChoiceSim(NPOP_SIM,FIRST_COH:LAST_COH_SIM,(FIRST_AGE-9):LAST_AGE))

    allocate(EducInit(NPOP,(LAST_YR-LAST_AGE):(LAST_YR-FIRST_AGE)))
    allocate(EducSim(NPOP_SIM,FIRST_COH:LAST_COH_SIM))

    allocate(GenderInit(NPOP,(LAST_YR-LAST_AGE):(LAST_YR-FIRST_AGE)))
    allocate(GenderSim(NPOP_SIM,FIRST_COH:LAST_COH_SIM))

    allocate(typInit(NPOP,(LAST_YR-LAST_AGE):(LAST_YR-FIRST_AGE)))

    allocate(eps(NPOP_SIM,36,0:NSECTORS))        
    allocate(eta(NPOP_SIM,36,0:NSECTORS))

    theta(1:7)      = param(1:7)
    beta(1,1:12)    = param(8:19)
    beta(2,1:12)    = param(20:31)
    beta(3,1:12)    = param(32:43)
    beta(4,1:12)    = param(44:55)
    beta(5,1:12)    = param(56:67)
    beta(6,1:12)    = param(68:79)
    beta(7,1:12)    = param(80:91)
    sigma(0:7)      = param(92:99)
    kappa(1:26)     = param(100:125)
    tau             = 0.0
    tau(2:7)        = param(126:131)
    SigmaPref       = param(132)
    omega(2,0:7)    = param(133:140)
    omega(3,0:7)    = param(141:148)
    lambda          = 0.0
    lambda(2:3)     = param(149:150)
    gamma(2,0:11)   = param(151:162)
    gamma(3,0:11)   = param(163:174)
    aa(1,0,0:1)     = param(175:176)
    aa(2,0,0:1)     = param(177:178)
    aa(3,0,0:1)     = param(179:180)
    aa(4,0,0:1)     = param(181:182)
    aa(5,0,0:1)     = param(183:184)
    aa(6,0,0:1)     = param(185:186)
    aa(7,0,0:1)     = param(187:188)
    aa(1,1,0:1)     = param(189:190)
    aa(2,1,0:1)     = param(191:192)
    aa(3,1,0:1)     = param(193:194)
    aa(4,1,0:1)     = param(195:196)
    aa(5,1,0:1)     = param(197:198)
    aa(6,1,0:1)     = param(199:200)
    aa(7,1,0:1)     = param(201:202)
    sigma_prod(1:7) = param(203:209)
    rsk_init(1:7,0) = param(210:216)
    rsk_init(1:7,1) = param(217:223)       


    ! ******************
    ! Consumption Shares
    ! ******************

    alpha_c(1) = 0.08
    alpha_c(2) = 0.07
    alpha_c(3) = 0.09
    alpha_c(4) = 0.05
    alpha_c(5) = 0.11
    alpha_c(6) = 0.12
    alpha_c(7) = 0.48


    do ed = 0, 1
        do s = 1, NSECTORS
            do t = FIRST_YR, LAST_YR
                alpha_prod(s,ed,t) = aa(s,ed,0) + aa(s,ed,1)*(t-FIRST_YR)
            end do
        end do
    end do


    ! ********************
    ! Convergence Criteria
    ! ********************

    crit1 = 0.00005
    crit2 = 0.00005
    crit3 = 0.00005


    ! *****************************************
    ! Initializing with "missing values" = -999
    ! *****************************************

    EducSim    = -999
    EducInit   = -999
    GenderSim  = -999
    GenderInit = -999
    ChoiceSim  = -999
    ChoiceInit = -999
    CohortWgt  = -999


    ! ******************************************************************
    ! Read Initial Conditions from file generated in the estimation code
    ! ******************************************************************

    open(unit = 1, file = 'Init_Simulation.csv')
    do 
        read(1,4321,IOSTAT=status2) n, t, coh, EducInit(n,coh), GenderInit(n,coh), &
                                    ChoiceInit(n,coh,t-coh)   , ChoiceInit(n,coh,t-coh-1), &
                                    ChoiceInit(n,coh,t-coh-2) , ChoiceInit(n,coh,t-coh-3), &
                                    ChoiceInit(n,coh,t-coh-4) , ChoiceInit(n,coh,t-coh-5), & 
                                    ChoiceInit(n,coh,t-coh-6) , ChoiceInit(n,coh,t-coh-7), &
                                    ChoiceInit(n,coh,t-coh-8) , ChoiceInit(n,coh,t-coh-9), &
                                    ed, typInit(n,coh), CohortWgt(coh,ed)
        if (status2 /= 0 ) exit
    end do
    close(1)

    4321 format(17(i5,','),f9.4,',')

    AllocEmp = 0.0
    do coh = FIRST_COH_SIM, LAST_COH
        do n = 1, NPOP
            tp = typInit(n,coh)
            gender = GenderInit(n,coh)
            educ = EducInit(n,coh)
            if (educ <= 2) then
                ed = 0
            elseif (educ >= 3) then
                ed = 1
            end if
            do s = 1, NSECTORS
                if (ChoiceInit(n,coh,LAST_YR-coh) == s) then
                    AllocEmp(s,tp,gender,educ) = AllocEmp(s,tp,gender,educ) + CohortWgt(coh,ed)
                end if
            end do
        end do
    end do  

    do educ = 1, CAT_EDUC
        do gender = 1, CAT_GENDER
            do tp = 1, NTYPES
                AllocEmp(:,tp,gender,educ) = AllocEmp(:,tp,gender,educ) / sum(AllocEmp(:,tp,gender,educ))
            end do
        end do
    end do

    SimExpansion = int(NPOP_SIM / NPOP)

    do coh = FIRST_COH_SIM, LAST_COH_SIM
        if (coh <= LAST_COH) then     
            do ed = 0, 1       
                CohortWgt(coh,ed) = (real(ExpansionFactorGlobal) / real(SimExpansion)) * &
                real(CohortSizeData(coh,ed)) / real(CohortSizeData(FIRST_COH,0))
            end do
        else
            do ed = 0, 1
                CohortWgt(coh,ed) = CohortWgt(LAST_COH,ed)
            end do
        end if
    end do 


    do coh = (LAST_YR-LAST_AGE), (LAST_YR-FIRST_AGE)
        do n = 1, NPOP_SIM
            n2               = int(((n-1)/SimExpansion) + 1)
            EducSim(n,coh)   = EducInit(n2,coh)
            GenderSim(n,coh) = GenderInit(n2,coh)
            ChoiceSim(n,coh,t-coh-9:t-coh) = ChoiceInit(n2,coh,t-coh-9:t-coh)
            typSim(n,coh) = typInit(n2,coh)
        end do
    end do



    ! **********************************************************
    ! From FIRST_YR_SIM+1 (2006) on, all cohorts look the same
    ! as the 1980 cohort in terms of education and gender
    ! but start with zero experience
    ! **********************************************************

    do coh = (FIRST_YR_SIM + 1 - FIRST_AGE), LAST_COH_SIM
        do n = 1, NPOP_SIM
            n2 = int(((n-1)/SimExpansion) + 1)
            EducSim(n,coh)   = EducInit(n2,1980)
            GenderSim(n,coh) = GenderInit(n2,1980)
            ChoiceSim(n,coh,FIRST_AGE-9:FIRST_AGE) = 0
            typSim(n,coh) = typInit(n2,1980)
        end do
    end do 


    ! *****************************************************************************
    ! Read equilibrium returns to skill, productivity (z's), Simulated Capital, etc
    ! that were generated in the estimation code at the optimal parameter vector
    ! *****************************************************************************

    open(unit = 1, file = 'Outcomes.csv')
    read(1,*)
    do 
        read(1,1212,IOSTAT=status2) t, z(1,t), z(2,t), z(3,t), z(4,t), z(5,t), z(6,t), z(7,t), &
                                    rsk_eq(1,0,t), rsk_eq(2,0,t), rsk_eq(3,0,t), rsk_eq(4,0,t), rsk_eq(5,0,t), rsk_eq(6,0,t), rsk_eq(7,0,t), &
                                    rsk_eq(1,1,t), rsk_eq(2,1,t), rsk_eq(3,1,t), rsk_eq(4,1,t), rsk_eq(5,1,t), rsk_eq(6,1,t), rsk_eq(7,1,t), &
                                    sk_sup_eq(1,0,t), sk_sup_eq(2,0,t), sk_sup_eq(3,0,t), sk_sup_eq(4,0,t), sk_sup_eq(5,0,t), sk_sup_eq(6,0,t), sk_sup_eq(7,0,t), & 
                                    sk_sup_eq(1,1,t), sk_sup_eq(2,1,t), sk_sup_eq(3,1,t), sk_sup_eq(4,1,t), sk_sup_eq(5,1,t), sk_sup_eq(6,1,t), sk_sup_eq(7,1,t), & 
                                    CapitalSim(1,t), CapitalSim(2,t), CapitalSim(3,t), CapitalSim(4,t), CapitalSim(5,t), CapitalSim(6,t), CapitalSim(7,t), CapitalData(t)
        if (status2 /= 0 ) exit
    end do
    close(1)

    1212 format(i5,',',43(f20.8,','))


    ! ***************************************************************
    ! This is the last year before the shock in the price of sector 1
    ! ***************************************************************

    allocate(WelfareSim(NPOP_SIM,(HALF_YR_SIM-40-LAST_AGE):(HALF_YR_SIM-FIRST_AGE),FIRST_AGE:LAST_AGE))
    allocate(WageSim(NPOP_SIM,(HALF_YR_SIM-40-LAST_AGE):(HALF_YR_SIM-FIRST_AGE),FIRST_AGE:LAST_AGE))
    allocate(ExperSim(NPOP_SIM,(HALF_YR_SIM-40-LAST_AGE):(HALF_YR_SIM-FIRST_AGE),FIRST_AGE:LAST_AGE)) 

    ! **************************************************
    ! prices of sectors 1 to 7 are fixed until the shock
    ! **************************************************

    price(1:7,FIRST_YR_SIM:LAST_YR_SIM) = 1.0
    do t = FIRST_YR_SIM, LAST_YR_SIM
        prod(:,t)         = z(:,LAST_YR)
        alpha_prod(:,:,t) = alpha_prod(:,:,LAST_YR)
    end do    



    ! ***********************************************
    ! price of sector 3 gets a once and for all shock
    ! ***********************************************
    price(1,HALF_YR_SIM+1:LAST_YR_SIM) = 1.0
    price(2,HALF_YR_SIM+1:LAST_YR_SIM) = 1.0
    price(3,HALF_YR_SIM+1:LAST_YR_SIM) = 0.7 

end if    
    
! ******************************************************
! Broadcast to other processors the vector of parameters
! ******************************************************

Call MPI_BCAST(param,NPARAM,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

count1 = size(typSim)
Call MPI_BCAST(typSim,count1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
count1 = size(AllocEmp)
Call MPI_BCAST(AllocEmp,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
count1 = size(rsk_eq)
Call MPI_BCAST(rsk_eq,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

! *************************************
! Compute Emax_Myopic in each processor
! *************************************

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

! ****************************************************
! Upper and Lower Bounds for draws of returns to skill
! Used in the calculation of Emax
! These are global variables
! ****************************************************

LOW_RSK(1,0) = 0.8
UP_RSK(1,0)  = 2.9

LOW_RSK(2,0) = 0.8
UP_RSK(2,0)  = 3.0

LOW_RSK(3,0) = 1.0
UP_RSK(3,0)  = 3.8

LOW_RSK(4,0) = 1.0
UP_RSK(4,0)  = 3.8

LOW_RSK(5,0) = 0.8
UP_RSK(5,0)  = 3.0

LOW_RSK(6,0) = 1.0
UP_RSK(6,0)  = 3.8

LOW_RSK(7,0) = 0.8
UP_RSK(7,0)  = 3.0



LOW_RSK(1,1) = 1.4
UP_RSK(1,1)  = 5.1

LOW_RSK(2,1) = 1.6
UP_RSK(2,1)  = 5.6

LOW_RSK(3,1) = 2.0
UP_RSK(3,1)  = 6.8

LOW_RSK(4,1) = 2.0
UP_RSK(4,1)  = 6.8

LOW_RSK(5,1) = 1.8
UP_RSK(5,1)  = 5.6

LOW_RSK(6,1) = 2.0
UP_RSK(6,1)  = 6.8

LOW_RSK(7,1) = 1.8
UP_RSK(7,1)  = 5.6


cut = LOW_RSK + (UP_RSK-LOW_RSK)/2.0

! ******************************
! Compute Emax in each processor
! ******************************

if (rank == 0) then
    print*
    print*, '************************'
    print*, 'Computing Emax Myopic...'
    print*, '************************'
    print*
end if

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

count1 = (LAST_AGE - FIRST_AGE)*(NSECTORS+1)*NREG
count2 = (LAST_AGE - FIRST_AGE)*(NSECTORS+1)


! *****************
! Division of Labor
! *****************

tp = int(rank/(CAT_EDUC*CAT_GENDER)) + 1
if (mod(rank,CAT_EDUC*CAT_GENDER) <= CAT_EDUC - 1) then
    gender = 1
else
    gender = 2
end if
educ = mod(mod(rank,CAT_EDUC*CAT_GENDER),CAT_EDUC) + 1

Call EmaxCoef_Myopic(param,cut,gender,educ,tp,A,B)

! **********************************************************
! Slave processors send Emax results to the Master processor
! **********************************************************

if (rank /= 0) then
    tag = 1
    Call MPI_SEND(A,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierr)
    tag = 2
    Call MPI_SEND(B,count2,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierr)
end if


! ************************************************************
! The Master processor receives the messages and stock results
! into the PI and R_sq arrays
! ************************************************************

if (rank == 0) then
    do age = FIRST_AGE+1, LAST_AGE
        do lag = 0, NSECTORS
            PI_Myopic(age,lag,gender,educ,tp,:) = A(age,lag,:)
            R_sq_Myopic(age,lag,gender,educ,tp) = B(age,lag)
        end do
    end do
    do proc = 1, numproc-1
        tp = int(proc/(CAT_EDUC*CAT_GENDER)) + 1
        if (mod(proc,CAT_EDUC*CAT_GENDER) <= CAT_EDUC - 1) then
            gender = 1
        else
            gender = 2
        end if
        educ = mod(mod(proc,CAT_EDUC*CAT_GENDER),CAT_EDUC) + 1
        tag = 1
        Call MPI_RECV(A,count1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
        tag = 2
        Call MPI_RECV(B,count2,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
        do age = FIRST_AGE+1, LAST_AGE
            do lag = 0, NSECTORS
                PI_Myopic(age,lag,gender,educ,tp,:) = A(age,lag,:)
                R_sq_Myopic(age,lag,gender,educ,tp) = B(age,lag)
            end do
        end do
    end do
end if

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

if (rank == 0) then
    print*
    print*, '***********************************'
    print*, 'Emax Myopic Computed'
    print*, '***********************************'
    print*
    print*, 'Min and Max R_sq for Emax:'
    print*, minval(R_sq_Myopic(:,1:NSECTORS,:,:,:)), minloc(R_sq_Myopic(:,1:NSECTORS,:,:,:))
    print*, maxval(R_sq_Myopic(:,1:NSECTORS,:,:,:)), maxloc(R_sq_Myopic(:,1:NSECTORS,:,:,:))  
    
end if            

count1 = size(PI_Myopic)
Call MPI_BCAST(PI_Myopic,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

PI_RE = 0.0

RE_loop = 1

avec = 1.0

RatExp_loop: do while (RE_loop <= 2)
    
    if (RE_loop >= 2) then

        if (rank == 0) then
            print*
            print*, '***************************************'
            print*, 'Computing Emax Rational Expectations...'
            print*, '***************************************'
            print*
        end if

        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

        ! *********************************
        ! Compute Emax_RE in each processor
        ! *********************************
        count1 = size(A_RE)
        count2 = size(B_RE)

        ! *****************
        ! Division of Labor
        ! *****************

        tp = int(rank/(CAT_EDUC*CAT_GENDER)) + 1
        if (mod(rank,CAT_EDUC*CAT_GENDER) <= CAT_EDUC - 1) then
            gender = 1
        else
            gender = 2
        end if
        educ = mod(mod(rank,CAT_EDUC*CAT_GENDER),CAT_EDUC) + 1


        Call EmaxCoef_RE(param,avec,cut,gender,educ,tp,A_RE,B_RE,PI_Myopic)


        ! **********************************************************
        ! Slave processors send Emax results to the Master processor
        ! **********************************************************

        if (rank /= 0) then
            tag = 1
            Call MPI_SEND(A_RE,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierr)
            tag = 2
            Call MPI_SEND(B_RE,count2,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierr)
        end if

        ! ************************************************************
        ! The Master processor receives the messages and stock results
        ! into the PI and R_sq arrays
        ! ************************************************************


        if (rank == 0) then
            do t = (YR_ANNMNT+1), (LAST_YR_SIM-1)
                do age = FIRST_AGE+1, LAST_AGE
                    do lag = 0, NSECTORS
                        PI_RE(t,age,lag,gender,educ,tp,:) = A_RE(t,age,lag,:)
                        R_sq_RE(t,age,lag,gender,educ,tp) = B_RE(t,age,lag)
                    end do
                end do
            end do
            do proc = 1, numproc-1
                tp = int(proc/(CAT_EDUC*CAT_GENDER)) + 1
                if (mod(proc,CAT_EDUC*CAT_GENDER) <= CAT_EDUC - 1) then
                    gender = 1
                else
                    gender = 2
                end if
                educ = mod(mod(proc,CAT_EDUC*CAT_GENDER),CAT_EDUC) + 1
                tag = 1
                Call MPI_RECV(A_RE,count1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
                tag = 2
                Call MPI_RECV(B_RE,count2,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
                do t = (YR_ANNMNT+1), (LAST_YR_SIM-1)
                    do age = FIRST_AGE+1, LAST_AGE
                        do lag = 0, NSECTORS
                            PI_RE(t,age,lag,gender,educ,tp,:) = A_RE(t,age,lag,:)
                            R_sq_RE(t,age,lag,gender,educ,tp) = B_RE(t,age,lag)
                        end do
                    end do
                end do
            end do
        end if

        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if (rank == 0) then
            print*
            print*, '***********************************'
            print*, 'Emax Rational Expectations Computed'
            print*, '***********************************'
            print*
            print*, 'Min and Max R_sq for Emax:'
            print*, minval(R_sq_RE(:,:,1:NSECTORS,:,:,:)), minloc(R_sq_RE(:,:,1:NSECTORS,:,:,:))
            print*, maxval(R_sq_RE(:,:,1:NSECTORS,:,:,:)), maxloc(R_sq_RE(:,:,1:NSECTORS,:,:,:))
        
        end if
        
        count1 = size(PI_RE)
        Call MPI_BCAST(PI_RE,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    end if

    if (rank == 0) then
        
        !*********************
        ! For the random draws
        !*********************

        ! Normal Draws

        brng         = VSL_BRNG_MT19937
        method_gauss = VSL_METHOD_DGAUSSIAN_BOXMULLER
        seed         = 112233445

        ! ***** Initializing *****
        errcode = vslnewstream(stream_gauss, brng, seed)
        ! ***** Warming Up *****
        do j = 1, WARMUP
            errcode = vdrnggaussian(method_gauss, stream_gauss, NSECTORS+1, &
                                    rvec_eps, real(0.0,DOUBLE), real(1.0,DOUBLE))
        end do 


        ! Gumbel Draws

        brng          = VSL_BRNG_MT19937
        method_gumbel = VSL_METHOD_DGUMBEL_ICDF
        seed          = 223344556

        ! ***** Initializing *****
        errcode = vslnewstream(stream_gumbel, brng, seed)
        loc   = 0.0
        scale = 1.0
    
        ! ***** Warming Up *****
        do j = 1, WARMUP
            errcode = vdrnggumbel(method_gumbel, stream_gumbel, NSECTORS+1, rvec_eta, loc, scale)
        end do
  
        ! ************
        ! Initializing
        ! ************

        p_sup(4:7,1) = 1.0
        p_dem(4:7,1) = 0.0

        rK_eq(FIRST_YR_SIM:LAST_YR_SIM)    = 0.0

        emp_eq    = 0    
        iter      = 0
        step      = 0.0
        step2     = 0.0

        ! *******************************************************
        ! The Total Capital Stock is throughout fixed at the 2005
        ! simulated value
        ! *******************************************************

        CapitalStock = sum(CapitalSim(:,LAST_YR))
        
    end if
    
    sk_sup_eq = 0.0
    rsk_eq(:,:,FIRST_YR_SIM:LAST_YR_SIM) = 0.0 
    rsk_sup        = 0.0
    rsk_sup(:,0,1) = rsk_eq(:,0,LAST_YR-1)
    rsk_sup(:,1,1) = rsk_eq(:,1,LAST_YR-1)

    year_loop: do t = FIRST_YR_SIM, LAST_YR_SIM

        year = t
        if (rank == 0) then
            print*, year
        end if
        
        if (rank == 0) then
    
            ! ***** Random Draws *****
            do gen = 1, 36
                do n = 1, NPOP_SIM
                    errcode = vdrnggaussian(method_gauss, stream_gauss, NSECTORS+1, &
                              rvec_eps, real(0.0,DOUBLE), real(1.0,DOUBLE))
                    do s = 0, NSECTORS
                        eps(n,gen,s) = sigma(s)*rvec_eps(s)
                    end do            
                    errcode = vdrnggumbel(method_gauss, stream_gumbel, NSECTORS+1, rvec_eta, loc, scale)
                    do s = 0, NSECTORS
                        eta(n,gen,s) = -sigmaPref*rvec_eta(s) - sigmaPref*0.577215665
                    end do
                end do
            end do  
        
        end if
      
    
        if (t > FIRST_YR_SIM) then
                
            ! Use the result of the previous year
            ! As the initial point for rsk
            rsk_sup        = 0.0
            rsk_sup(:,:,1) = rsk_eq(:,:,t-1)
            if (rank == 0) then
                rsk_dem        = 0.0
                rK_sup(1)      = rK_eq(t-1)
                rK_dem         = 0.0
                p_sup          = 0.0
                p_dem          = 0.0
                p_sup(4:7,1)   = price(4:7,t-1)
            end if
                
        end if


        it        = 1
        flag      = 0
        converged = .false.
    
    
        while_loop: do while (.not. converged .and. it <= MAXIT)
                
            sk_sup               = 0.0
            emp(:,it)            = 0.0
            emp_by_skill(:,:,t)  = 0.0
            Welfare(t)           = 0.0
            Welfare_Educ(:,t)    = 0.0
            Real_Wages(t)        = 0.0
            Real_Wages_Educ(:,t) = 0.0
            
            ! Compute aggregate supply of skills  
        
            do coh_batch = 2, NCOH_BATCH
                
                if (rank == 0) then
                    
                    older_coh = (year-LAST_AGE) + (coh_batch-1)*COH_PROC
                    newer_coh = (year-LAST_AGE) + coh_batch*COH_PROC - 1
                    older_gen = older_coh - (year-LAST_AGE) + 1
                    newer_gen = newer_coh - (year-LAST_AGE) + 1
                        
                    do coh = older_coh, newer_coh
                        age = t - coh
                        ind = coh - older_coh + 1
                        ChoiceSimAux(:,ind,-9:-1) = ChoiceSim(:,coh,age-9:age-1)
                    end do
                
                    EducSimAux = EducSim(:,older_coh:newer_coh)
                    GenderSimAux = GenderSim(:,older_coh:newer_coh)                
                    CohortWgtAux = CohortWgt(older_coh:newer_coh,:)
                
                    epsAux = eps(:,older_gen:newer_gen,:)
                    etaAux = eta(:,older_gen:newer_gen,:)

                    tag = coh_batch*100 + 1                        
                    count1 = size(ChoiceSimAux) 
                    Call MPI_SEND(ChoiceSimAux,count1,MPI_INTEGER,coh_batch-1,tag,MPI_COMM_WORLD,status,ierr)
                        
                    tag = coh_batch*100 + 2                      
                    count1 = size(EducSimAux) 
                    Call MPI_SEND(EducSimAux,count1,MPI_INTEGER,coh_batch-1,tag,MPI_COMM_WORLD,status,ierr)
                        
                    tag = coh_batch*100 + 3                    
                    count1 = size(GenderSimAux) 
                    Call MPI_SEND(GenderSimAux,count1,MPI_INTEGER,coh_batch-1,tag,MPI_COMM_WORLD,status,ierr)
                    
                    tag = coh_batch*100 + 4                    
                    count1 = size(CohortWgtAux) 
                    Call MPI_SEND(CohortWgtAux,count1,MPI_DOUBLE_PRECISION,coh_batch-1,tag,MPI_COMM_WORLD,status,ierr)

                    tag = coh_batch*100 + 5
                    count1 = size(epsAux)
                    Call MPI_SEND(epsAux,count1,MPI_DOUBLE_PRECISION,coh_batch-1,tag,MPI_COMM_WORLD,status,ierr)

                    tag = coh_batch*100 + 6
                    count1 = size(etaAux)
                    Call MPI_SEND(etaAux,count1,MPI_DOUBLE_PRECISION,coh_batch-1,tag,MPI_COMM_WORLD,status,ierr)

                end if

                if (rank == coh_batch-1) then

                    tag = coh_batch*100 + 1
                    count1 = size(ChoiceSimAux)
                    Call MPI_RECV(ChoiceSimAux,count1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
                        
                    tag = coh_batch*100 + 2                      
                    count1 = size(EducSimAux) 
                    Call MPI_RECV(EducSimAux,count1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
                        
                    tag = coh_batch*100 + 3                    
                    count1 = size(GenderSimAux) 
                    Call MPI_RECV(GenderSimAux,count1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
                        
                    tag = coh_batch*100 + 4                    
                    count1 = size(CohortWgtAux) 
                    Call MPI_RECV(CohortWgtAux,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)

                    tag = coh_batch*100 + 5
                    count1 = size(epsAux)
                    Call MPI_RECV(epsAux,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)

                    tag = coh_batch*100 + 6
                    count1 = size(etaAux)
                    Call MPI_RECV(etaAux,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)

                end if

            end do


            if (rank == 0) then
                    
                coh_batch = 1
                older_coh = (year-LAST_AGE) + (coh_batch-1)*COH_PROC
                newer_coh = (year-LAST_AGE) + coh_batch*COH_PROC - 1
                older_gen = older_coh - (year-LAST_AGE) + 1
                newer_gen = newer_coh - (year-LAST_AGE) + 1
                    
                do coh = older_coh, newer_coh
                    age = t - coh
                    ind = coh - older_coh + 1
                    ChoiceSimAux(:,ind,-9:-1) = ChoiceSim(:,coh,age-9:age-1)
                end do
                
                EducSimAux = EducSim(:,older_coh:newer_coh)
                GenderSimAux = GenderSim(:,older_coh:newer_coh)                
                CohortWgtAux = CohortWgt(older_coh:newer_coh,:)
                
                epsAux = eps(:,older_gen:newer_gen,:)
                etaAux = eta(:,older_gen:newer_gen,:)
                
            end if
                
            older_coh = (year-LAST_AGE) + rank*COH_PROC
            newer_coh = (year-LAST_AGE) + (rank+1)*COH_PROC - 1

            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            
            if (rank <= NCOH_BATCH-1) then
                Call ParallelCohorts(RE_loop, older_coh, newer_coh, t, param, cut, avec, EducSimAux, GenderSimAux, typSim, AllocEmp, rsk_sup(:,:,it), &
                epsAux, etaAux, PI_Myopic, PI_RE(t+1,:,:,:,:,:,:), &
                ChoiceSimAux, CohortWgtAux, ChoiceSimOut, &
                ExperSimOut, WageSimOut, WelfareSimOut, &
                empOut, emp_by_skillOut, sk_supOut, WelfareOut, Welfare_EducOut, Real_WagesOut, Real_Wages_EducOut)
            end if
                
            if (rank /= 0 .and. rank <= NCOH_BATCH-1) then
                
                tag = (rank+1)*100 + 1
                count1 = size(ChoiceSimOut)
                Call MPI_SEND(ChoiceSimOut,count1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
                    
                tag = (rank+1)*100 + 2
                count1 = size(ExperSimOut)
                Call MPI_SEND(ExperSimOut,count1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
                    
                tag = (rank+1)*100 + 3
                count1 = size(WageSimOut)
                Call MPI_SEND(WageSimOut,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)
                    
                tag = (rank+1)*100 + 4
                count1 = size(WelfareSimOut)
                Call MPI_SEND(WelfareSimOut,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)
                    
                tag = (rank+1)*100 + 5
                count1 = size(empOut)
                Call MPI_SEND(empOut,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)
                    
                tag = (rank+1)*100 + 6
                count1 = size(emp_by_skillOut)
                Call MPI_SEND(emp_by_skillOut,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)                    
                    
                tag = (rank+1)*100 + 7
                count1 = size(sk_supOut)
                Call MPI_SEND(sk_supOut,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)
                    
                tag = (rank+1)*100 + 8
                count1 = 1
                Call MPI_SEND(WelfareOut,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)
                    
                tag = (rank+1)*100 + 9
                count1 = size(Welfare_EducOut)
                Call MPI_SEND(Welfare_EducOut,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)
                    
                tag = (rank+1)*100 + 10
                count1 = 1
                Call MPI_SEND(Real_WagesOut,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)
                    
                tag = (rank+1)*100 + 11
                count1 = size(Real_Wages_EducOut)
                Call MPI_SEND(Real_Wages_EducOut,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)
                
            end if

            if (rank == 0) then
                    
                older_coh = (year-LAST_AGE) + rank*COH_PROC
                newer_coh = (year-LAST_AGE) + (rank+1)*COH_PROC - 1
                    
                do coh = older_coh, newer_coh
                    age = t - coh
                    ind = coh - older_coh + 1
                    ChoiceSim(:,coh,age) = ChoiceSimOut(:,ind)
                end do
                
                do coh = older_coh, newer_coh
                    if (coh >= HALF_YR_SIM-40-LAST_AGE .and. coh <= HALF_YR_SIM-FIRST_AGE) then
                        age = t - coh
                        ind = coh - older_coh + 1
                        ExperSim(:,coh,age)   = ExperSimOut(:,ind)
                        WelfareSim(:,coh,age) = WelfareSimOut(:,ind)
                        WageSim(:,coh,age) = WageSimOut(:,ind)
                    end if
                end do
                    
                sk_sup               = sk_sup + sk_supOut
                emp(:,it)            = emp(:,it) + empOut
                emp_by_skill(:,:,t)  = emp_by_skill(:,:,t) + emp_by_skillOut
                Welfare(t)           = Welfare(t) + WelfareOut
                Welfare_Educ(:,t)    = Welfare_Educ(:,t) + Welfare_EducOut
                Real_Wages(t)        = Real_Wages(t) + Real_WagesOut
                Real_Wages_Educ(:,t) = Real_Wages_Educ(:,t) + Real_Wages_EducOut
                    
                do proc = 1, NCOH_BATCH-1
                    
                    older_coh = (year-LAST_AGE) + proc*COH_PROC
                    newer_coh = (year-LAST_AGE) + (proc+1)*COH_PROC - 1

                    tag = (proc+1)*100 + 1
                    count1 = size(ChoiceSimOut)
                    Call MPI_RECV(ChoiceSimOut,count1,MPI_INTEGER,proc,tag,MPI_COMM_WORLD,status,ierr)
                    do coh = older_coh, newer_coh
                        age = t - coh
                        ind = coh - older_coh + 1
                        ChoiceSim(:,coh,age) = ChoiceSimOut(:,ind)
                    end do

                    tag = (proc+1)*100 + 2
                    count1 = size(ExperSimOut)
                    Call MPI_RECV(ExperSimOut,count1,MPI_INTEGER,proc,tag,MPI_COMM_WORLD,status,ierr)
                    do coh = older_coh, newer_coh
                        if (coh >= HALF_YR_SIM-40-LAST_AGE .and. coh <= HALF_YR_SIM-FIRST_AGE) then
                            age = t - coh
                            ind = coh - older_coh + 1
                            ExperSim(:,coh,age)   = ExperSimOut(:,ind)
                        end if
                    end do                      
                       
                    tag = (proc+1)*100 + 3
                    count1 = size(WageSimOut)
                    Call MPI_RECV(WageSimOut,count1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
                    do coh = older_coh, newer_coh
                        if (coh >= HALF_YR_SIM-40-LAST_AGE .and. coh <= HALF_YR_SIM-FIRST_AGE) then
                            age = t - coh
                            ind = coh - older_coh + 1
                            WageSim(:,coh,age)   = WageSimOut(:,ind)
                        end if
                    end do
                        
                    tag = (proc+1)*100 + 4
                    count1 = size(WelfareSimOut)
                    Call MPI_RECV(WelfareSimOut,count1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
                    do coh = older_coh, newer_coh
                        if (coh >= HALF_YR_SIM-40-LAST_AGE .and. coh <= HALF_YR_SIM-FIRST_AGE) then
                            age = t - coh
                            ind = coh - older_coh + 1
                            WelfareSim(:,coh,age)   = WelfareSimOut(:,ind)
                        end if
                    end do
                        
                    tag = (proc+1)*100 + 5
                    count1 = size(empOut)
                    Call MPI_RECV(empOut,count1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
                    emp(:,it) = emp(:,it) + empOut
                    
                    tag = (proc+1)*100 + 6
                    count1 = size(emp_by_skillOut)
                    Call MPI_RECV(emp_by_skillOut,count1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
                    emp_by_skill(:,:,t) = emp_by_skill(:,:,t) + emp_by_skillOut                        
                        
                    tag = (proc+1)*100 + 7
                    count1 = size(sk_supOut)
                    Call MPI_RECV(sk_supOut,count1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
                    sk_sup = sk_sup + sk_supOut
                        
                    tag = (proc+1)*100 + 8
                    count1 = 1
                    Call MPI_RECV(WelfareOut,count1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
                    Welfare(t) = Welfare(t) + WelfareOut
                        
                    tag = (proc+1)*100 + 9
                    count1 = size(Welfare_EducOut)
                    Call MPI_RECV(Welfare_EducOut,count1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
                    Welfare_Educ(:,t) = Welfare_Educ(:,t) + Welfare_EducOut
                    
                    tag = (proc+1)*100 + 10
                    count1 = 1
                    Call MPI_RECV(Real_WagesOut,count1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
                    Real_Wages(t) = Real_Wages(t) + Real_WagesOut
                        
                    tag = (proc+1)*100 + 11
                    count1 = size(Real_Wages_EducOut)
                    Call MPI_RECV(Real_Wages_EducOut,count1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
                    Real_Wages_Educ(:,t) = Real_Wages_Educ(:,t) + Real_Wages_EducOut
                    
                end do

            end if
            
            if (rank == 0) then
        
                Pindex(t) = (price(1,t)**alpha_c(1))  * (price(2,t)**alpha_c(2))  * & 
                            (price(3,t)**alpha_c(3))  * (p_sup(4,it)**alpha_c(4)) * &
                            (p_sup(5,it)**alpha_c(5)) * (p_sup(6,it)**alpha_c(6)) * &
                            (p_sup(7,it)**alpha_c(7)) 

                z(1:3,t) = price(1:3,t)*prod(1:3,t)*Pindex(LAST_YR)/Pindex(t)
                z(4:7,t) = p_sup(4:7,it)*prod(4:7,t)*Pindex(LAST_YR)/Pindex(t)
        
        
                ! *******************************************************************
                ! Compute the allocation of capital and the rate of return of capital
                ! *******************************************************************
        
                ! I allow total capital in the economy to grow
                ! With the total skill supply until the economy attains its steady state
                ! For that, returns to capital are fixed until the SS is attained
        
                if (t < HALF_YR_SIM - 10) then
                    do s = 1, NSECTORS
                        Const1 = alpha_prod(s,0,LAST_YR)*sk_sup(s,0)**sigma_prod(s) + alpha_prod(s,1,LAST_YR)*sk_sup(s,1)**sigma_prod(s)
                        Const2 = rKData(LAST_YR)/z(s,t)
                        Const3 = 1 - alpha_prod(s,0,LAST_YR) - alpha_prod(s,1,LAST_YR)
                        if (t == FIRST_YR_SIM) then
                            init = CapitalSim(s,LAST_YR)
                        else
                            init = CapitalSim(s,t-1)
                        end if
                        Call NewtonSolver(Const1,Const2,Const3,sigma_prod(s),init,solution)
                        CapitalSim(s,t) = solution
                    end do
                    CapitalStock = sum(CapitalSim(:,t))
                end if
            
                ! *******************************************************************
                ! rK equilibrium and Capital that results from rsk_sup and from p_sup
                ! *******************************************************************
            
                if (t == FIRST_YR_SIM) then
                    !rK_sup     = 0.0
                    !rK_dem     = 0.0
                    rK_sup(1)  = rKData(LAST_YR)
                else
                    !rK_sup     = 0.0
                    !rK_dem     = 0.0
                    rK_sup(1)  = rK_eq(t-1)
                    if (it > 1) then
                        rK_sup(1) = rK_sup(it2-1)
                    end if
                end if
        
                it2        = 1
                flag2      = 0
                converged2 = .false.
        
                do while (.not. converged2 .and. it2 <= MAXIT)
                
                    do s = 1, NSECTORS-1
                        Const1 = alpha_prod(s,0,LAST_YR)*sk_sup(s,0)**sigma_prod(s) + alpha_prod(s,1,LAST_YR)*sk_sup(s,1)**sigma_prod(s)
                        Const2 = rK_sup(it2)/z(s,t)
                        Const3 = 1 - alpha_prod(s,0,LAST_YR) - alpha_prod(s,1,LAST_YR)
                        if (t == FIRST_YR_SIM) then
                            init = CapitalSim(s,LAST_YR)
                        else
                            init = CapitalSim(s,t-1)
                            if (it > 1) then
                                init = Cap(s)
                            end if
                        end if
                        Call NewtonSolver(Const1,Const2,Const3,sigma_prod(s),init,solution)
                        Cap(s) = solution
                    end do

                    Cap(7) = max(CapitalStock - sum(Cap(1:6)),0.00001)

                    rK_dem(it2) = (1 - alpha_prod(7,0,LAST_YR) - alpha_prod(7,1,LAST_YR))*z(7,t)* &
                                  ((alpha_prod(7,0,LAST_YR)*sk_sup(7,0)**sigma_prod(7) + alpha_prod(7,1,LAST_YR)*sk_sup(7,1)**sigma_prod(7) + &
                                   (1 - alpha_prod(7,0,LAST_YR) - alpha_prod(7,1,LAST_YR))*Cap(7)**sigma_prod(7))**(1/sigma_prod(7) - 1))* &
                                   Cap(7)**(sigma_prod(7)-1)
           
                    check2 = abs(log(rK_sup(it2)) - log(rK_dem(it2)))
            
                    converged2 = (check2 <= crit2)
            
                    if ( .not. converged2 .and. it2 < MAXIT ) then
                
                        if (it2 < 30) then
                            step2 = 1.0/10.0
                        else 
                            if (SIGN(1.0, rK_dem(it2)-rK_dem(it2-1)) .ne.  & 
                                SIGN (1.0, rK_dem(it2-1)-rK_dem(it2-2)) .and. flag2 == 0) then
                                    flag2 = it2
                                    step2 = step2 / 1.3
                            end if
                            if (SIGN(1.0, rK_dem(it2)-rK_dem(it2-1)) .ne.  & 
                                SIGN (1.0, rK_dem(it2-1)-rK_dem(it2-2)) .and. it2 >= flag2+5) then
                                    flag2 = it2
                                    step2 = step2 / 1.3                     
                            end if
                        end if
                
                        ! Updating the return to capital
                        rK_sup(it2+1) = rK_sup(it2) - step2*(rK_sup(it2) - rK_dem(it2))
                
                    else
            
                        rK_eq(t) = rK_dem(it2)
                        iter2(t) = it2
                
                    end if            
            
                    it2 = it2 + 1
        
                end do

        
                do s = 1, NSECTORS
                    output(s,t) = z(s,t)*(alpha_prod(s,0,LAST_YR)*sk_sup(s,0)**sigma_prod(s) + alpha_prod(s,1,LAST_YR)*sk_sup(s,1)**sigma_prod(s) + &
                                          (1 - alpha_prod(s,0,LAST_YR) - alpha_prod(s,1,LAST_YR))*Cap(s)**sigma_prod(s))**(1/sigma_prod(s))
                end do
            
                ! ***************************************************************************************
                ! Given rsk_sup, rK_eq(rsk_sup) and Capital(rsk_sup) we compute the equilibrium NT prices
                ! ***************************************************************************************
            
                do s = 4, 7
                    p_dem(s,it) = (alpha_c(s)/(1-alpha_c(s))) * (sum(output(1:7,t))-output(s,t)) / &
                                   (prod(s,t)*(alpha_prod(s,0,LAST_YR)*sk_sup(s,0)**sigma_prod(s) + alpha_prod(s,1,LAST_YR)*sk_sup(s,1)**sigma_prod(s) + &
                                   (1 - alpha_prod(s,0,LAST_YR) - alpha_prod(s,1,LAST_YR))*Cap(s)**sigma_prod(s))**(1/sigma_prod(s)))
                    p_dem(s,it) = p_dem(s,it)*(Pindex(t)/Pindex(LAST_YR))
                end do                    
                
                do s = 1, NSECTORS
                    do ed = 0, 1
                        if (sk_sup(s,ed) > 0) then
                            rsk_dem(s,ed,it) = min(( alpha_prod(s,ed,t) / (1 - alpha_prod(s,0,t) - alpha_prod(s,1,t)) )*rK_eq(t)*((Cap(s)/sk_sup(s,ed))**(1-sigma_prod(s))),1000.0)
                        else
                            rsk_dem(s,ed,it) = min(2*rsk_sup(s,ed,it),1000.0)
                        end if
                    end do
                end do
    
                check1(:,:) = abs(log(rsk_sup(:,:,it)) - log(rsk_dem(:,:,it)))
                check3      = abs(log(p_sup(4:7,it)) - log(p_dem(4:7,it)))
        
        
                converged = (check1(1,0) <= crit1 .and. check1(2,0) <= crit1 .and. check1(3,0) <= crit1 .and. &
                             check1(4,0) <= crit1 .and. check1(1,1) <= crit1 .and. check1(2,1) <= crit1 .and. &
                             check1(3,1) <= crit1 .and. check1(4,1) <= crit1 .and. check3(4)   <= crit3 .and. &
                             check3(5)   <= crit3 .and. check3(6)   <= crit3 .and. check3(7)   <= crit3)            
           
             
                if ( .not. converged .and. it < MAXIT ) then
               
                    ! Adjust the step for updating the skill price
                    ! If rsk oscillates, decrease the step            
                    if (it < 30) then
                        step  = 1.0/6.0
                        step3 = 1.0/6.0
                    else 
                        do s = 1, NSECTORS
                            do ed = 0, 1
                                if (SIGN(1.0, rsk_dem(s,ed,it)-rsk_dem(s,ed,it-1)) .ne.  & 
                                    SIGN (1.0, rsk_dem(s,ed,it-1)-rsk_dem(s,ed,it-2)) .and. flag(s,ed) == 0) then
                                    flag(s,ed) = it
                                    step(s,ed) = step(s,ed) / 1.3                     
                                end if
                                if (SIGN(1.0, rsk_dem(s,ed,it)-rsk_dem(s,ed,it-1)) .ne.  & 
                                    SIGN (1.0, rsk_dem(s,ed,it-1)-rsk_dem(s,ed,it-2)) .and. it >= flag(s,ed)+5) then
                                    flag(s,ed) = it
                                    step(s,ed) = step(s,ed) / 1.3                     
                                end if
                            end do
                        end do
                        do s = 4, 7
                            if (SIGN(1.0, p_dem(s,it)-p_dem(s,it-1)) .ne.  & 
                                SIGN (1.0, p_dem(s,it-1)-p_dem(s,it-2)) .and. flag3 == 0) then
                                flag3 = it
                                step3(s) = step3(s) / 1.3
                            end if
                            if (SIGN(1.0, p_dem(s,it)-p_dem(s,it-1)) .ne.  & 
                                SIGN (1.0, p_dem(s,it-1)-p_dem(s,it-2)) .and. it >= flag3+5) then
                                flag3 = it
                                step3(s) = step3(s) / 1.3
                            end if                            
                        end do
                    end if
           
                    ! Updating the return to skill
                    do s = 1, NSECTORS
                        do ed = 0, 1
                            rsk_sup(s,ed,it+1) = rsk_sup(s,ed,it) - step(s,ed)*(rsk_sup(s,ed,it) - rsk_dem(s,ed,it))    
                        end do
                    end do 
                    
                    do s = 4, 7
                        p_sup(s,it+1) = p_sup(s,it) - step3(s)*(p_sup(s,it) - p_dem(s,it))
                    end do
            
                else 
        
                    ! Recording Equilibrium Quantities
                    rsk_eq(:,:,t)    = rsk_dem(:,:,it)   
                    CapitalSim(:,t)  = Cap
                    sk_sup_eq(:,:,t) = sk_sup        
                    emp_eq(:,t)      = emp(:,it)            
                    iter(t)          = it
                    price(4:7,t)     = p_dem(4:7,it)
                    Welfare(t)       = Welfare(t) + rK_eq(t)*sum(CapitalSim(:,t))
    
                end if
                
            end if

            count1 = size(rsk_sup)
            Call MPI_BCAST(rsk_sup,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            count1 = size(rsk_eq)
            Call MPI_BCAST(rsk_eq,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
                
            it = it + 1        
     
        end do while_loop
    
        if (rank == 0) then
            
            print*
            print*
            print*, 'RE Iteration:', RE_loop
            print*
            print*, 'rK:'
            print*, rK_eq(t)
            print*
            print*, 'rsk Ed=0:'
            print*, rsk_eq(1,0,t), rsk_eq(2,0,t), rsk_eq(3,0,t), rsk_eq(4,0,t), rsk_eq(5,0,t), rsk_eq(6,0,t), rsk_eq(7,0,t)
            print*
            print*, 'rsk Ed=1:'
            print*, rsk_eq(1,1,t), rsk_eq(2,1,t), rsk_eq(3,1,t), rsk_eq(4,1,t), rsk_eq(5,1,t), rsk_eq(6,1,t), rsk_eq(7,1,t)
            print*
            print*, 'Employment:'
            print*, emp_eq(1,t)/real(sum(emp_eq(1:7,t))), emp_eq(2,t)/real(sum(emp_eq(1:7,t))), &
                    emp_eq(3,t)/real(sum(emp_eq(1:7,t))), emp_eq(4,t)/real(sum(emp_eq(1:7,t))), &
                    emp_eq(5,t)/real(sum(emp_eq(1:7,t))), emp_eq(6,t)/real(sum(emp_eq(1:7,t))), &
                    emp_eq(7,t)/real(sum(emp_eq(1:7,t)))
            print*
            print*, 'Fraction ZERO Sector:'
            print*, emp_eq(0,t)/real(sum(emp_eq(0:7,t)))          
            print*
            print*, 'Capital Stock:'
            print*, CapitalStock
            print*
            print*, 'Welfare:'
            print*, Welfare(t)
            print*
            print*, 'Prices NT:'
            print*, price(4,t), price(5,t), price(6,t), price(7,t)
            print*
            print*, 'iter:'
            print*, iter(t)
            print*
            print*, 'iter2:'
            print*, iter2(t)
            print*
            print*, 'check1:'
            print*, check1
            print*
            print*, 'check2:'
            print*, check2
            print*
            print*, 'check3:'
            print*, check3
            print*
            print*
            print*, '================================================================'
            print*, '================================================================'
            
        end if
            
    end do year_loop

    if (rank == 0) then
        
        avec = 1.0

        do t = YR_ANNMNT, LAST_YR_SIM-1
            do ed = 0,1
                do s = 1, NSECTORS
                    avec(s,ed,t+1) = rsk_eq(s,ed,t+1) / rsk_eq(s,ed,t)
                end do
            end do
        end do

        if (RE_loop == 1) then
            open(unit = 1, file = 'RSK_EQ1.csv')
            do ed = 0, 1
                do t = FIRST_YR_SIM, LAST_YR_SIM
                    write(1,3476) ed, t, rsk_eq(1,ed,t), rsk_eq(2,ed,t), rsk_eq(3,ed,t), rsk_eq(4,ed,t), rsk_eq(5,ed,t), rsk_eq(6,ed,t), rsk_eq(7,ed,t), & 
                                  avec(1,ed,t), avec(2,ed,t), avec(3,ed,t), avec(4,ed,t), avec(5,ed,t), avec(6,ed,t), avec(7,ed,t)
                end do
            end do    
            close(1)
        elseif (RE_loop == 2) then
            open(unit = 1, file = 'RSK_EQ2.csv')
            do ed = 0, 1
                do t = FIRST_YR_SIM, LAST_YR_SIM
                    write(1,3476) ed, t, rsk_eq(1,ed,t), rsk_eq(2,ed,t), rsk_eq(3,ed,t), rsk_eq(4,ed,t), rsk_eq(5,ed,t), rsk_eq(6,ed,t), rsk_eq(7,ed,t), & 
                                  avec(1,ed,t), avec(2,ed,t), avec(3,ed,t), avec(4,ed,t), avec(5,ed,t), avec(6,ed,t), avec(7,ed,t)
                end do
            end do
            close(1)
        elseif (RE_loop == 3) then
            open(unit = 1, file = 'RSK_EQ3.csv')
            do ed = 0, 1
                do t = FIRST_YR_SIM, LAST_YR_SIM
                    write(1,3476) ed, t, rsk_eq(1,ed,t), rsk_eq(2,ed,t), rsk_eq(3,ed,t), rsk_eq(4,ed,t), rsk_eq(5,ed,t), rsk_eq(6,ed,t), rsk_eq(7,ed,t), & 
                                  avec(1,ed,t), avec(2,ed,t), avec(3,ed,t), avec(4,ed,t), avec(5,ed,t), avec(6,ed,t), avec(7,ed,t)
                end do
            end do
            close(1)  
        end if    
        
        3476 format(2(i10,','),14(f20.6,','))
             
    end if
    
    ! ********************************
    ! End of processing by processor 0
    ! ********************************

    RE_loop = RE_loop + 1
    
    ! Broadcast avec to the other processors
    count1 = size(avec)
    Call MPI_BCAST(avec,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

end do RatExp_loop
    
             


if (rank == 0) then


    ! ***** Deinitialize *****
    errcode = vsldeletestream(stream_gauss)
    errcode = vsldeletestream(stream_gumbel)


    ! **********
    ! Pre Shock
    ! **********

    Size2 = 0
    do coh = HALF_YR_SIM - 40 - LAST_AGE, HALF_YR_SIM - 40 - FIRST_AGE
        age = HALF_YR_SIM - 40 - coh
        do n = 1, NPOP_SIM
            if (ChoiceSim(n,coh,age) == AFFSECTOR) then
                Size2 = Size2 + 1
            end if
        end do
    end do   

    allocate(coh_index(Size2))
    allocate(n_index(Size2))

    ind = 1
    do coh = HALF_YR_SIM - 40 - LAST_AGE, HALF_YR_SIM - 40 - FIRST_AGE
        age = HALF_YR_SIM - 40 - coh
        do n = 1, NPOP_SIM
           if (ChoiceSim(n,coh,age) == AFFSECTOR) then
                coh_index(ind) = coh
                n_index(ind)   = n
                ind            = ind + 1
            end if
        end do
     end do

        open(unit = 1, file = 'AffectedWorkers_Pre1.csv')
        write(1,4738) 'ID','Cohort','Year','Educ','Choice','Exper','Welfare', 'Welfare2'
        do ind = 1, Size2 
            coh = coh_index(ind)
            n   = n_index(ind)
            do t = HALF_YR_SIM - 40, HALF_YR_SIM - 40 + 35
                age = t - coh
                if (age >= FIRST_AGE .and. age <= LAST_AGE) then
                    write(1,8374) n, coh, coh+age, EducSim(n,coh), ChoiceSim(n,coh,age), & 
                                  ExperSim(n,coh,age), WelfareSim(n,coh,age)
                end if
            end do
        end do
        close(1)

    deallocate(coh_index,n_index)




    ! **********
    ! Post Shock
    ! **********

    Size2 = 0
    do coh = HALF_YR_SIM - LAST_AGE, HALF_YR_SIM - FIRST_AGE
        age = HALF_YR_SIM - coh
        do n = 1, NPOP_SIM
            if (ChoiceSim(n,coh,age) == AFFSECTOR) then
                Size2 = Size2 + 1
            end if
        end do
    end do   

    allocate(coh_index(Size2))
    allocate(n_index(Size2))

    ind = 1
    do coh = HALF_YR_SIM - LAST_AGE, HALF_YR_SIM - FIRST_AGE
        age = HALF_YR_SIM - coh
        do n = 1, NPOP_SIM
            if (ChoiceSim(n,coh,age) == AFFSECTOR) then
                coh_index(ind) = coh
                n_index(ind)   = n
                ind            = ind + 1
            end if
        end do
    end do

    open(unit = 1, file = 'AffectedWorkers_Post1.csv')
    write(1,4738) 'ID','Cohort','Year','Educ','Choice','Exper','Welfare','Welfare2'
    do ind = 1, Size2 
        coh = coh_index(ind)
        n   = n_index(ind)
        do t = HALF_YR_SIM, HALF_YR_SIM + 35
            age = t - coh
            if (age >= FIRST_AGE .and. age <= LAST_AGE) then
                write(1,8374) n, coh, coh+age, EducSim(n,coh), ChoiceSim(n,coh,age), & 
                              ExperSim(n,coh,age), WelfareSim(n,coh,age)
            end if
        end do
    end do
    close(1)
    
    4738 format(8(a20,','))
    8374 format(6(i10,','),1(f20.6,','))
    
    
    Call ComputeWelfareChanges(EducSim,ChoiceSim,ExperSim,CohortWgt,WelfareSim) 



    open(unit = 1, file = 'Simulation1.csv')
        write(1,9876) 'Year', 'PRICE1', 'PRICE2', 'PRICE3', 'PRICE4', 'PRICE5', 'PRICE6', 'PRICE7', 'PRICE INDEX', &
                      'Z_REAL1', 'Z_REAL2', 'Z_REAL3', 'Z_REAL4', 'Z_REAL5', 'Z_REAL6', 'Z_REAL7', &
                      'RSK1_0', 'RSK2_0', 'RSK3_0', 'RSK4_0', 'RSK5_0', 'RSK6_0', 'RSK7_0', &
                      'RSK1_1', 'RSK2_1', 'RSK3_1', 'RSK4_1', 'RSK5_1', 'RSK6_1', 'RSK7_1', &
                      'SK1_0', 'SK2_0', 'SK3_0', 'SK4_0', 'SK5_0', 'SK6_0', 'SK7_0', &
                      'SK1_1', 'SK2_1', 'SK3_1', 'SK4_1', 'SK5_1', 'SK6_1', 'SK7_1', &
                      'RK_1', 'RK_2', 'RK_3', 'RK_4', 'RK_5', 'RK_6', 'RK_7', & 
                      'CAPITAL1', 'CAPITAL2', 'CAPITAL3', 'CAPITAL4', 'CAPITAL5', 'CAPITAL6', 'CAPITAL7', &
                      'OUTPUT1', 'OUTPUT2', 'OUTPUT3', 'OUTPUT4', 'OUTPUT5', 'OUTPUT6', 'OUTPUT7', 'TOT_OUTPUT',&
                      'EMP0', 'EMP1', 'EMP2', 'EMP3', 'EMP4', 'EMP5', 'EMP6', 'EMP7', &
                      'EMP0_0', 'EMP1_0', 'EMP2_0', 'EMP3_0', 'EMP4_0', 'EMP5_0', 'EMP6_0', 'EMP7_0', &
                      'EMP0_1', 'EMP1_1', 'EMP2_1', 'EMP3_1', 'EMP4_1', 'EMP5_1', 'EMP6_1', 'EMP7_1', &
                      'FRACTION_0', &
                      'PERC_EMP1', 'PERC_EMP2', 'PERC_EMP3', 'PERC_EMP4', 'PERC_EMP5', 'PERC_EMP6', 'PERC_EMP7', &
                      'WELFARE', 'WELFARE0', 'WELFARE1', &
                      'RWAGE', 'RWAGE0', 'RWAGE1', &
                      'ITER', 'ITER2'
        do t = FIRST_YR_SIM, LAST_YR_SIM
            write(1,3333) t, price(1,t), price(2,t), price(3,t), price(4,t), price(5,t), price(6,t), price(7,t), Pindex(t), &
                          z(1,t), z(2,t), z(3,t), z(4,t), z(5,t), z(6,t), z(7,t), &
                          rsk_eq(1,0,t), rsk_eq(2,0,t), rsk_eq(3,0,t), rsk_eq(4,0,t), rsk_eq(5,0,t), rsk_eq(6,0,t), rsk_eq(7,0,t), & 
                          rsk_eq(1,1,t), rsk_eq(2,1,t), rsk_eq(3,1,t), rsk_eq(4,1,t), rsk_eq(5,1,t), rsk_eq(6,1,t), rsk_eq(7,1,t), &
                          sk_sup_eq(1,0,t), sk_sup_eq(2,0,t), sk_sup_eq(3,0,t), sk_sup_eq(4,0,t), sk_sup_eq(5,0,t), sk_sup_eq(6,0,t), sk_sup_eq(7,0,t), &
                          sk_sup_eq(1,1,t), sk_sup_eq(2,1,t), sk_sup_eq(3,1,t), sk_sup_eq(4,1,t), sk_sup_eq(5,1,t), sk_sup_eq(6,1,t), sk_sup_eq(7,1,t), &
                          rK_eq(t), rK_eq(t), rK_eq(t), rK_eq(t), rK_eq(t), rK_eq(t), rK_eq(t), &
                          CapitalSim(1,t), CapitalSim(2,t), CapitalSim(3,t), CapitalSim(4,t), CapitalSim(5,t), CapitalSim(6,t), CapitalSim(7,t), & 
                          output(1,t), output(2,t), output(3,t), output(4,t), output(5,t), output(6,t), output(7,t), sum(output(:,t)), &
                          emp_eq(0,t), emp_eq(1,t), emp_eq(2,t), emp_eq(3,t), emp_eq(4,t), emp_eq(5,t), emp_eq(6,t), emp_eq(7,t), &
                          emp_by_skill(0,0,t), emp_by_skill(1,0,t), emp_by_skill(2,0,t), emp_by_skill(3,0,t), emp_by_skill(4,0,t), emp_by_skill(5,0,t), emp_by_skill(6,0,t), emp_by_skill(7,0,t), &
                          emp_by_skill(0,1,t), emp_by_skill(1,1,t), emp_by_skill(2,1,t), emp_by_skill(3,1,t), emp_by_skill(4,1,t), emp_by_skill(5,1,t), emp_by_skill(6,1,t), emp_by_skill(7,1,t), &
                          emp_eq(0,t)/real(sum(emp_eq(0:7,t))), &
                          emp_eq(1,t)/real(sum(emp_eq(1:7,t))), emp_eq(2,t)/real(sum(emp_eq(1:7,t))), &
                          emp_eq(3,t)/real(sum(emp_eq(1:7,t))), emp_eq(4,t)/real(sum(emp_eq(1:7,t))), &
                          emp_eq(5,t)/real(sum(emp_eq(1:7,t))), emp_eq(6,t)/real(sum(emp_eq(1:7,t))), &
                          emp_eq(7,t)/real(sum(emp_eq(1:7,t))), &
                          Welfare(t), Welfare_Educ(0,t), Welfare_Educ(1,t), &
                          Real_Wages(t), Real_Wages_Educ(0,t), Real_Wages_Educ(1,t), &
                          iter(t), iter2(t)
        end do
    close(1)

    9876 format(106(a15,','))
    3333 format(i5,',',103(f18.5,','),2(i5,','))


    ! ****************
    ! Adjustment Costs
    ! ****************

    ! ******
    ! Output
    ! ******

    ! Steady State Real Ouput Before the Shock
    output_inf1 = sum(output(:,HALF_YR_SIM-29:HALF_YR_SIM)) / 30.0
    ! Steady State Real Ouput After the Shock
    output_inf2 = sum(output(:,LAST_YR_SIM-29:LAST_YR_SIM)) / 30.0

    ! PV Gains from Liberalization, if Reallocation is immediate
    Gains1 = (output_inf2 - output_inf1) / 0.05


    output_A = 0.0
    do t = HALF_YR_SIM + 1, LAST_YR_SIM
        output_A = output_A + (0.95**(t-HALF_YR_SIM-1))*sum(output(:,t))
    end do
    output_A = output_A + (output_inf2 / 0.05)*(0.95**(LAST_YR_SIM-HALF_YR_SIM))

    ! PV Gains from Liberalization along the adjustment path
    Gains2 = output_A - output_inf1 / 0.05

    ! Adjustment Cost
    AdjCost = 100.0*(Gains1 - Gains2)/Gains1
   

    ! *******
    ! Welfare
    ! *******

    ! ******************
    ! **** Overall *****
    ! ******************

    ! Steady State Welfare Before the Shock
    Welfare_inf1 = sum(Welfare(HALF_YR_SIM-29:HALF_YR_SIM)) / 30.0
    ! Steady State Welfare After the Shock
    Welfare_inf2 = sum(Welfare(LAST_YR_SIM-29:LAST_YR_SIM)) / 30.0

    ! PV Gains from Liberalization, if Reallocation is immediate
    WGains1 = (Welfare_inf2 - Welfare_inf1) / 0.05

    Welfare_A = 0.0
    do t = HALF_YR_SIM + 1, LAST_YR_SIM
        Welfare_A = Welfare_A + (0.95**(t-HALF_YR_SIM-1))*Welfare(t)
    end do
    Welfare_A = Welfare_A + (Welfare_inf2 / 0.05)*(0.95**(LAST_YR_SIM-HALF_YR_SIM))

    ! PV Gains from Liberalization along the adjustment path
    WGains2 = Welfare_A - Welfare_inf1 / 0.05

    ! Adjustment Cost
    WAdjCost = 100.0*(WGains1 - WGains2)/WGains1



    ! *********************************
    ! **** Educated / NonEducated *****
    ! *********************************

    ! Steady State Welfare Before the Shock
    Welfare_inf_Educ1 = sum(Welfare_Educ(1,HALF_YR_SIM-29:HALF_YR_SIM)) / 30.0
    ! Steady State Welfare After the Shock
    Welfare_inf_Educ2 = sum(Welfare_Educ(1,LAST_YR_SIM-29:LAST_YR_SIM)) / 30.0

    ! PV Gains from Liberalization, if Reallocation is immediate
    WGains_Educ1 = (Welfare_inf_Educ2 - Welfare_inf_Educ1) / 0.05

    Welfare_EducA = 0.0
    do t = HALF_YR_SIM + 1, LAST_YR_SIM
        Welfare_EducA = Welfare_EducA + (0.95**(t-HALF_YR_SIM-1))*Welfare_Educ(1,t)
    end do
    Welfare_EducA = Welfare_EducA + (Welfare_inf_Educ2 / 0.05)*(0.95**(LAST_YR_SIM-HALF_YR_SIM))

    ! PV Gains from Liberalization along the adjustment path
    WGains_Educ2 = Welfare_EducA - Welfare_inf_Educ1 / 0.05

    ! Adjustment Cost
    WAdjCost_Educ = 100.0*(WGains_Educ1 - WGains_Educ2)/WGains_Educ1



    ! Steady State Welfare Before the Shock
    Welfare_inf_NonEduc1 = sum(Welfare_Educ(0,HALF_YR_SIM-29:HALF_YR_SIM)) / 30.0
    ! Steady State Welfare After the Shock
    Welfare_inf_NonEduc2 = sum(Welfare_Educ(0,LAST_YR_SIM-29:LAST_YR_SIM)) / 30.0

    ! PV Gains from Liberalization, if Reallocation is immediate
    WGains_NonEduc1 = (Welfare_inf_NonEduc2 - Welfare_inf_NonEduc1) / 0.05

    Welfare_NonEducA = 0.0
    do t = HALF_YR_SIM + 1, LAST_YR_SIM
        Welfare_NonEducA = Welfare_NonEducA + (0.95**(t-HALF_YR_SIM-1))*Welfare_Educ(0,t)
    end do
    Welfare_NonEducA = Welfare_NonEducA + (Welfare_inf_NonEduc2 / 0.05)*(0.95**(LAST_YR_SIM-HALF_YR_SIM))

    ! PV Gains from Liberalization along the adjustment path
    WGains_NonEduc2 = Welfare_NonEducA - Welfare_inf_NonEduc1 / 0.05

    ! Adjustment Cost
    WAdjCost_NonEduc = 100.0*(WGains_NonEduc1 - WGains_NonEduc2)/WGains_NonEduc1



    ! **********
    ! Real Wages
    ! **********

    ! ******************
    ! **** Overall *****
    ! ******************

    ! Steady State Welfare Before the Shock
    Wage_inf1 = sum(Real_Wages(HALF_YR_SIM-29:HALF_YR_SIM)) / 30.0
    ! Steady State Welfare After the Shock
    Wage_inf2 = sum(Real_Wages(LAST_YR_SIM-29:LAST_YR_SIM)) / 30.0

    ! PV Gains from Liberalization, if Reallocation is immediate
    WageGains1 = (Wage_inf2 - Wage_inf1) / 0.05

    Wage_A = 0.0
    do t = HALF_YR_SIM + 1, LAST_YR_SIM
        Wage_A = Wage_A + (0.95**(t-HALF_YR_SIM-1))*Real_Wages(t)
    end do
    Wage_A = Wage_A + (Wage_inf2 / 0.05)*(0.95**(LAST_YR_SIM-HALF_YR_SIM))

    ! PV Gains from Liberalization along the adjustment path
    WageGains2 = Wage_A - Wage_inf1 / 0.05

    ! Adjustment Cost
    WageAdjCost = 100.0*(WageGains1 - WageGains2)/WageGains1



    ! *********************************
    ! **** Educated / NonEducated *****
    ! *********************************

    ! Steady State Welfare Before the Shock
    Wage_inf_Educ1 = sum(Real_Wages_Educ(1,HALF_YR_SIM-29:HALF_YR_SIM)) / 30.0
    ! Steady State Welfare After the Shock
    Wage_inf_Educ2 = sum(Real_Wages_Educ(1,LAST_YR_SIM-29:LAST_YR_SIM)) / 30.0

    ! PV Gains from Liberalization, if Reallocation is immediate
    WageGains_Educ1 = (Wage_inf_Educ2 - Wage_inf_Educ1) / 0.05

    Wage_EducA = 0.0
    do t = HALF_YR_SIM + 1, LAST_YR_SIM
        Wage_EducA = Wage_EducA + (0.95**(t-HALF_YR_SIM-1))*Real_Wages_Educ(1,t)
    end do
    Wage_EducA = Wage_EducA + (Wage_inf_Educ2 / 0.05)*(0.95**(LAST_YR_SIM-HALF_YR_SIM))

    ! PV Gains from Liberalization along the adjustment path
    WageGains_Educ2 = Wage_EducA - Wage_inf_Educ1 / 0.05

    ! Adjustment Cost
    WageAdjCost_Educ = 100.0*(WageGains_Educ1 - WageGains_Educ2)/WageGains_Educ1



    ! Steady State Welfare Before the Shock
    Wage_inf_NonEduc1 = sum(Real_Wages_Educ(0,HALF_YR_SIM-29:HALF_YR_SIM)) / 30.0
    ! Steady State Welfare After the Shock
    Wage_inf_NonEduc2 = sum(Real_Wages_Educ(0,LAST_YR_SIM-29:LAST_YR_SIM)) / 30.0

    ! PV Gains from Liberalization, if Reallocation is immediate
    WageGains_NonEduc1 = (Wage_inf_NonEduc2 - Wage_inf_NonEduc1) / 0.05

    Wage_NonEducA = 0.0
    do t = HALF_YR_SIM + 1, LAST_YR_SIM
        Wage_NonEducA = Wage_NonEducA + (0.95**(t-HALF_YR_SIM-1))*Real_Wages_Educ(0,t)
    end do
    Wage_NonEducA = Wage_NonEducA + (Wage_inf_NonEduc2 / 0.05)*(0.95**(LAST_YR_SIM-HALF_YR_SIM))

    ! PV Gains from Liberalization along the adjustment path
    WageGains_NonEduc2 = Wage_NonEducA - Wage_inf_NonEduc1 / 0.05

    ! Adjustment Cost
    WageAdjCost_NonEduc = 100.0*(WageGains_NonEduc1 - WageGains_NonEduc2)/WageGains_NonEduc1



    open(unit = 1, file = 'AdjCosts1.csv')
        write(1,9090) ' ', 'Real Output', 'Welfare', 'RWages'
        write(1,8967) 'PV Welfare_inf1:', output_inf1 / 0.05, Welfare_inf1 / 0.05, Wage_inf1 / 0.05
        write(1,8967) 'PV Welfare_inf2:', output_inf2 / 0.05, Welfare_inf2 / 0.05, Wage_inf2 / 0.05
        write(1,8967) 'PV Welfare_A:', output_A, Welfare_A, Wage_A
        write(1,8967) 'Gains1:', Gains1, WGains1, WageGains1
        write(1,8967) 'Gains2:', Gains2, WGains2, WageGains2
        write(1,8967) 'Long Term Welfare Gain:', 100.0*(output_inf2 - output_inf1)/output_inf1, &
                                    100.0*(Welfare_inf2 - Welfare_inf1)/Welfare_inf1, &
                                    100.0*(Wage_inf2 - Wage_inf1)/Wage_inf1
        write(1,8967) 'AdjCost:', AdjCost, WAdjCost, WageAdjCost
        write(1,8967) 
        write(1,8967) 'Educated'
        write(1,8967) 'PV Welfare_inf1:', 0.0, Welfare_inf_Educ1 / 0.05, Wage_inf_Educ1 / 0.05
        write(1,8967) 'PV Welfare_inf2:', 0.0, Welfare_inf_Educ2 / 0.05, Wage_inf_Educ2 / 0.05
        write(1,8967) 'PV Welfare_A:', 0.0, Welfare_EducA, Wage_EducA
        write(1,8967) 'Gains1:', 0.0, WGains_Educ1, WageGains_Educ1
        write(1,8967) 'Gains2:', 0.0, WGains_Educ2, WageGains_Educ2
        write(1,8967) 'Long Term Welfare Gain:', 0.0, &
                      100.0*(Welfare_inf_Educ2 - Welfare_inf_Educ1)/Welfare_inf_Educ1, &
                      100.0*(Wage_inf_Educ2 - Wage_inf_Educ1)/Wage_inf_Educ1
        write(1,8967) 'AdjCost:', 0.0, WAdjCost_Educ, WageAdjCost_Educ
        write(1,8967)
        write(1,8967) 'NonEducated'
        write(1,8967) 'PV Welfare_inf1:', 0.0, Welfare_inf_NonEduc1 / 0.05, Wage_inf_NonEduc1 / 0.05
        write(1,8967) 'PV Welfare_inf2:', 0.0, Welfare_inf_NonEduc2 / 0.05, Wage_inf_NonEduc2 / 0.05
        write(1,8967) 'PV Welfare_A:', 0.0, Welfare_NonEducA, Wage_NonEducA
        write(1,8967) 'Gains1:', 0.0, WGains_NonEduc1, WageGains_NonEduc1
        write(1,8967) 'Gains2:', 0.0, WGains_NonEduc2, WageGains_NonEduc2
        write(1,8967) 'Long Term Welfare Gain:', 0.0, &
                      100.0*(Welfare_inf_NonEduc2 - Welfare_inf_NonEduc1)/Welfare_inf_NonEduc1, &
                      100.0*(Wage_inf_NonEduc2 - Wage_inf_NonEduc1)/Wage_inf_NonEduc1
        write(1,8967) 'AdjCost:', 0.0, WAdjCost_NonEduc, WageAdjCost_NonEduc
    close(1)    

    9090 format(4(a20,','))
    8967 format(a20,',',3(f20.6,','))



    ! ************
    ! Reallocation
    ! ************

    ! Looking at Real Output
    do t = HALF_YR_SIM, LAST_YR_SIM
        perc_complete(t) = 100.0*(sum(output(:,t)) - output_inf1) / (output_inf2 - output_inf1) 
    end do 



    ! Looking at employment at the affected sector
    emp1 = 0.0
    emp2 = 0.0

    do t = HALF_YR_SIM-29, HALF_YR_SIM
         emp1 = emp1 + emp_eq(AFFSECTOR,t)    
    end do
    emp1 = emp1 / 30.0


    do t = LAST_YR_SIM-29, LAST_YR_SIM
       emp2 = emp2 + emp_eq(AFFSECTOR,t)    
    end do
    emp2 = emp2 / 30.0

    do t = HALF_YR_SIM, LAST_YR_SIM
        perc_complete2(t) = 100.0*(emp_eq(AFFSECTOR,t) - emp1) / (emp2 - emp1)
    end do    


    open(unit = 1, file = 'Reallocation1.csv')
        write(1,8888) 'Year', 'Output', 'Affected Sector'
        do t = HALF_YR_SIM, LAST_YR_SIM 
            write(1,7777) t, perc_complete(t), perc_complete2(t)
        end do 
    close(1)

    8888 format(3(a20,','))
    7777 format(i5,',',2(f10.6,','))



    deallocate(coh_index,n_index)
    deallocate(PI_Myopic, R_sq_Myopic, PI_RE, R_sq_RE)
    deallocate(ChoiceInit,ChoiceSim)
    deallocate(EducInit,EducSim)
    deallocate(GenderInit,GenderSim)
    deallocate(eps)
    deallocate(WelfareSim,ExperSim)
    
end if    
                                            
end subroutine Simulation1


end module Simulation_MOD1

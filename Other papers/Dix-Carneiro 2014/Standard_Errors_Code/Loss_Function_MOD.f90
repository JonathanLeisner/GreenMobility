MODULE Loss_Function_MOD


Contains


! **************************************************************************
! SMM_OBJ_FCN simulates the economy, calculates simulated moments and
! in the end calculates the Loss Function, which is the distance between the 
! Data Moments and the Simulated Moments
! WriteMomentsToFile is a routine that write the Simulated Moments vs 
! The actual data moments in a file
! Summarizing, this module contains 3 subroutines/functions:
! subroutine SMM_OBJ_FCN(ParamScaled,F)
! function   Loss_Function(Dataset,choice_index,lag_choice_index,&
!                          sector_size,lag_sector_size)
! subroutine WriteCoefToFile(betaSim,SigmaWageSim,gammaSim,phiSim,&
!                            Loss_Function)
! subroutine Write_Sim_Data(ChoiceSim,WageSim,ExperSim,CohortWgt, &
!                           CapitalSim,z,sk_sup_eq,rsk_eq)
! **************************************************************************

Subroutine SMM_OBJ_FCN(ParamScaled,F)

USE Global_data
USE Emax_MOD
USE MPI
USE Newton
USE ParallelCohorts_MOD

implicit none

real(KIND=DOUBLE) test1, test2, test3, test4, test5, test6
 
real(KIND=DOUBLE) , intent(in)  :: ParamScaled(:)
real(KIND=DOUBLE) , intent(out) :: F


! Model parameters
real(KIND=DOUBLE) param(NPARAM), theta(7), beta(NSECTORS,12), & 
                  sigma(0:NSECTORS), kappa(26), tau(0:NSECTORS), SigmaPref, &
                  omega(2:NTYPES,0:NSECTORS), lambda(NTYPES), gamma(2:NTYPES,0:NSECTORS+4), &
                  aa(NSECTORS,0:1,0:1), alpha_prod(NSECTORS,0:1,FIRST_YR:LAST_YR), sigma_prod(NSECTORS),rsk_init(NSECTORS,0:1),&
                  alpha_EXP(NSECTORS,0:1), beta_EXP(NSECTORS,0:1), sigma2_EXP(NSECTORS,0:1), alpha_CB(NSECTORS,0:1,FIRST_YR:LAST_YR)


integer           i, ii, j, n, t, k, l, s, s1, s2, it, year, age, gender, ed, &
                  coh, obs, educ, educ2, educ3, educ4, RE_loop, tp, &
                  size1998, size2000, size2005, info, dim1, dim2, FLAG_STOP, FLAG_EMAX, SUM_FLAG_EMAX, coh_batch, older_coh, &
                  newer_coh
                  
        
! Iteration Loop for the Equilibrium 
real(KIND=DOUBLE) rsk_sup(NSECTORS,0:1,MAXIT), &
                  rsk_dem(NSECTORS,0:1,MAXIT), &
                  rsk_eq(NSECTORS,0:1,FIRST_YR:LAST_YR), &
                  rsk_ratio(NSECTORS,0:1,FIRST_YR:LAST_YR-1), &
                  rsk_tomorrow(NSECTORS), &
                  rsk_lag(NSECTORS,0:1) , &
                  new_rsk_tom(NSECTORS,0:1), &
                  sk_sup(NSECTORS,0:1), &
                  sk_supOut(NSECTORS,0:1), &
                  sk_sup_coh(NSECTORS,0:1), &
                  sk_sup_eq(NSECTORS,0:1,FIRST_YR:LAST_YR), &
                  cost(0:NSECTORS,0:NSECTORS), & 
                  Emax(0:NSECTORS), & 
                  Emax0(0:NSECTORS), & 
                  w(0:NSECTORS), & 
                  V(0:NSECTORS), & 
                  VMAX, &
                  check(NSECTORS,0:1,FIRST_YR:LAST_YR), & 
                  step(NSECTORS,0:1), & 
                  CohortWgt(FIRST_COH:LAST_COH,0:1), &
                  CohortWgtOut(COH_PROC,0:1), &
                  weight, & 
                  eps_aux(0:NSECTORS), &
                  eta_aux(0:NSECTORS), &
                  crit, &
                  cut(NSECTORS,0:1), &
                  term2, &
                  term3, &
                  prob(3), &
                  COEF(2,1), &
                  COEF_EXP(NSECTORS,0:1,2,1), &
                  SST, SSE, &
                  VY(9,1), &
                  VX(9,2), &
                  sigma_EXP(NSECTORS,0:1), &
                  covar_EXP(NSECTORS,NSECTORS,0:1), &
                  resid(NSECTORS,0:1,9), &
                  const1, const2, const3, const4, &
                  Capital(NSECTORS), &
                  solution, &
                  init, &
                  LaborSharesSim(NSECTORS,0:1,FIRST_YR:LAST_YR), &
                  CapitalSharesSim(NSECTORS,FIRST_YR:LAST_YR)
                  
real(KIND=DOUBLE) A(FIRST_AGE+1:LAST_AGE,0:NSECTORS,NREG), &
                  B(FIRST_AGE+1:LAST_AGE,0:NSECTORS), &
                  A_RE(FIRST_AGE+1:LAST_AGE,0:NSECTORS,NREG), &
                  B_RE(FIRST_AGE+1:LAST_AGE,0:NSECTORS)                  
                  
integer           emp(0:NSECTORS,0:1,MAXIT), &
                  empOut(0:NSECTORS,0:1), &
                  emp_coh(0:NSECTORS,0:1), &
                  emp_eq(0:NSECTORS,0:1,FIRST_YR:LAST_YR), &
                  iter(FIRST_YR:LAST_YR), & 
                  flag(NSECTORS,0:1), & 
                  lag, & 
                  exper(NSECTORS), &
                  ExperTomorrow(NSECTORS), &
                  dummy(NSECTORS,0:1), &
                  typ(NPOP,FIRST_COH:LAST_COH) 
    
real(KIND=DOUBLE) AggChoices(0:NSECTORS), AggChoices2(0:NSECTORS), &
                  AggLogWages(NSECTORS), AggTr(0:NSECTORS,0:NSECTORS), &
                  SumWeightChoices
                  

! Vectors of sizes and indices, used in the estimation
! Of the auxiliary models
integer           Size2, &
                  ind(0:NSECTORS), &
                  lag_ind(0:NSECTORS), &
                  switch_ind(0:NSECTORS), &
                  switch_size(0:NSECTORS), &
                  ind2(0:NSECTORS), &
                  sector_size(0:NSECTORS), &
                  lag_sector_size(0:NSECTORS), &
                  wagedif_size(0:NSECTORS)      
                  

logical           converged                                            


real(KIND=DOUBLE) CapitalSim(NSECTORS,FIRST_YR:LAST_YR), &
                  z(NSECTORS,FIRST_YR:LAST_YR)                  
          
integer          , allocatable, dimension(:,:,:)   :: ChoiceSim, ChoiceSimAux, ChoiceSimOut
integer          , allocatable, dimension(:,:,:,:) :: ExperSim, ExperSimAux, ExperSimOut
real(KIND=DOUBLE), allocatable, dimension(:,:,:)   :: WageSim, WageSimAux, WageSimOut
real(KIND=DOUBLE), allocatable, dimension(:,:,:)   :: Cost1Sim, Cost2Sim, VSim, & 
                                                      Cost1SimOut, Cost2SimOut, VSimOut
real(KIND=DOUBLE), allocatable, dimension(:,:,:,:) :: eps, epsAux       
real(KIND=DOUBLE), allocatable, dimension(:,:,:,:) :: eta, etaAux 
                          
                  
real(KIND=DOUBLE), allocatable, dimension(:,:)   :: Dataset
integer          , allocatable, dimension(:,:)   :: choice_index, & 
                                                    lag_choice_index, &
                                                    switch_index, &
                                                    index_wagedif
integer          , allocatable, dimension(:,:,:) :: index_init
real(KIND=DOUBLE), allocatable, dimension(:,:,:) :: frac
     
! Coefficients for the calculation of Emax_hat
real(KIND=DOUBLE), allocatable, dimension(:,:,:,:,:,:)   :: PI_Myopic
real(KIND=DOUBLE), allocatable, dimension(:,:,:,:,:)     :: R_sq_Myopic
real(KIND=DOUBLE), allocatable, dimension(:,:,:,:,:,:)   :: PI_RE 
real(KIND=DOUBLE), allocatable, dimension(:,:,:,:,:)     :: R_sq_RE
                  
real(KIND=DOUBLE) PI_COEF(NREG)                    
          
! Elapsed time variables
integer           time_array_0(8), time_array_1(8)
real(KIND=DOUBLE) start_time, finish_time, &
                  start_time_emax, finish_time_emax, &
                  start_time_eq, finish_time_eq, &
                  start_time_loss, finish_time_loss

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
allocate(PI_RE(FIRST_AGE+1:LAST_AGE,0:NSECTORS,1:CAT_GENDER,1:CAT_EDUC,1:NTYPES,NREG))
allocate(R_sq_RE(FIRST_AGE+1:LAST_AGE,0:NSECTORS,1:CAT_GENDER,1:CAT_EDUC,1:NTYPES))    

allocate(ChoiceSimAux(NPOP,COH_PROC,-9:-1)) 
allocate(ExperSimAux(NSECTORS,NPOP,COH_PROC,1))
allocate(WageSimAux(NPOP,COH_PROC,1))
allocate(epsAux(NPOP,COH_PROC,1,0:NSECTORS))
allocate(etaAux(NPOP,COH_PROC,1,0:NSECTORS))

allocate(ChoiceSimOut(NPOP,COH_PROC,1))
allocate(ExperSimOut(NSECTORS,NPOP,COH_PROC,1))
allocate(WageSimOut(NPOP,COH_PROC,1))
allocate(Cost1SimOut(NPOP,COH_PROC,1))
allocate(Cost2SimOut(NPOP,COH_PROC,1))
allocate(VSimOut(NPOP,COH_PROC,1))


FLAG_STOP = 0
FLAG_EMAX = 0
SUM_FLAG_EMAX = 0
    

! *****************************
! Use only the Master processor
! *****************************

rankcond1: if (rank == 0) then        
                  
    !***********************
    ! Initializing the clock
    !***********************

    call date_and_time(values=time_array_0)
    start_time = time_array_0(5) * 3600 + time_array_0(6) * 60 &
               + time_array_0(7) + 0.001 * time_array_0(8)


    allocate(ChoiceSim(NPOP,FIRST_COH:LAST_COH,(FIRST_YR-9):LAST_YR))  
    allocate(Cost1Sim(NPOP,FIRST_COH:LAST_COH,FIRST_YR:LAST_YR))
    allocate(Cost2Sim(NPOP,FIRST_COH:LAST_COH,FIRST_YR:LAST_YR))
    allocate(VSim(NPOP,FIRST_COH:LAST_COH,FIRST_YR:LAST_YR))
    allocate(ExperSim(NSECTORS,NPOP,FIRST_COH:LAST_COH,FIRST_YR:LAST_YR))
    allocate(WageSim(NPOP,FIRST_COH:LAST_COH,FIRST_YR:LAST_YR))
    allocate(eps(NPOP,FIRST_COH:LAST_COH,FIRST_YR:LAST_YR,0:NSECTORS))
    allocate(eta(NPOP,FIRST_COH:LAST_COH,FIRST_YR:LAST_YR,0:NSECTORS))

    ! ****************************
    ! Function evaluations counter
    ! ****************************

    FUNCTION_ITER  = FUNCTION_ITER + 1

    ! *************************
    ! Re-scaling the parameters
    ! *************************

    param = ParamScaled

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
    rsk_init(1:7,0) = 1.0
    rsk_init(1:7,1) = 1.0

    do ed = 0, 1
        do s = 1, NSECTORS
            do t = FIRST_YR, LAST_YR
                alpha_prod(s,ed,t) = aa(s,ed,0) + aa(s,ed,1)*(t-FIRST_YR)
                if (alpha_prod(s,ed,t) <= 0 .or. alpha_prod(s,ed,t) >= 1) then
                    FLAG_STOP = 1
                end if
            end do
        end do
    end do


    if (FLAG_STOP == 1) then
    
        F = 999999.9999
    
        print*, 'FUNCTION ITERATION: ', FUNCTION_ITER
        print*
        print*, '*****************'
        print*, 'Loss Function'
        print*, F
        print*, '*****************'
        print*
        print*, 'Negative Labor Share'
        print*
        print* , '====================================================================='
        print* , '====================================================================='
        print*
        print*
    
    else

        ChoiceSim   = -999 
        WageSim     = -999 
        ExperSim    = -999

        ! Assign initial conditions (FIRST_YR = 1995) to ChoiceSim
        do coh = FIRST_COH, (LAST_YR - FIRST_AGE)
            do n = 1, NPOP
                age = FirstYearData(n,coh) - coh
                do i = 1, 9
                    ChoiceSim(n,coh,coh+age-i) = ChoiceData(n,coh,coh+age-i)
                end do 
            end do
        end do
    
    
        ! *****************************************
        ! Idiosynchractic shocks for the simulation
        ! *****************************************
    
        do s = 0, NSECTORS
            eps(:,:,:,s) = sigma(s)*eps_global(:,:,:,s)
            eta(:,:,:,s) = -sigmaPref*eta_global(:,:,:,s) - sigmaPref*0.577215665
        end do  
        
    end if
    
end if rankcond1

! ******************************************************
! Broadcast to other processors the vector of parameters
! ******************************************************

Call MPI_BCAST(param,NPARAM,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)    

Call MPI_BCAST(FLAG_STOP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)



flagstop: if (FLAG_STOP == 0) then
    ! Otherwise jumps to the very end of the routine    

    ! ****************************************************
    ! Upper and Lower Bounds for draws of returns to skill
    ! Used in the calculation of Emax
    ! These are global variables
    ! ****************************************************


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

    ! **********************************************************************************************************************************
    ! EMAX Myopic : First Interation
    ! **********************************************************************************************************************************

    ! ******************************
    ! Compute Emax in each processor
    ! ******************************

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

    if (rank == 0) then
        call date_and_time(values=time_array_0)
        start_time_emax = time_array_0(5) * 3600 + time_array_0(6) * 60 &
                        + time_array_0(7) + 0.001 * time_array_0(8)
        print*
        print*, '************************'
        print*, 'Computing Emax Myopic...'
        print*, '************************'
        print*
    end if

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
        
        call date_and_time(values=time_array_0)
        finish_time_emax = time_array_0(5) * 3600 + time_array_0(6) * 60 &
                         + time_array_0(7) + 0.001 * time_array_0(8)
        print*
        print*, 'Emax Time     :', real(finish_time_emax - start_time_emax,4)
        print*
        print*, '********************'
        print*, 'Emax Myopic Computed'
        print*, '********************'
        print*   
        
        open(unit = 1, file = 'Emax_coef.csv')
            write(1,4646) 'Age, Lag, Gender, Educ, Type, Reg, PI_Myopic'
            do age = FIRST_AGE+1, LAST_AGE
                do lag = 0, NSECTORS
                    do gender = 1, CAT_GENDER
                        do educ = 1, CAT_EDUC
                            do tp = 1, NTYPES
                                do i = 1, NREG
                                    write(1,4545) age, lag, gender, educ, tp, i, PI_Myopic(age,lag,gender,educ,tp,i)
                                end do
                            end do
                        end do
                    end do
                end do
            end do        
        close(1)        
        4646 format(a50)
        4545 format(6(i6, ','),f18.8)
             
    end if   

    count1 = size(PI_Myopic)
    Call MPI_BCAST(PI_Myopic,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    ! **********************************************************************************************************************************
    ! End EMAX Myopic 
    ! **********************************************************************************************************************************


    if (rank == 0) then
        
        call date_and_time(values=time_array_0)
        start_time_eq = time_array_0(5) * 3600 + time_array_0(6) * 60 &
                      + time_array_0(7) + 0.001 * time_array_0(8)
    
        do coh = FIRST_COH, LAST_COH
            do n = 1, NPOP
                educ    = EducData(n,coh)
                educ2 = 0
                educ3 = 0
                educ4 = 0
                if (educ == 2) then
                    educ2 = 1
                elseif (educ == 3) then
                    educ3 = 1
                elseif (educ == 4) then
                    educ4 = 1
                end if
                gender  = GenderData(n,coh)
                age   = FirstYearData(n,coh) - coh
                exper = 0
                do i = 1, 9
                    do s = 1, NSECTORS                    
                        if (ChoiceData(n,coh,coh+age-i) == s) then
                            exper(s) = exper(s) + 1
                        end if                    
                    end do
                end do        
                term2 = exp(gamma(2,0) + gamma(2,1)*exper(1) + gamma(2,2)*exper(2) + &
                            gamma(2,3)*exper(3) + gamma(2,4)*exper(4) + gamma(2,5)*exper(5) + &
                            gamma(2,6)*exper(6) + gamma(2,7)*exper(7) + gamma(2,8)*(gender-1) + &
                            gamma(2,9)*educ2 + gamma(2,10)*educ3 + gamma(2,11)*educ4)
                term3 = exp(gamma(3,0) + gamma(3,1)*exper(1) + gamma(3,2)*exper(2) + &
                            gamma(3,3)*exper(3) + gamma(3,4)*exper(4) + gamma(3,5)*exper(5) + &
                            gamma(3,6)*exper(6) + gamma(3,7)*exper(7) + gamma(3,8)*(gender-1) + &
                            gamma(3,9)*educ2 + gamma(3,10)*educ3 + gamma(3,11)*educ4)
                prob(1) = 1.0    / (1.0 + term2 + term3)
                prob(2) = term2  / (1.0 + term2 + term3)
                prob(3) = term3  / (1.0 + term2 + term3)
                if (unif_global(n,coh,1) < prob(1)) then
                    typ(n,coh) = 1
                else if (unif_global(n,coh,1) >= prob(1) .and. unif_global(n,coh,1) < prob(1) + prob(2)) then
                    typ(n,coh) = 2
                else if (unif_global(n,coh,1) >= prob(1) + prob(2) .and. unif_global(n,coh,1) <= 1) then
                    typ(n,coh) = 3
                end if
            end do
        end do

    end if   

    count1 = size(typ)
    Call MPI_BCAST(typ,count1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 


    RE_loop = 1

    RatExp_loop: do while (RE_loop <= 2)

        reloop_gt_2: if (RE_loop >= 2) then

            ! **********************************************************************************************************************************
            ! EMAX Rational Expectations
            ! **********************************************************************************************************************************
            
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

            if (rank == 0) then
                call date_and_time(values=time_array_0)
                start_time_emax = time_array_0(5) * 3600 + time_array_0(6) * 60 &
                                + time_array_0(7) + 0.001 * time_array_0(8)
                print*
                print*, '***************************************'
                print*, 'Computing Emax Rational Expectations...'
                print*, '***************************************'
                print*
            end if
            
            ! *********************************
            ! Compute Emax_RE in each processor
            ! *********************************
            count1 = size(A)
            count2 = size(B)

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
        
            Call EmaxCoef_RE(param,COEF_EXP,covar_EXP,cut,gender,educ,tp,A,B,FLAG_EMAX)
        
            if (rank /= 0) then
                tag = 1
                Call MPI_SEND(FLAG_EMAX,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,ierr)
            end if
        
            if (rank == 0) then
                SUM_FLAG_EMAX = FLAG_EMAX
                do proc = 1, numproc-1
                    tag = 1
                    Call MPI_RECV(FLAG_EMAX,1,MPI_INTEGER,proc,tag,MPI_COMM_WORLD,status,ierr)
                    SUM_FLAG_EMAX = SUM_FLAG_EMAX + FLAG_EMAX
                end do
                print*
                print*, 'SUM_FLAG_EMAX =', SUM_FLAG_EMAX
                print*
            end if
        
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
            Call MPI_BCAST(SUM_FLAG_EMAX,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            
            flag_emax_cond1: if(SUM_FLAG_EMAX == 0) then

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
                            PI_RE(age,lag,gender,educ,tp,:) = A(age,lag,:)
                            R_sq_RE(age,lag,gender,educ,tp) = B(age,lag)
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
                                PI_RE(age,lag,gender,educ,tp,:) = A(age,lag,:)
                                R_sq_RE(age,lag,gender,educ,tp) = B(age,lag)
                            end do
                        end do
                    end do
                    call date_and_time(values=time_array_0)
                    finish_time_emax = time_array_0(5) * 3600 + time_array_0(6) * 60 &
                                     + time_array_0(7) + 0.001 * time_array_0(8)
                    print*
                    print*, 'Emax Time     :', real(finish_time_emax - start_time_emax,4)
                    print*
                    print*, '***********************************'
                    print*, 'Emax Rational Expectations Computed'
                    print*, '***********************************'
                    print*

                end if
                ! End if rank == 0

                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
                count1 = size(PI_RE)
                Call MPI_BCAST(PI_RE,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    
            end if flag_emax_cond1
            ! End of SUM_FLAG_EMAX == 0 Conditional
    
        end if reloop_gt_2
        ! End of Conditional RE_loop >= 2
        
    
        ! If EMAX was in trouble for any type,
        ! Then exit the while do statements
        if (SUM_FLAG_EMAX > 0) then
            exit
            ! It goes to the end of RatExp_loop
        end if
    
   
        if (rank == 0) then    
            if (FUNCTION_ITER == 1) then
                rsk_global(:,0) = exp(sum(betaData(:,1:(LAST_YR-FIRST_YR+1)),2) / real((LAST_YR-FIRST_YR+1)))
                rsk_global(:,1) = exp(sum(betaData(:,1:(LAST_YR-FIRST_YR+1)),2) / real((LAST_YR-FIRST_YR+1)))
            end if            
        end if
        
        Call MPI_BCAST(rsk_global,size(rsk_global),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
        rsk_sup         = 0.0
        rsk_dem         = 0.0
        rsk_sup(:,:,1)  = rsk_global 
        
        crit = 0.001
        
        rsk_eq    = 0.0        
        sk_sup_eq = 0.0
        emp_eq    = 0    
        iter      = 0
        step      = 0.0    

        year_loop: do t = FIRST_YR, LAST_YR

            ! ********************************
            ! Initialize the iteration counter   
            ! ********************************
            it = 1
    
         
            if (t > FIRST_YR) then
                ! Use the result of the previous year
                ! As the initial point for rsk
                rsk_sup        = 0.0
                rsk_dem        = 0.0        
                rsk_sup(:,0,1) = min(max(rsk_eq(:,0,t-1),LOW_RSK(:,0)),UP_RSK(:,0))
                rsk_sup(:,1,1) = min(max(rsk_eq(:,1,t-1),LOW_RSK(:,1)),UP_RSK(:,1))
            end if


            ! These flag variables are used in order to check whether 
            ! we need to change the rsk updating step
            flag = 0
    
    
    
            converged = .false.

            while_loop: do while (.not. converged .and. it <= MAXIT)
                
                year        = t
                sk_sup      = 0.0
                emp(:,:,it) = 0
              
                do coh_batch = 2, NCOH_BATCH
                
                    if (rank == 0) then
                    
                        older_coh = (year-LAST_AGE) + (coh_batch-1)*COH_PROC
                        newer_coh = (year-LAST_AGE) + coh_batch*COH_PROC - 1
                        
                        ChoiceSimAux = ChoiceSim(:,older_coh:newer_coh,year-9:year-1)
                        epsAux = eps(:,older_coh:newer_coh,year:year,:)
                        etaAux = eta(:,older_coh:newer_coh,year:year,:)

                        tag = coh_batch*10 + 1                        
                        count1 = size(ChoiceSimAux) 
                        Call MPI_SEND(ChoiceSimAux,count1,MPI_INTEGER,coh_batch-1,tag,MPI_COMM_WORLD,status,ierr)

                        tag = coh_batch*10 + 2
                        count1 = size(epsAux)
                        Call MPI_SEND(epsAux,count1,MPI_DOUBLE_PRECISION,coh_batch-1,tag,MPI_COMM_WORLD,status,ierr)

                        tag = coh_batch*10 + 3
                        count1 = size(etaAux)
                        Call MPI_SEND(etaAux,count1,MPI_DOUBLE_PRECISION,coh_batch-1,tag,MPI_COMM_WORLD,status,ierr)

                    end if

                    if (rank == coh_batch-1) then

                        tag = coh_batch*10 + 1
                        count1 = size(ChoiceSimAux)
                        Call MPI_RECV(ChoiceSimAux,count1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)

                        tag = coh_batch*10 + 2
                        count1 = size(epsAux)
                        Call MPI_RECV(epsAux,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)

                        tag = coh_batch*10 + 3
                        count1 = size(etaAux)
                        Call MPI_RECV(etaAux,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)

                    end if

                end do


                if (rank == 0) then
                    
                    coh_batch = 1
                    older_coh = (year-LAST_AGE) + (coh_batch-1)*COH_PROC
                    newer_coh = (year-LAST_AGE) + coh_batch*COH_PROC - 1
                    
                    ChoiceSimAux = ChoiceSim(:,older_coh:newer_coh,year-9:year-1)
                    epsAux = eps(:,older_coh:newer_coh,year:year,:)
                    etaAux = eta(:,older_coh:newer_coh,year:year,:)
                
                end if
                
                older_coh = (year-LAST_AGE) + rank*COH_PROC
                newer_coh = (year-LAST_AGE) + (rank+1)*COH_PROC - 1

                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                
                if (RE_loop == 1 .and. rank <= NCOH_BATCH-1) then
                    Call ParallelCohorts_Myopic(older_coh, newer_coh, year, param, cut, typ, rsk_sup(:,:,it), &
                    epsAux, etaAux, PI_Myopic, &
                    ChoiceSimAux, ChoiceSimOut, &
                    ExperSimOut, WageSimOut, Cost1SimOut, &
                    Cost2SimOut, VSimOut, CohortWgtOut, &
                    empOut, sk_supOut) 
                else if (RE_loop == 2 .and. rank <= NCOH_BATCH-1) then    
                    if (year > FIRST_YR) then
                        rsk_lag = rsk_eq(:,:,year-1)
                    else
                        rsk_lag = 0.0
                    end if
                    Call ParallelCohorts_RE(older_coh, newer_coh, year, param, cut, typ, rsk_sup(:,:,it), rsk_lag, &
                    epsAux, etaAux, COEF_EXP, covar_EXP, PI_RE, &
                    ChoiceSimAux, ChoiceSimOut, &
                    ExperSimOut, WageSimOut, Cost1SimOut, &
                    Cost2SimOut, VSimOut, CohortWgtOut, &
                    empOut, sk_supOut)   
                end if
                
                if (rank /= 0 .and. rank <= NCOH_BATCH-1) then
                
                    tag = (rank+1)*10 + 1
                    count1 = size(ChoiceSimOut)
                    Call MPI_SEND(ChoiceSimOut,count1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
                    
                    tag = (rank+1)*10 + 2
                    count1 = size(ExperSimOut)
                    Call MPI_SEND(ExperSimOut,count1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
                    
                    tag = (rank+1)*10 + 3
                    count1 = size(WageSimOut)
                    Call MPI_SEND(WageSimOut,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)
                    
                    tag = (rank+1)*10 + 4
                    count1 = size(Cost1SimOut)
                    Call MPI_SEND(Cost1SimOut,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)
                    
                    tag = (rank+1)*10 + 5
                    count1 = size(Cost2SimOut)
                    Call MPI_SEND(Cost2SimOut,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)
                    
                    tag = (rank+1)*10 + 6
                    count1 = size(VSimOut)
                    Call MPI_SEND(VSimOut,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)
                    
                    tag = (rank+1)*10 + 7
                    count1 = size(empOut)
                    Call MPI_SEND(empOut,count1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
                    
                    tag = (rank+1)*10 + 8
                    count1 = size(sk_supOut)
                    Call MPI_SEND(sk_supOut,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)
                    
                    tag = (rank+1)*10 + 9
                    count1 = size(CohortWgtOut)
                    Call MPI_SEND(CohortWgtOut,count1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)
                
                end if
                
                if (rank == 0) then
                    
                    older_coh = (year-LAST_AGE) + rank*COH_PROC
                    newer_coh = (year-LAST_AGE) + (rank+1)*COH_PROC - 1

                    ChoiceSim(:,older_coh:newer_coh,year:year) = ChoiceSimOut
                    ExperSim(:,:,older_coh:newer_coh,year:year) = ExperSimOut
                    WageSim(:,older_coh:newer_coh,year:year) = WageSimOut
                    Cost1Sim(:,older_coh:newer_coh,year:year) = Cost1SimOut
                    Cost2Sim(:,older_coh:newer_coh,year:year) = Cost2SimOut
                    VSim(:,older_coh:newer_coh,year:year) = VSimOut
                    CohortWgt(older_coh:newer_coh,0:1) = CohortWgtOut
                    
                    sk_sup = sk_sup + sk_supOut
                    emp(:,:,it) = emp(:,:,it) + empOut
                    
                    do proc = 1, NCOH_BATCH-1
                    
                        older_coh = (year-LAST_AGE) + proc*COH_PROC
                        newer_coh = (year-LAST_AGE) + (proc+1)*COH_PROC - 1

                        tag = (proc+1)*10 + 1
                        count1 = size(ChoiceSimOut)
                        Call MPI_RECV(ChoiceSimOut,count1,MPI_INTEGER,proc,tag,MPI_COMM_WORLD,status,ierr)
                        ChoiceSim(:,older_coh:newer_coh,year:year) = ChoiceSimOut

                        tag = (proc+1)*10 + 2
                        count1 = size(ExperSimOut)
                        Call MPI_RECV(ExperSimOut,count1,MPI_INTEGER,proc,tag,MPI_COMM_WORLD,status,ierr)
                        ExperSim(:,:,older_coh:newer_coh,year:year) = ExperSimOut                        
                       
                        tag = (proc+1)*10 + 3
                        count1 = size(WageSimOut)
                        Call MPI_RECV(WageSimOut,count1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
                        WageSim(:,older_coh:newer_coh,year:year) = WageSimOut
                        
                        tag = (proc+1)*10 + 4
                        count1 = size(Cost1SimOut)
                        Call MPI_RECV(Cost1SimOut,count1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
                        Cost1Sim(:,older_coh:newer_coh,year:year) = Cost1SimOut
                    
                        tag = (proc+1)*10 + 5
                        count1 = size(Cost2SimOut)
                        Call MPI_RECV(Cost2SimOut,count1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
                        Cost2Sim(:,older_coh:newer_coh,year:year) = Cost2SimOut
                    
                        tag = (proc+1)*10 + 6
                        count1 = size(VSimOut)
                        Call MPI_RECV(VSimOut,count1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
                        VSim(:,older_coh:newer_coh,year:year) = VSimOut
                       
                        tag = (proc+1)*10 + 7
                        count1 = size(empOut)
                        Call MPI_RECV(empOut,count1,MPI_INTEGER,proc,tag,MPI_COMM_WORLD,status,ierr)
                        emp(:,:,it) = emp(:,:,it) + empOut
                        
                        tag = (proc+1)*10 + 8
                        count1 = size(sk_supOut)
                        Call MPI_RECV(sk_supOut,count1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
                        sk_sup = sk_sup + sk_supOut
                        
                        tag = (proc+1)*10 + 9
                        count1 = size(CohortWgtOut)
                        Call MPI_RECV(CohortWgtOut,count1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierr)
                        CohortWgt(older_coh:newer_coh,0:1) = CohortWgtOut
                    
                    end do

                end if
                

                if (rank == 0) then
                    
                    ! Start Newton's method to solve for K
                    do s = 1, NSECTORS
                        const1 = OutputData(s,t)*(1-alpha_prod(s,0,t)-alpha_prod(s,1,t))
                        const2 = alpha_prod(s,0,t)*sk_sup(s,0)**sigma_prod(s) + &
                                 alpha_prod(s,1,t)*sk_sup(s,1)**sigma_prod(s)
                        const3 = (1-alpha_prod(s,0,t)-alpha_prod(s,1,t))
                        const4 = rKData(t)
                        init = (1 - LaborSharesData(s,0,t) - LaborSharesData(s,1,t))*OutputData(s,t) / rKData(t)                
                        call NewtonSolver(const1,const2,const3,const4,sigma_prod(s),init,solution)
                        Capital(s) = solution
                        test1 = rKData(t)
                        test2 = ( OutputData(s,t) / ( alpha_prod(s,0,t)*sk_sup(s,0)**sigma_prod(s) + &
                                                      alpha_prod(s,1,t)*sk_sup(s,1)**sigma_prod(s) + & 
                                                     (1-alpha_prod(s,0,t)-alpha_prod(s,1,t))*Capital(s)**sigma_prod(s) ) )*(1-alpha_prod(s,0,t)-alpha_prod(s,1,t))*Capital(s)**(sigma_prod(s)-1)
                        if (abs(test1 - test2) > 0.001) then
                            print*, 'Warning on Newton Method'
                        end if
                    end do

                    do s = 1, NSECTORS
                        do ed = 0, 1
                            if (sk_sup(s,ed) > 0) then
                                rsk_dem(s,ed,it) = min(( alpha_prod(s,ed,t) / (1 - alpha_prod(s,0,t) - alpha_prod(s,1,t)) )*rKData(t)*((Capital(s)/sk_sup(s,ed))**(1-sigma_prod(s))),1000.0)
                            else
                                rsk_dem(s,ed,it) = min(2*rsk_sup(s,ed,it),1000.0)
                            end if
                        end do            
                    end do

        
                    check(:,:,t) = abs(log(rsk_sup(:,:,it)) - log(rsk_dem(:,:,it)))
            
        
                    converged = (check(1,0,t) <= crit .and. check(2,0,t) <= crit .and. check(3,0,t) <= crit .and. &
                                 check(4,0,t) <= crit .and. check(5,0,t) <= crit .and. check(6,0,t) <= crit .and. &
                                 check(7,0,t) <= crit .and. check(1,1,t) <= crit .and. check(2,1,t) <= crit .and. &
                                 check(3,1,t) <= crit .and. check(4,1,t) <= crit .and. check(5,1,t) <= crit .and. &
                                 check(6,1,t) <= crit .and. check(7,1,t) <= crit)
                
             
                    if ( .not. converged .and. it < MAXIT ) then
           
                        ! Adjust the step for updating the skill price
                        ! If rsk oscillates, decrease the step            
                        if (it < 30) then
                            step = 1.0/6.0
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
                        end if
           
                        ! Updating the return to skill
                        do s = 1, NSECTORS
                            do ed = 0, 1
                                rsk_sup(s,ed,it+1) = rsk_sup(s,ed,it) - step(s,ed)*(rsk_sup(s,ed,it) - rsk_dem(s,ed,it))    
                            end do
                        end do 
            
                    else 
        
                        ! Recording Equilibrium Quantities
        
                        rsk_eq(:,:,t)     = rsk_dem(:,:,it)     
                        sk_sup_eq(:,:,t)  = sk_sup  
                        CapitalSim(:,t)   = Capital
                        emp_eq(:,:,t)     = emp(:,:,it)            
                        iter(t)           = it
    
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
            
                do s = 1, NSECTORS
                    do ed = 0, 1
                        LaborSharesSim(s,ed,t) = (rsk_eq(s,ed,t)*sk_sup_eq(s,ed,t))/OutputData(s,t)
                    end do
                end do
            
                do s = 1, NSECTORS
                    CapitalSharesSim(s,t) = rKData(t)*CapitalSim(s,t) / OutputData(s,t)
                end do
            
            end if
            
        end do year_loop
        
        if (rank == 0) then
            
            ! Estimating the Stochastic Processes for HC prices
            dim1 = 9
            dim2 = 2
        
            do ed = 0,1
                do s = 1, NSECTORS
                    do i = 1, 9
                        VY(i,1) = log(rsk_eq(s,ed,FIRST_YR+1+i))-log(rsk_eq(s,ed,FIRST_YR+1+i-1))
                        VX(i,1) = 1.0
                        VX(i,2) = log(rsk_eq(s,ed,FIRST_YR+1+i-1))-log(rsk_eq(s,ed,FIRST_YR+1+i-2))
                    end do
                    Call LinReg(VY,VX,dim1,dim2,COEF,SST,SSE,info)  
                    COEF_EXP(s,ed,:,:) = COEF(:,:)
                    sigma_EXP(s,ed) = sqrt(SSE/7)
                    resid(s,ed,:) = VY(:,1) - COEF(1,1)*VX(:,1) - COEF(2,1)*VX(:,2)
                end do
            end do
        
            covar_EXP = 0.0
            do ed = 0,1
                do s1 = 1, NSECTORS
                    do s2 = 1, NSECTORS
                        do t = 1, 9
                            covar_EXP(s1,s2,ed) = covar_EXP(s1,s2,ed) + resid(s1,ed,t)*resid(s2,ed,t) / 7
                        end do
                    end do
                end do
            end do
        
            if (RE_loop == 1) then
                open(unit = 1, file = 'RSK_loop1.csv')
                write(1,987) 'Year', 'Educ', 'RSK_EQ1', 'RSK_EQ2', 'RSK_EQ3', 'RSK_EQ4', 'RSK_EQ5', 'RSK_EQ6', 'RSK_EQ7'
                do ed = 0, 1
                    do t = FIRST_YR, LAST_YR
                        write(1,781) t, ed, (rsk_eq(s,ed,t), s = 1, NSECTORS)
                    end do
                end do
                close(1)
                open(unit = 1, file = 'RSK_Regressions_loop1.csv')
                write(1,2121) 'Sector', 'Education', 'Intercept', 'Slope'
                do ed = 0,1
                    do s = 1, NSECTORS
                        write(1,2323) s, ed, COEF_EXP(s,ed,1,1), COEF_EXP(s,ed,2,1)
                    end do
                end do
                close(1)        
                open(unit = 1, file = 'RSK_COVAR0_loop1.csv')
                    do s1 = 1, NSECTORS
                        write(1,2929) (covar_EXP(s1,s2,0), s2 = 1, NSECTORS)
                    end do
                close(1)
                open(unit = 1, file = 'RSK_COVAR1_loop1.csv')
                    do s1 = 1, NSECTORS
                        write(1,2929) (covar_EXP(s1,s2,1), s2 = 1, NSECTORS)
                    end do
                close(1)
                elseif (RE_loop == 2) then
                    open(unit = 1, file = 'RSK_loop2.csv')
                    write(1,987) 'Year', 'Educ', 'RSK_EQ1', 'RSK_EQ2', 'RSK_EQ3', 'RSK_EQ4', 'RSK_EQ5', 'RSK_EQ6', 'RSK_EQ7'
                    do ed = 0, 1
                        do t = FIRST_YR, LAST_YR
                            write(1,781) t, ed, (rsk_eq(s,ed,t), s = 1, NSECTORS)
                        end do
                    end do
                    close(1)
                    open(unit = 1, file = 'RSK_Regressions_loop2.csv')
                    write(1,2121) 'Sector', 'Education', 'Intercept', 'Slope'
                    do ed = 0,1
                        do s = 1, NSECTORS
                            write(1,2323) s, ed, COEF_EXP(s,ed,1,1), COEF_EXP(s,ed,2,1)
                        end do
                    end do
                close(1)        
                open(unit = 1, file = 'RSK_COVAR0_loop2.csv')
                    do s1 = 1, NSECTORS
                        write(1,2929) (covar_EXP(s1,s2,0), s2 = 1, NSECTORS)
                    end do
                close(1)
                open(unit = 1, file = 'RSK_COVAR1_loop2.csv')
                    do s1 = 1, NSECTORS
                        write(1,2929) (covar_EXP(s1,s2,1), s2 = 1, NSECTORS)
                    end do
                close(1)
            elseif (RE_loop == 3) then
                open(unit = 1, file = 'RSK_loop3.csv')
                write(1,987) 'Year', 'Educ', 'RSK_EQ1', 'RSK_EQ2', 'RSK_EQ3', 'RSK_EQ4', 'RSK_EQ5', 'RSK_EQ6', 'RSK_EQ7'
                do ed = 0, 1
                    do t = FIRST_YR, LAST_YR
                        write(1,781) t, ed, (rsk_eq(s,ed,t), s = 1, NSECTORS)
                    end do
                end do
                close(1)
                open(unit = 1, file = 'RSK_Regressions_loop3.csv')
                write(1,2121) 'Sector', 'Education', 'Intercept', 'Slope'
                do ed = 0,1
                    do s = 1, NSECTORS
                        write(1,2323) s, ed, COEF_EXP(s,ed,1,1), COEF_EXP(s,ed,2,1)
                    end do
                end do
                close(1)        
                open(unit = 1, file = 'RSK_COVAR0_loop3.csv')
                    do s1 = 1, NSECTORS
                        write(1,2929) (covar_EXP(s1,s2,0), s2 = 1, NSECTORS)
                    end do
                close(1)
                open(unit = 1, file = 'RSK_COVAR1_loop3.csv')
                do s1 = 1, NSECTORS
                    write(1,2929) (covar_EXP(s1,s2,1), s2 = 1, NSECTORS)
                end do
                close(1)
            end if
        
            987  format(9(a14,','))
            781  format(2(i5,','),4(f20.8,','))        
            2121 format(4(a14,','))
            2323 format(2(i5,','),2(f20.8,','))
            2929 format(7(f20.8,','))  
             
        end if
        ! ********************************
        ! End of processing by processor 0
        ! ********************************
    
        count1 = size(COEF_EXP)
        Call MPI_BCAST(COEF_EXP,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  
        count2 = size(covar_EXP)
        Call MPI_BCAST(covar_EXP,count2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  
    
        RE_loop = RE_loop + 1

    end do RatExp_loop
    ! *************************************************
    ! If SUM_FLAG_EMAX > 0, the code jumps here
    ! *************************************************
    

    rankcond2: if (rank == 0) then
    
        if(SUM_FLAG_EMAX == 0) then  
        ! Otherwise just set F to 100000.000

            call date_and_time(values=time_array_0)
            finish_time_eq = time_array_0(5) * 3600 + time_array_0(6) * 60 &
                           + time_array_0(7) + 0.001 * time_array_0(8)

            print*
            print*, 'Time Computing Equilibrium:   ', real(finish_time_eq - start_time_eq,4)
            print*

            ! *****************************************************************************
            ! Construct array Dataset, that will be used in the auxiliary model regressions
            ! For the computation of the Loss Function
            ! *****************************************************************************

            Size2 = NPOP*(LAST_AGE-FIRST_AGE+1)*(LAST_YR-FIRST_YR+1)

            allocate(Dataset(Size2,15+NSECTORS))
            allocate(choice_index(Size2,0:NSECTORS))
            allocate(frac(0:NSECTORS,NPOP,FIRST_COH:LAST_COH))
            allocate(index_init(NPOP,FIRST_COH:LAST_COH,4))
            allocate(lag_choice_index(Size2,0:NSECTORS))
            allocate(switch_index(Size2,0:NSECTORS))
            allocate(index_wagedif(Size2,0:NSECTORS))

            obs             = 0
            sector_size     = 0
            lag_sector_size = 0
            ind             = 0
            ind2            = 0
            lag_ind         = 0
            switch_ind      = 0
            switch_size     = 0    
            wagedif_size    = 0
            size1998        = 0
            size2000        = 0
            size2005        = 0
            frac            = 0.0

            do t = FIRST_YR, LAST_YR
                do coh = (t-LAST_AGE), (t-FIRST_AGE)
                    do n = 1, NPOP
       
                        obs = obs + 1
            
                        educ = EducData(n,coh)
            
                        if (educ <= 2) then
                            ed = 0
                        else if (educ >= 3) then
                            ed = 1
                        end if
            
                        weight = sqrt(CohortWgt(coh,ed))
            
                        educ2 = 0
                        educ3 = 0
                        educ4 = 0
                        if (EducData(n,coh) == 2) then
                            educ2 = 1
                        elseif (EducData(n,coh) == 3) then
                            educ3 = 1
                        elseif (EducData(n,coh) == 4) then
                            educ4 = 1
                        end if
        
                        Dataset(obs,1)  = n
                        Dataset(obs,2)  = coh
                        Dataset(obs,3)  = t
                        Dataset(obs,4)  = GenderData(n,coh) - 1
                        Dataset(obs,5)  = educ2
                        Dataset(obs,6)  = educ3
                        Dataset(obs,7)  = educ4
                        Dataset(obs,8)  = (t-coh) - 25
                        Dataset(obs,9)  = ((t-coh) - 25)**2
                        Dataset(obs,10) = ExperSim(1,n,coh,t)
                        Dataset(obs,11) = ExperSim(2,n,coh,t)
                        Dataset(obs,12) = ExperSim(3,n,coh,t)
                        Dataset(obs,13) = ExperSim(4,n,coh,t)
                        Dataset(obs,14) = ExperSim(5,n,coh,t)
                        Dataset(obs,15) = ExperSim(6,n,coh,t)
                        Dataset(obs,16) = ExperSim(7,n,coh,t)
                        Dataset(obs,17) = ChoiceSim(n,coh,t)
                        ! lagged choice
                        Dataset(obs,18) = ChoiceSim(n,coh,t-1) 
                        ! lagged twice choice
                        Dataset(obs,19) = ChoiceSim(n,coh,t-2) 
                        Dataset(obs,20) = log(WageSim(n,coh,t))
                        ! lagged wage
                        if (t > FIRST_YR) then
                            Dataset(obs,21) = log(WageSim(n,coh,t-1))
                        else
                           Dataset(obs,21) = log(-1.0)
                        end if
                        Dataset(obs,22) = weight
            

                        if (t == 1995) then
                            index_init(n,coh,1) = obs
                        else if (t == 1998) then
                            index_init(n,coh,2) = obs
                            if ((t-coh) >= 28 .and. (t-coh) <= 60) then
                                size1998 = size1998 + 1
                            end if
                        else if (t == 2000) then
                            index_init(n,coh,3) = obs
                            if ((t-coh) >= 30 .and. (t-coh) <= 60) then
                                size2000 = size2000 + 1
                            end if
                        else if (t == 2005) then
                            index_init(n,coh,4) = obs
                            if ((t-coh) >= 35 .and. (t-coh) <= 60) then
                               size2005 = size2005 + 1
                            end if
                        end if
            
                        do s = 0, NSECTORS
                            if (int(ChoiceSim(n,coh,t)) == s) then
                                frac(s,n,coh)          = frac(s,n,coh) + 1.0
                                ind(s)                 = ind(s) + 1
                                choice_index(ind(s),s) = obs
                                sector_size(s)         = sector_size(s) + 1
                            end if
                            if (int(ChoiceSim(n,coh,t-1)) == s) then
                                lag_ind(s)                     = lag_ind(s) + 1
                                lag_choice_index(lag_ind(s),s) = obs
                                lag_sector_size(s)             = lag_sector_size(s) + 1
                            end if
                            if (int(ChoiceSim(n,coh,t-2)) == s .and. int(ChoiceSim(n,coh,t-1)) /= s) then
                                switch_ind(s)                 = switch_ind(s) + 1
                                switch_index(switch_ind(s),s) = obs
                                switch_size(s)                = switch_size(s) + 1
                            end if
                            if(t > FIRST_YR .and. int(ChoiceSim(n,coh,t)) == s .and. int(ChoiceSim(n,coh,t-1)) == s .and. (t-coh) >= 26) then
                                ind2(s) = ind2(s) + 1
                                index_wagedif(ind2(s),s) = obs
                                wagedif_size(s) = wagedif_size(s) + 1
                            end if
                        end do            
            
                    end do
                end do 
            end do


            do t = FIRST_YR, LAST_YR
                do s = 1, NSECTORS
                    z(s,t) = OutputData(s,t) / &
                            ( alpha_prod(s,0,t)*sk_sup_eq(s,0,t)**sigma_prod(s) + &
                              alpha_prod(s,1,t)*sk_sup_eq(s,1,t)**sigma_prod(s) + &
                             (1 - alpha_prod(s,0,t) - alpha_prod(s,1,t))*CapitalSim(s,t)**sigma_prod(s) )**(1/sigma_prod(s) )
                end do
            end do    
    
            ! Write Generated Data to File
            ! Call Write_Sim_Data(ChoiceSim,WageSim,ExperSim,typ,VSim,Cost1Sim,Cost2Sim,CohortWgt,CapitalSim,z,sk_sup_eq,rsk_eq)


            ! The equilibrium value for rsk will be used
            ! As the starting guess for the next iteration
            ! Of the optimization algorithm
            rsk_global(:,0) = min(max(rsk_eq(:,0,FIRST_YR),LOW_RSK(:,0)),UP_RSK(:,0))
            rsk_global(:,1) = min(max(rsk_eq(:,1,FIRST_YR),LOW_RSK(:,1)),UP_RSK(:,1))



            ! ****************************
            ! Compute aggregate statistics
            ! ****************************

            AggChoices = 0.0
            AggLogWages = 0.0
            AggTr = 0.0
            SumWeightChoices = 0.0

            do t = FIRST_YR, LAST_YR
                do coh = (t-LAST_AGE), (t-FIRST_AGE)
                    do n = 1, NPOP
                        if (EducData(n,coh) <= 2) then
                            ed = 0
                        else if (EducData(n,coh) >= 3) then
                            ed = 1
                        end if
                        age = t - coh
                        SumWeightChoices = SumWeightChoices + CohortWgt(coh,ed)

                        do s = 0, NSECTORS
                            if (ChoiceSim(n,coh,coh+age) == s) then
                                AggChoices(s) = AggChoices(s) + CohortWgt(coh,ed)
                                if (s > 0) then
                                    AggLogWages(s) = AggLogWages(s) + log(WageSim(n,coh,coh+age))*CohortWgt(coh,ed)
                                end if
                            end if
                        end do

                        do s1 = 0, NSECTORS
                            do s2 = 0, NSECTORS
                                if (ChoiceSim(n,coh,coh+age-1) == s1 .and. ChoiceSim(n,coh,coh+age) == s2) then
                                    AggTr(s1,s2) = AggTr(s1,s2) + CohortWgt(coh,ed)
                                end if
                            end do
                        end do

                    end do
                end do
            end do

            AggChoices2 = 0.0

            do s1 = 0, NSECTORS
                do s2 = 0, NSECTORS
                    AggChoices2(s1) = AggChoices2(s1) + AggTr(s1,s2)
                end do
            end do
    
            do s1 = 0, NSECTORS
                AggTr(s1,:) = AggTr(s1,:) / AggChoices2(s1)
            end do

            do s = 0, NSECTORS
                if (s > 0) then
                    AggLogWages(s) = AggLogWages(s) / AggChoices(s)
                end if
                AggChoices(s) = AggChoices(s) / SumWeightChoices
            end do

            AggChoices = 100.0*AggChoices
            AggTr = 100.0*AggTr

    
            F = 0.0

    
            ! Compute the Loss Function
            F = Loss_Function(Dataset,choice_index,lag_choice_index,sector_size,lag_sector_size, &
                              switch_index, switch_size, &
                              index_wagedif,wagedif_size,index_init,size1998,size2000,size2005,frac,&
                              LaborSharesSim,CapitalSharesSim)

            !**********************
            ! End of Iteration Time
            !**********************

            call date_and_time(values=time_array_1)
            finish_time = time_array_1(5) * 3600 + time_array_1(6) * 60 &
                        + time_array_1(7) + 0.001 * time_array_1(8) 
        

            do t = FIRST_YR, LAST_YR-1
                do s = 1, NSECTORS
                    do ed = 0, 1
                        rsk_ratio(s,ed,t) = rsk_eq(s,ed,t) / rsk_eq(s,ed,t+1)
                    end do
                end do
            end do

            print*, 'BATCH =', BATCH
            print*
            if (ALGORITHM == 1) then
                print*, 'ALGORITHM: NEWUOA'
            else if (ALGORITHM == 2) then
                print*, 'ALGORITHM: NELDER MEAD'
            end if
            print*
            print*, 'FUNCTION ITERATION: ', FUNCTION_ITER
            print*
            print*, 'theta:'
            write(*,90), (theta(i), i = 1, 7) 
            print*
            print*, 'wage equation 1'
            write(*,91), (beta(1,i), i = 1, 5)
            write(*,912), (beta(1,i), i = 6, 12)
            print*
            print*, 'wage equation 2'
            write(*,91), (beta(2,i), i = 1, 5)
            write(*,912), (beta(2,i), i = 6, 12)
            print*
            print*, 'wage equation 3'
            write(*,91), (beta(3,i), i = 1, 5)
            write(*,912), (beta(3,i), i = 6, 12)
            print*
            print*, 'wage equation 4'
            write(*,91), (beta(4,i), i = 1, 5)
            write(*,912), (beta(4,i), i = 6, 12)
            print*
            print*, 'wage equation 5'
            write(*,91), (beta(5,i), i = 1, 5)
            write(*,912), (beta(5,i), i = 6, 12)
            print*
            print*, 'wage equation 6'
            write(*,91), (beta(6,i), i = 1, 5)
            write(*,912), (beta(6,i), i = 6, 12)
            print*
            print*, 'wage equation 7'
            write(*,91), (beta(7,i), i = 1, 5)
            write(*,912), (beta(7,i), i = 6, 12)
            print*
            print*, 'sigma:'
            write(*,92), sigma
            print*  
            print*, 'cost:'
            write(*,93), kappa(1), kappa(2), kappa(3), kappa(4), kappa(5), kappa(6), kappa(7)
            write(*,93), 0.00     , kappa(8), kappa(9), kappa(10), kappa(11), kappa(12), kappa(13)
            write(*,93), kappa(14), kappa(15), kappa(16), kappa(17), kappa(18), kappa(19), kappa(20)
            write(*,94), kappa(21:26)            
            print*
            print*, 'lambda:'
            write(*,95), (lambda(i), i = 1, NTYPES)
            print*
            print*, 'tau:'
            write(*,96), (tau(i), i = 0, NSECTORS)
            print*
            print*, 'preference shock:'
            write(*,97), sigmaPref
            print*
            print*, 'omega 2:'
            write(*,98), (omega(2,i), i = 0, NSECTORS)
            print*
            print*, 'omega 3:'
            write(*,98), (omega(3,i), i = 0, NSECTORS)
            print*
            print*, 'gamma 2:'
            write(*,99), (gamma(2,i), i = 0, NSECTORS)
            write(*,991), (gamma(2,i), i = NSECTORS+1, NSECTORS+4)
            print*
            print*, 'gamma 3:'
            write(*,99), (gamma(3,i), i = 0, NSECTORS)
            write(*,991), (gamma(3,i), i = NSECTORS+1, NSECTORS+4)
            print*
            print*, 'sigma_prod:'
            write(*,900), sigma_prod
            print*
            print*, 'alpha_prod 1:'
            write(*,901), alpha_prod(1,0,FIRST_YR), alpha_prod(1,0,LAST_YR)
            write(*,901), alpha_prod(1,1,FIRST_YR), alpha_prod(1,1,LAST_YR)
            print*
            print*, 'alpha_prod 2:'
            write(*,901), alpha_prod(2,0,FIRST_YR), alpha_prod(2,0,LAST_YR)
            write(*,901), alpha_prod(2,1,FIRST_YR), alpha_prod(2,1,LAST_YR)
            print*
            print*, 'alpha_prod 3:'
            write(*,901), alpha_prod(3,0,FIRST_YR), alpha_prod(3,0,LAST_YR)
            write(*,901), alpha_prod(3,1,FIRST_YR), alpha_prod(3,1,LAST_YR)
            print*
            print*, 'alpha_prod 4:'
            write(*,901), alpha_prod(4,0,FIRST_YR), alpha_prod(4,0,LAST_YR)
            write(*,901), alpha_prod(4,1,FIRST_YR), alpha_prod(4,1,LAST_YR)
            print*
            print*, 'alpha_prod 5:'
            write(*,901), alpha_prod(5,0,FIRST_YR), alpha_prod(5,0,LAST_YR)
            write(*,901), alpha_prod(5,1,FIRST_YR), alpha_prod(5,1,LAST_YR)
            print*
            print*, 'alpha_prod 6:'
            write(*,901), alpha_prod(6,0,FIRST_YR), alpha_prod(6,0,LAST_YR)
            write(*,901), alpha_prod(6,1,FIRST_YR), alpha_prod(6,1,LAST_YR)
            print*
            print*, 'alpha_prod 7:'
            write(*,901), alpha_prod(7,0,FIRST_YR), alpha_prod(7,0,LAST_YR)
            write(*,901), alpha_prod(7,1,FIRST_YR), alpha_prod(7,1,LAST_YR)
            print*
            print*, 'rsk_init Ed = 0:'
            write(*,909) (rsk_init(s,0), s = 1, NSECTORS)
            print*
            print*, 'rsk_init Ed = 1:'
            write(*,909) (rsk_init(s,1), s = 1, NSECTORS)
            print*
            print*, 'Skill Prices Ed = 0:'
            do year = FIRST_YR, LAST_YR
                write(*,902) year, rsk_eq(:,0,year)
            end do
            print*
            print*, 'Skill Prices Ed = 1:'
            do year = FIRST_YR, LAST_YR
                write(*,902) year, (rsk_eq(s,1,year), s = 1, NSECTORS)
            end do
            print*
            print*, 'Iterations: '
            print*, iter
            print*
            print*, 'Check Ed = 0: '
            do year = FIRST_YR, LAST_YR
                write(*,903) year, (check(s,0,year), s = 1, NSECTORS)
            end do
            print*
            print*, 'Check Ed = 1: '
            do year = FIRST_YR, LAST_YR
                write(*,903) year, (check(s,1,year), s = 1, NSECTORS)
            end do
            print*
            print*, 'Min and Max R_sq for Emax:'
            if (RE_loop == 2) then
                write(*,904), minval(R_sq_Myopic(:,1:NSECTORS,:,:,:)), minloc(R_sq_Myopic(:,1:NSECTORS,:,:,:))
                write(*,904), maxval(R_sq_Myopic(:,1:NSECTORS,:,:,:)), maxloc(R_sq_Myopic(:,1:NSECTORS,:,:,:))
            else
                write(*,904), minval(R_sq_RE(:,1:NSECTORS,:,:,:)), minloc(R_sq_RE(:,1:NSECTORS,:,:,:))
                write(*,904), maxval(R_sq_RE(:,1:NSECTORS,:,:,:)), maxloc(R_sq_RE(:,1:NSECTORS,:,:,:))
            end if
            print*
            print*, 'Min and Max rsk_ratio:'
            write(*,906), minval(rsk_ratio), minloc(rsk_ratio)
            write(*,906), maxval(rsk_ratio), maxloc(rsk_ratio)
            print*
            print*, 'Sectoral Choices:'
            write(*,907), (AggChoices(s), s = 0, NSECTORS) 
            print*
            print*, 'Log Wages:'
            write(*,908), (AggLogWages(s), s = 1, NSECTORS)
            print*
            print*, 'Transitions:'
            write(*,907), (AggTr(0,s), s = 0, NSECTORS)
            write(*,907), (AggTr(1,s), s = 0, NSECTORS)
            write(*,907), (AggTr(2,s), s = 0, NSECTORS)
            write(*,907), (AggTr(3,s), s = 0, NSECTORS)
            write(*,907), (AggTr(4,s), s = 0, NSECTORS)
            write(*,907), (AggTr(5,s), s = 0, NSECTORS)
            write(*,907), (AggTr(6,s), s = 0, NSECTORS)
            write(*,907), (AggTr(7,s), s = 0, NSECTORS)
            print*  
            print*, 'Labor Shares 0 (2005):'
            write(*,908), (LaborSharesSim(s,0,LAST_YR), s = 1, NSECTORS)
            print*
            print*, 'Labor Shares 1 (2005):'
            write(*,908), (LaborSharesSim(s,1,LAST_YR), s = 1, NSECTORS)
            print*
            print*, 'Capital Shares (2005):'
            write(*,908), (CapitalSharesSim(s,LAST_YR), s = 1, NSECTORS)
            print*
            print*, '*****************'
            print*, 'Loss Function'
            print*, F
            print*, '*****************'
            print*
            print*, 'Iteration Time     :', real(finish_time - start_time,4)
            print*
            print* , '====================================================================='
            print* , '====================================================================='
            print*
            print*

            90  format(7(f7.4,'  '))
            91  format(5(f7.4,'  '))
            912 format(7(f7.4,'  '))    
            92  format(8(f5.2,'  '))
            93  format(7(f5.2,'  '))
            94  format(6(f7.4,'  '))
            95  format(3(f7.4,'  '))
            96  format(8(f7.4,'  '))
            97  format(f7.2,'  ')
            98  format(8(f7.4,'  '))
            99  format(8(f7.4,'  '))
            991 format(4(f7.4,'  '))    
            900 format(7(f5.2,'  '))
            901 format(7(f5.2,'  '))
            902 format(i5,'  ',7(f5.2,'  '))
            903 format(i5,'  ',7(f7.4,'  '))
            904 format(f5.2,'  ',5(i4))
            906 format(f5.2,'  ',i4,'  ',i4,'  ',i4,'  ')
            907 format(8(f5.2,'  '))
            908 format(7(f5.2,'  '))
            909 format(7(f6.3,'  '))
    
            deallocate(PI_Myopic,R_sq_Myopic,ChoiceSim,ExperSim,WageSim,eps,eta,PI_RE,R_sq_RE,Dataset, &
                       choice_index,frac,index_init,lag_choice_index,switch_index,index_wagedif)     

        elseif (SUM_FLAG_EMAX > 0) then
    
            F = 999999.9999

            print*, 'FUNCTION ITERATION: ', FUNCTION_ITER
            print*
            print*, '*****************'
            print*, 'Loss Function'
            print*, F
            print*, '*****************'
            print*
            print*, 'SUM_FLAG_EMAX > 0'
            print*
            print* , '====================================================================='
            print* , '====================================================================='
            print*
            print*

        end if  
        ! **********************************
        ! End of FLAG_EMAX Conditional
        ! **********************************

    end if rankcond2
    ! **********************************
    ! End of Processing over processor 0
    ! **********************************  
    
end if flagstop
! **********************
! End if(FLAG_STOP == 0)
! **********************

! Broadcast the Loss Function to all processors
Call MPI_BCAST(F,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)     
   
    
end subroutine SMM_OBJ_FCN








Function Loss_Function(Dataset,choice_index,lag_choice_index,sector_size,lag_sector_size, &
                       switch_index, switch_size, &
                       index_wagedif,wagedif_size,index_init,size1998,size2000,size2005,frac, &
                       LaborSharesSim,CapitalSharesSim)

! *****************************************************************************************
! Indirect Inference Loss Function
! *****************************************************************************************

USE Global_data
USE LinReg_MOD

implicit none

integer , intent(in) :: sector_size(0:), &
                        lag_sector_size(0:), &
                        choice_index(:,0:), &
                        lag_choice_index(:,0:), &
                        switch_index(:,0:), &
                        switch_size(0:), &
                        index_wagedif(:,0:), &
                        wagedif_size(0:), &
                        index_init(:,FIRST_COH:,:), &
                        size1998, &
                        size2000, &
                        size2005

real(KIND=DOUBLE), intent(in) ::  Dataset(:,:), frac(0:,:,FIRST_COH:), LaborSharesSim(:,0:,FIRST_YR:), &
                                  CapitalSharesSim(:,FIRST_YR:)

integer           i, j, n, t, age, coh, s, s1, s2, ed, &
                  index, info, obs, Size2, status_Wage(NSECTORS), &
                  status_Emp(0:NSECTORS), status_Tr(0:NSECTORS,0:NSECTORS), &
                  status_Ret(0:NSECTORS)

real(KIND=DOUBLE) betaSim(NSECTORS,NREG_Wage), &
                  sigmaWageSim(1:NSECTORS), &
                  gammaSim(0:NSECTORS,NREG_Emp), &
                  sigmaEmpSim(0:NSECTORS), &
                  phiSim(0:NSECTORS,0:NSECTORS,NREG_Tr), &
                  sigmaTrSim(0:NSECTORS,0:NSECTORS), &
                  rhoSim(0:NSECTORS,NREG_Return), &
                  sigmaRetSim(0:NSECTORS), &
                  SigmaWageDifSim(NSECTORS), &
                  xsi1998Sim(0:NSECTORS,NREG_Pers), &
                  xsi2000Sim(0:NSECTORS,NREG_Pers), &
                  xsi2005Sim(0:NSECTORS,NREG_Pers), &
                  etaSim(0:NSECTORS,NREG_Freq), &
                  SSE, SST

real(KIND=DOUBLE) Loss_Function, &
                  Loss_Function_Wage, &
                  Loss_Function_Emp, & 
                  Loss_Function_Tr, &
                  Loss_Function_Return, &
                  Loss_Function_Sigma, &
                  Loss_Function_SigmaDif, &
                  Loss_Function_Pers1998, &
                  Loss_Function_Pers2000, &
                  Loss_Function_Pers2005, &
                  Loss_Function_Freq, &
                  Loss_Function_Labor, &
                  Loss_Function_Capital, &
                  A(NREG_Wage,NREG_Wage), & 
                  B(NREG_Emp,NREG_Emp), &
                  C(NREG_Tr,NREG_Tr), &
                  D(NREG_Pers,NREG_Pers), &
                  E(NREG_Freq,NREG_Freq), &
                  F(NREG_Return,NREG_Return), &
                  wgt


! Elapsed time variables
integer           time_array_0(8), time_array_1(8)
real(KIND=DOUBLE) start_time, finish_time 

real(KIND=DOUBLE) , allocatable, dimension(:,:) :: XW, XXW, invXXW, Wage, coefw, eps, &
                                                   XEmp, XXEmp, invXXEmp, Dummy, coefe, &
                                                   XTr, XXTr, invXXTr, coeft, WageDif, coefwd, &
                                                   XRet, XXRet, invXXRet, DummyRet, coefret, &
                                                   Xpers1998, XXPers1998, invXXPers1998, coefpers1998, &
                                                   Xpers2000, XXPers2000, invXXPers2000, coefpers2000, &
                                                   Xpers2005, XXPers2005, invXXPers2005, coefpers2005, &
                                                   Y, XF, XXF, invXXF, coeffreq
                                                   
! *****************************************************************************************                                                   
! *****************************************************************************************

status_Wage = 0
status_Emp  = 0
status_Tr   = 0

Size2 = NPOP*(LAST_AGE-FIRST_AGE+1)*(LAST_YR-FIRST_YR+1)


! ****************
! Wage Regressions
! ****************

allocate(coefw(NREG_Wage,NSECTORS))
        
do s = 1, NSECTORS
        
    allocate(Wage(sector_size(s),1))
    allocate(eps(sector_size(s),1))
    allocate(XW(sector_size(s),NREG_Wage))
    allocate(XXW(NREG_Wage,NREG_Wage))
    allocate(invXXW(NREG_Wage,NREG_Wage))

    do t = FIRST_YR, LAST_YR
        do i = 1, sector_size(s)
            if (int(Dataset(choice_index(i,s),3)) == t) then
                XW(i,t-FIRST_YR+1) = 1
            else
                XW(i,t-FIRST_YR+1) = 0
            end if       
        end do
    end do
    
    XW(:,12:NREG_Wage) = Dataset(choice_index(1:sector_size(s),s),4:16)

    do i = 1, NREG_Wage       
        XW(:,i)   = XW(:,i)*Dataset(choice_index(1:sector_size(s),s),22)        
    end do

    Wage(:,1) = Dataset(choice_index(1:sector_size(s),s),20)
    Wage(:,1) = Wage(:,1)*Dataset(choice_index(1:sector_size(s),s),22) 

    call get_invXX(XW,sector_size(s),NREG_Wage,XXW,invXXW,info)
    
    status_Wage(s) = info
    
    call LinReg2(invXXW,XW,Wage,sector_size(s),NREG_Wage,coefw(:,s:s),SSE,SST,1,eps) 
    
    betaSim(s,:) = coefw(:,s)
    
    sigmaWageSim(s)  = sqrt(SSE / sum(Dataset(choice_index(1:sector_size(s),s),22)**2))
    
    deallocate(Wage,eps,XW,XXW,invXXW)

end do       

deallocate(coefw)

! ******************************************************************************
! ******************************************************************************



! ******************************************************
! Employment Regressions
! ******************************************************
 
allocate(XEmp(Size2,NREG_Emp))
allocate(XXEmp(NREG_Emp,NREG_Emp))
allocate(invXXEmp(NREG_Emp,NREG_Emp))
allocate(Dummy(Size2,0:NSECTORS))
allocate(coefe(NREG_Emp,0:NSECTORS))


do t = FIRST_YR, LAST_YR
    where (int(Dataset(:,3)) == t)
        XEmp(:,t-FIRST_YR+1) = 1
    elsewhere
        XEmp(:,t-FIRST_YR+1) = 0
    end where       
end do

XEmp(:,12:NREG_Emp) = Dataset(:,4:16)


! *********
! Weighting
! *********

do i = 1, NREG_Emp       
    XEmp(:,i)   = XEmp(:,i)*Dataset(:,22)        
end do


call get_invXX(XEmp,Size2,NREG_Emp,XXEmp,invXXEmp,info)

status_Emp = info

do s = 0, NSECTORS
    
    where (int(Dataset(:,17)) == s)
        Dummy(:,s) = 1
    elsewhere
        Dummy(:,s) = 0
    end where 

    Dummy(:,s) = Dummy(:,s)*Dataset(:,22)
    
end do    

do s = 0, NSECTORS
    call LinReg2(invXXEmp,XEmp,Dummy(:,s:s),Size2,NREG_Emp,coefe(:,s:s),SSE,SST,0)
end do        

do s = 0, NSECTORS
    gammaSim(s,:) = coefe(:,s)
end do

deallocate(XEmp,XXEmp,invXXEmp,Dummy,coefe)



! ******************************************************************************
! ******************************************************************************



! ******************************************************
! Transition Regressions
! ******************************************************

do s1 = 0, NSECTORS

    allocate(XTr(lag_sector_size(s1),NREG_Tr))  
    allocate(XXTr(NREG_Tr,NREG_Tr)) 
    allocate(invXXTr(NREG_Tr,NREG_Tr))
    allocate(Dummy(lag_sector_size(s1),0:NSECTORS)) 
    allocate(coeft(NREG_Tr,0:NSECTORS))
   
    do t = FIRST_YR, LAST_YR
        do i = 1, lag_sector_size(s1)
            if (int(Dataset(lag_choice_index(i,s1),3)) == t) then
                XTr(i,t-FIRST_YR+1) = 1
            else
                XTr(i,t-FIRST_YR+1) = 0
            end if
        end do       
    end do      
    
    XTr(:,12:NREG_Tr) = Dataset(lag_choice_index(1:lag_sector_size(s1),s1),4:16)

    do i = 1, NREG_Tr
        XTr(:,i) = XTr(:,i)*Dataset(lag_choice_index(1:lag_sector_size(s1),s1),22)
    end do

    call get_invXX(XTr,lag_sector_size(s1),NREG_Tr,XXTr,invXXTr,info)

    status_Tr(s1,:) = info

    do s2 = 0, NSECTORS

        do i = 1, lag_sector_size(s1)
            if (int(Dataset(lag_choice_index(i,s1),17)) == s2) then
                Dummy(i,s2) = 1
            else
                Dummy(i,s2) = 0
            end if
        end do  

        Dummy(:,s2) = Dummy(:,s2)*Dataset(lag_choice_index(1:lag_sector_size(s1),s1),22)
        
        call LinReg2(invXXTr,XTr,Dummy(:,s2:s2),lag_sector_size(s1),NREG_Tr,coeft(:,s2:s2),SSE,SST,0)
    
    end do    
    
    do s2 = 0, NSECTORS
        phiSim(s1,s2,:) = coeft(:,s2)
    end do
    
    deallocate(XTr,XXTr,invXXTr,Dummy,coeft)      
     
end do    

! ******************************************************************************
! ******************************************************************************



! ******************************************************
! Return
! ******************************************************

do s = 0, NSECTORS

    allocate(XRet(switch_size(s),NREG_Return))  
    allocate(XXRet(NREG_Return,NREG_Return)) 
    allocate(invXXRet(NREG_Return,NREG_Return))
    allocate(DummyRet(switch_size(s),1)) 
    allocate(coefret(NREG_Return,1))
   
    do t = FIRST_YR, LAST_YR
        do i = 1, switch_size(s)
            if (int(Dataset(switch_index(i,s),3)) == t) then
                XRet(i,t-FIRST_YR+1) = 1
            else
                XRet(i,t-FIRST_YR+1) = 0
            end if
        end do       
    end do      
    
    XRet(:,12:NREG_Return) = Dataset(switch_index(1:switch_size(s),s),4:16)

    do i = 1, NREG_Return
        XRet(:,i) = XRet(:,i)*Dataset(switch_index(1:switch_size(s),s),22)
    end do

    call get_invXX(XRet,switch_size(s),NREG_Return,XXRet,invXXRet,info)

    status_Ret(s) = info

    do i = 1, switch_size(s)
        if (int(Dataset(switch_index(i,s),17)) == s) then
            DummyRet(i,1) = 1
        else
            DummyRet(i,1) = 0
        end if
    end do

    DummyRet(:,1) = DummyRet(:,1)*Dataset(switch_index(1:switch_size(s),s),22)

    call LinReg2(invXXRet,XRet,DummyRet(:,1:1),switch_size(s),NREG_Return,coefret(:,1:1),SSE,SST,0)
    
    rhoSim(s,:) = coefret(:,1)
    
    deallocate(XRet,XXRet,invXXRet,DummyRet,coefret)      
     
end do    

! ******************************************************************************
! ******************************************************************************




! **************************************
! 1st Difference in Log Wage Regressions
! **************************************

allocate(coefwd(NREG_WageDif,NSECTORS))

do s = 1, NSECTORS

    allocate(WageDif(wagedif_size(s),1))
    allocate(XW(wagedif_size(s),NREG_WageDif))
    allocate(XXW(NREG_WageDif,NREG_WageDif))
    allocate(invXXW(NREG_WageDif,NREG_WageDif))

    do t = FIRST_YR + 1, LAST_YR
        do i = 1, wagedif_size(s)
            if (int(Dataset(index_wagedif(i,s),3)) == t) then
                XW(i,t-FIRST_YR) = 1
            else
                XW(i,t-FIRST_YR) = 0
            end if       
        end do
    end do

    XW(:,11) = Dataset(index_wagedif(1:wagedif_size(s),s),8)

    do i = 1, NREG_WageDif       
        XW(:,i)   = XW(:,i)*Dataset(index_wagedif(1:wagedif_size(s),s),22)        
    end do

    WageDif(:,1) = Dataset(index_wagedif(1:wagedif_size(s),s),20) - Dataset(index_wagedif(1:wagedif_size(s),s),21)
    WageDif(:,1) = WageDif(:,1)*Dataset(index_wagedif(1:wagedif_size(s),s),22) 

    call get_invXX(XW,wagedif_size(s),NREG_WageDif,XXW,invXXW,info)
    
    call LinReg2(invXXW,XW,WageDif,wagedif_size(s),NREG_WageDif,coefwd(:,s:s),SSE,SST,1) 
    
    SigmaWageDifSim(s) = sqrt(SSE / sum(Dataset(index_wagedif(1:wagedif_size(s),s),22)**2))
    
    deallocate(WageDif,XW,XXW,invXXW)

end do

deallocate(coefwd)

! ******************************************************************************
! ******************************************************************************


! ******************************************************
! Persistence Regressions 
! ******************************************************

! ****
! 1998
! ****
 
allocate(XPers1998(size1998,NREG_Pers))
allocate(XXPers1998(NREG_Pers,NREG_Pers))
allocate(invXXPers1998(NREG_Pers,NREG_Pers))
allocate(Dummy(size1998,0:NSECTORS))
allocate(coefpers1998(NREG_Pers,0:NSECTORS))

Dummy     = 0.0
XPers1998 = 0.0
i         = 0        
        
        
do coh = FIRST_YR - 57, FIRST_YR - 25
    do n = 1, NPOP
        obs            = index_init(n,coh,1)
        wgt            = Dataset(obs,22)
        i              = i + 1
        if (Dataset(index_init(n,coh,2),17) == 0) then
            Dummy(i,0) = 1.0*wgt
        else if (Dataset(index_init(n,coh,2),17) == 1) then
            Dummy(i,1) = 1.0*wgt
        else if (Dataset(index_init(n,coh,2),17) == 2) then
            Dummy(i,2) = 1.0*wgt
        else if (Dataset(index_init(n,coh,2),17) == 3) then
            Dummy(i,3) = 1.0*wgt
        else if (Dataset(index_init(n,coh,2),17) == 4) then
            Dummy(i,4) = 1.0*wgt
        else if (Dataset(index_init(n,coh,2),17) == 5) then
            Dummy(i,5) = 1.0*wgt
        else if (Dataset(index_init(n,coh,2),17) == 6) then
            Dummy(i,6) = 1.0*wgt
        else if (Dataset(index_init(n,coh,2),17) == 7) then
            Dummy(i,7) = 1.0*wgt
        end if
        
        
        if (Dataset(obs,18) == 0) then      
            XPers1998(i,1)         = 1.0*wgt
        else if (Dataset(obs,18) == 1) then
            XPers1998(i,2)         = 1.0*wgt
        else if (Dataset(obs,18) == 2) then
            XPers1998(i,3)         = 1.0*wgt
        else if (Dataset(obs,18) == 3) then
            XPers1998(i,4)         = 1.0*wgt
        else if (Dataset(obs,18) == 4) then
            XPers1998(i,5)         = 1.0*wgt
        else if (Dataset(obs,18) == 5) then
            XPers1998(i,6)         = 1.0*wgt
        else if (Dataset(obs,18) == 6) then
            XPers1998(i,7)         = 1.0*wgt
        else if (Dataset(obs,18) == 7) then
            XPers1998(i,8)         = 1.0*wgt
        end if
        XPers1998(i,9)         = Dataset(obs,4)*wgt
        XPers1998(i,10)        = Dataset(obs,5)*wgt
        XPers1998(i,11)        = Dataset(obs,6)*wgt
        XPers1998(i,12)        = Dataset(obs,7)*wgt
        XPers1998(i,13)        = Dataset(obs,8)*wgt
        XPers1998(i,14)        = Dataset(obs,9)*wgt
        XPers1998(i,15)        = Dataset(obs,10)*wgt
        XPers1998(i,16)        = Dataset(obs,11)*wgt
        XPers1998(i,17)        = Dataset(obs,12)*wgt
        XPers1998(i,18)        = Dataset(obs,13)*wgt
        XPers1998(i,19)        = Dataset(obs,14)*wgt
        XPers1998(i,20)        = Dataset(obs,15)*wgt
        XPers1998(i,21)        = Dataset(obs,16)*wgt
    end do  
end do

call get_invXX(XPers1998,size1998,NREG_Pers,XXPers1998,invXXPers1998,info)


do s = 0, NSECTORS
    call LinReg2(invXXPers1998,XPers1998,Dummy(:,s:s),size1998,NREG_Pers,coefpers1998(:,s:s),SSE,SST,0)
end do


do s = 0, NSECTORS
    xsi1998Sim(s,:) = coefpers1998(:,s)
end do

    
deallocate(XPers1998,XXPers1998,invXXPers1998,Dummy,coefpers1998)



! ****
! 2000
! ****
 
allocate(XPers2000(size2000,NREG_Pers))
allocate(XXPers2000(NREG_Pers,NREG_Pers))
allocate(invXXPers2000(NREG_Pers,NREG_Pers))
allocate(Dummy(size2000,0:NSECTORS))
allocate(coefpers2000(NREG_Pers,0:NSECTORS))

Dummy     = 0.0
XPers2000 = 0.0
i         = 0        
        
        
do coh = FIRST_YR - 55, FIRST_YR - 25
    do n = 1, NPOP
        obs            = index_init(n,coh,1)
        wgt            = Dataset(obs,22)
        i              = i + 1
        if (Dataset(index_init(n,coh,3),17) == 0) then
            Dummy(i,0) = 1.0*wgt
        else if (Dataset(index_init(n,coh,3),17) == 1) then
            Dummy(i,1) = 1.0*wgt
        else if (Dataset(index_init(n,coh,3),17) == 2) then
            Dummy(i,2) = 1.0*wgt
        else if (Dataset(index_init(n,coh,3),17) == 3) then
            Dummy(i,3) = 1.0*wgt
        else if (Dataset(index_init(n,coh,3),17) == 4) then
            Dummy(i,4) = 1.0*wgt
        else if (Dataset(index_init(n,coh,3),17) == 5) then
            Dummy(i,5) = 1.0*wgt
        else if (Dataset(index_init(n,coh,3),17) == 6) then
            Dummy(i,6) = 1.0*wgt
        else if (Dataset(index_init(n,coh,3),17) == 7) then
            Dummy(i,7) = 1.0*wgt
        end if
        if (Dataset(obs,18) == 0) then      
            XPers2000(i,1)         = 1.0*wgt
        else if (Dataset(obs,18) == 1) then
            XPers2000(i,2)         = 1.0*wgt
        else if (Dataset(obs,18) == 2) then
            XPers2000(i,3)         = 1.0*wgt
        else if (Dataset(obs,18) == 3) then
            XPers2000(i,4)         = 1.0*wgt
        else if (Dataset(obs,18) == 4) then
            XPers2000(i,5)         = 1.0*wgt
        else if (Dataset(obs,18) == 5) then
            XPers2000(i,6)         = 1.0*wgt
        else if (Dataset(obs,18) == 6) then
            XPers2000(i,7)         = 1.0*wgt
        else if (Dataset(obs,18) == 7) then
            XPers2000(i,8)         = 1.0*wgt
        end if
        XPers2000(i,9)         = Dataset(obs,4)*wgt
        XPers2000(i,10)        = Dataset(obs,5)*wgt
        XPers2000(i,11)        = Dataset(obs,6)*wgt
        XPers2000(i,12)        = Dataset(obs,7)*wgt
        XPers2000(i,13)        = Dataset(obs,8)*wgt
        XPers2000(i,14)        = Dataset(obs,9)*wgt
        XPers2000(i,15)        = Dataset(obs,10)*wgt
        XPers2000(i,16)        = Dataset(obs,11)*wgt
        XPers2000(i,17)        = Dataset(obs,12)*wgt
        XPers2000(i,18)        = Dataset(obs,13)*wgt
        XPers2000(i,19)        = Dataset(obs,14)*wgt
        XPers2000(i,20)        = Dataset(obs,15)*wgt
        XPers2000(i,21)        = Dataset(obs,16)*wgt
    end do  
end do

call get_invXX(XPers2000,size2000,NREG_Pers,XXPers2000,invXXPers2000,info)


do s = 0, NSECTORS
    call LinReg2(invXXPers2000,XPers2000,Dummy(:,s:s),size2000,NREG_Pers,coefpers2000(:,s:s),SSE,SST,0)
end do    

do s = 0, NSECTORS
    xsi2000Sim(s,:) = coefpers2000(:,s)
end do

deallocate(XPers2000,XXPers2000,invXXPers2000,Dummy,coefpers2000)



! ****
! 2005
! ****
 
allocate(XPers2005(size2005,NREG_Pers))
allocate(XXPers2005(NREG_Pers,NREG_Pers))
allocate(invXXPers2005(NREG_Pers,NREG_Pers))
allocate(Dummy(size2005,0:NSECTORS))
allocate(coefpers2005(NREG_Pers,0:NSECTORS))

Dummy     = 0.0
XPers2005 = 0.0
i         = 0        
        
        
do coh = FIRST_YR - 50, FIRST_YR - 25
    do n = 1, NPOP
        obs            = index_init(n,coh,1)
        wgt            = Dataset(obs,22)
        i              = i + 1
        if (Dataset(index_init(n,coh,4),17) == 0) then
            Dummy(i,0) = 1.0*wgt
        else if (Dataset(index_init(n,coh,4),17) == 1) then
            Dummy(i,1) = 1.0*wgt
        else if (Dataset(index_init(n,coh,4),17) == 2) then
            Dummy(i,2) = 1.0*wgt
        else if (Dataset(index_init(n,coh,4),17) == 3) then
            Dummy(i,3) = 1.0*wgt
        else if (Dataset(index_init(n,coh,4),17) == 4) then
            Dummy(i,4) = 1.0*wgt
        else if (Dataset(index_init(n,coh,4),17) == 5) then
            Dummy(i,5) = 1.0*wgt
        else if (Dataset(index_init(n,coh,4),17) == 6) then
            Dummy(i,6) = 1.0*wgt
        else if (Dataset(index_init(n,coh,4),17) == 7) then
            Dummy(i,7) = 1.0*wgt            
        end if
        if (Dataset(obs,18) == 0) then      
            XPers2005(i,1)         = 1.0*wgt
        else if (Dataset(obs,18) == 1) then
            XPers2005(i,2)         = 1.0*wgt
        else if (Dataset(obs,18) == 2) then
            XPers2005(i,3)         = 1.0*wgt
        else if (Dataset(obs,18) == 3) then
            XPers2005(i,4)         = 1.0*wgt
        else if (Dataset(obs,18) == 4) then
            XPers2005(i,5)         = 1.0*wgt
        else if (Dataset(obs,18) == 5) then
            XPers2005(i,6)         = 1.0*wgt
        else if (Dataset(obs,18) == 6) then
            XPers2005(i,7)         = 1.0*wgt
        else if (Dataset(obs,18) == 7) then
            XPers2005(i,8)         = 1.0*wgt
        end if
        XPers2005(i,9)         = Dataset(obs,4)*wgt
        XPers2005(i,10)        = Dataset(obs,5)*wgt
        XPers2005(i,11)        = Dataset(obs,6)*wgt
        XPers2005(i,12)        = Dataset(obs,7)*wgt
        XPers2005(i,13)        = Dataset(obs,8)*wgt
        XPers2005(i,14)        = Dataset(obs,9)*wgt
        XPers2005(i,15)        = Dataset(obs,10)*wgt
        XPers2005(i,16)        = Dataset(obs,11)*wgt
        XPers2005(i,17)        = Dataset(obs,12)*wgt
        XPers2005(i,18)        = Dataset(obs,13)*wgt
        XPers2005(i,19)        = Dataset(obs,14)*wgt
        XPers2005(i,20)        = Dataset(obs,15)*wgt
        XPers2005(i,21)        = Dataset(obs,16)*wgt
    end do  
end do

call get_invXX(XPers2005,size2005,NREG_Pers,XXPers2005,invXXPers2005,info)

do s = 0, NSECTORS
    call LinReg2(invXXPers2005,XPers2005,Dummy(:,s:s),size2005,NREG_Pers,coefpers2005(:,s:s),SSE,SST,0)
end do

do s = 0, NSECTORS
    xsi2005Sim(s,:) = coefpers2005(:,s)
end do

deallocate(XPers2005,XXPers2005,invXXPers2005,Dummy,coefpers2005)

! ******************************************************************************
! ******************************************************************************


! *************************
! Frequency Regressions
! *************************

allocate(Y(NPOP*(50-25+1),0:NSECTORS))
allocate(XF(NPOP*(50-25+1),NREG_Freq))
allocate(XXF(NREG_Freq,NREG_Freq))
allocate(invXXF(NREG_Freq,NREG_Freq))
allocate(coeffreq(NREG_Freq,0:7))

XF = 0.0
i  = 0

do coh = FIRST_YR - 50, FIRST_YR - 25
    do n = 1, NPOP
        obs             = index_init(n,coh,1)
        i               = i + 1
        wgt             = Dataset(obs,22)
        do s = 0, NSECTORS
            Y(i,s) = frac(s,n,coh)*wgt 
        end do
        if (Dataset(obs,18) == 0) then      
            XF(i,1)         = 1.0*wgt
        else if (Dataset(obs,18) == 1) then
            XF(i,2)         = 1.0*wgt
        else if (Dataset(obs,18) == 2) then
            XF(i,3)         = 1.0*wgt
        else if (Dataset(obs,18) == 3) then
            XF(i,4)         = 1.0*wgt
        else if (Dataset(obs,18) == 4) then
            XF(i,5)         = 1.0*wgt
        else if (Dataset(obs,18) == 5) then
            XF(i,6)         = 1.0*wgt
        else if (Dataset(obs,18) == 6) then
            XF(i,7)         = 1.0*wgt
        else if (Dataset(obs,18) == 7) then
            XF(i,8)         = 1.0*wgt
        end if
        XF(i,9)         = Dataset(obs,4)*wgt
        XF(i,10)        = Dataset(obs,5)*wgt
        XF(i,11)        = Dataset(obs,6)*wgt
        XF(i,12)        = Dataset(obs,7)*wgt
        XF(i,13)        = Dataset(obs,8)*wgt
        XF(i,14)        = Dataset(obs,9)*wgt
        XF(i,15)        = Dataset(obs,10)*wgt
        XF(i,16)        = Dataset(obs,11)*wgt
        XF(i,17)        = Dataset(obs,12)*wgt
        XF(i,18)        = Dataset(obs,13)*wgt
        XF(i,19)        = Dataset(obs,14)*wgt
        XF(i,20)        = Dataset(obs,15)*wgt
        XF(i,21)        = Dataset(obs,16)*wgt
    end do  
end do

call get_invXX(XF,NPOP*(50-25+1),NREG_Freq,XXF,invXXF,info)

Size2 = NPOP*(50-25+1)

do s = 0, NSECTORS
    call LinReg2(invXXF,XF,Y(:,s:s),Size2,NREG_Freq,coeffreq(:,s:s),SSE,SST,0)
end do

do s = 0, NSECTORS
    etaSim(s,:) = coeffreq(:,s)
end do

! ******************************************************************************
! ******************************************************************************


! ****************************
! We compute the Loss Function
! ****************************

Loss_Function_Wage     = 0.0
Loss_Function_Emp      = 0.0
Loss_Function_Tr       = 0.0
Loss_Function_Return   = 0.0
Loss_Function_Sigma    = 0.0
Loss_Function_SigmaDif = 0.0
Loss_Function_Pers1998 = 0.0
Loss_Function_Pers2000 = 0.0
Loss_Function_Pers2005 = 0.0
Loss_Function_Freq     = 0.0
Loss_Function_Labor    = 0.0
Loss_Function_Capital  = 0.0


do s = 1, NSECTORS

    A = CPWageData(s,:,:)
    
    call dpotrf('L', NREG_Wage, A, NREG_Wage, info)
    if (info .ne. 0) then
        print*, 'CPWageData is not symmetric positive definite.'
        print*, 'Aborting program.'
        pause
        stop
    end if


    do i = 1, NREG_Wage
        do j = 1, NREG_Wage
            Loss_Function_Wage = Loss_Function_Wage + &
            (betaSim(s,i)-betaData(s,i))* & 
            invCOVbetaData(s,i,j)*(betaSim(s,j)-betaData(s,j))
        end do
    end do

end do  


do s = 0, NSECTORS

    B = CPEmpData(s,:,:)
    
    call dpotrf('L', NREG_Emp, B, NREG_Emp, info)
    if (info .ne. 0) then
        print*, 'CPEmpData is not symmetric positive definite.'
        print*, 'Aborting program.'
        pause
        stop
    end if


    do i = 1, NREG_Emp    
        do j = 1, NREG_Emp
            Loss_Function_Emp = Loss_Function_Emp + &
            (gammaSim(s,i)-gammaData(s,i))* & 
            invCOVgammaData(s,i,j)*(gammaSim(s,j)-gammaData(s,j))
        end do
    end do
    
    
end do   


do s1 = 0, NSECTORS
    do s2 = 0, NSECTORS
    
        C = CPTrData(s1,s2,:,:)
    
        call dpotrf('L', NREG_Tr, C, NREG_Tr, info)
        if (info .ne. 0) then
            print*, 'CPTrData is not symmetric positive definite.'
            print*, 'Aborting program.'
            pause
            stop
        end if
    


        do i = 1, NREG_Tr
            do j = 1, NREG_Tr
                Loss_Function_Tr = Loss_Function_Tr + &
                (phiSim(s1,s2,i)-phiData(s1,s2,i))* & 
                invCOVphiData(s1,s2,i,j)*(phiSim(s1,s2,j)-phiData(s1,s2,j))
            end do
        end do
        
    end do    
end do    


do s = 1, NSECTORS

    F = CPReturnData(s,:,:)
    
    call dpotrf('L', NREG_Return, F, NREG_Return, info)
    if (info .ne. 0) then
        print*, 'CPReturnData is not symmetric positive definite.'
        print*, 'Aborting program.'
        pause
        stop
    end if


    do i = 1, NREG_Return
        do j = 1, NREG_Return
            Loss_Function_Return = Loss_Function_Return + &
            (rhoSim(s,i)-rhoData(s,i))* & 
            invCOVrhoData(s,i,j)*(rhoSim(s,j)-rhoData(s,j))
        end do
    end do

end do  


do s = 0, NSECTORS

    D = CPPers1998Data(s,:,:)
    
    call dpotrf('L', NREG_Pers, D, NREG_Pers, info)
    if (info .ne. 0) then
        print*, 'CPPers1998Data is not symmetric positive definite.'
        print*, 'Aborting program.'
        pause
        stop
    end if


    do i = 1, NREG_Pers    
        do j = 1, NREG_Pers
            Loss_Function_Pers1998 = Loss_Function_Pers1998 + &
            (xsi1998Sim(s,i)-xsi1998Data(s,i))* & 
            invCOVxsi1998Data(s,i,j)*(xsi1998Sim(s,j)-xsi1998Data(s,j))
        end do
    end do
    
end do



do s = 0, NSECTORS

    D = CPPers2000Data(s,:,:)
    
    call dpotrf('L', NREG_Pers, D, NREG_Pers, info)
    if (info .ne. 0) then
        print*, 'CPPers2000Data is not symmetric positive definite.'
        print*, 'Aborting program.'
        pause
        stop
    end if


    do i = 1, NREG_Pers    
        do j = 1, NREG_Pers
            Loss_Function_Pers2000 = Loss_Function_Pers2000 + &
            (xsi2000Sim(s,i)-xsi2000Data(s,i))* & 
            invCOVxsi2000Data(s,i,j)*(xsi2000Sim(s,j)-xsi2000Data(s,j))
        end do
    end do
    
end do



do s = 0, NSECTORS

    D = CPPers2005Data(s,:,:)
    
    call dpotrf('L', NREG_Pers, D, NREG_Pers, info)
    if (info .ne. 0) then
        print*, 'CPPers2005Data is not symmetric positive definite.'
        print*, 'Aborting program.'
        pause
        stop
    end if


    do i = 1, NREG_Pers    
        do j = 1, NREG_Pers
            Loss_Function_Pers2005 = Loss_Function_Pers2005 + &
            (xsi2005Sim(s,i)-xsi2005Data(s,i))* & 
            invCOVxsi2005Data(s,i,j)*(xsi2005Sim(s,j)-xsi2005Data(s,j))
        end do
    end do
    
end do




do s = 0, NSECTORS

    E = CPFreqData(s,:,:)
    
    call dpotrf('L', NREG_Freq, E, NREG_Freq, info)
    if (info .ne. 0) then
        print*, 'CPFreqData is not symmetric positive definite.'
        print*, 'Aborting program.'
        pause
        stop
    end if


    do i = 1, NREG_Freq    
        do j = 1, NREG_Freq
            Loss_Function_Freq = Loss_Function_Freq + &
            (etaSim(s,i)-etaData(s,i))* & 
            invCOVetaData(s,i,j)*(etaSim(s,j)-etaData(s,j))
        end do
    end do
    
end do





do s = 1, NSECTORS
    Loss_Function_Sigma = Loss_Function_Sigma + & 
    (sigmaWageSim(s)**2-sigmaWageData(s)**2)**2 / VarSigma2Data(s)
end do

do s = 1, NSECTORS
    Loss_Function_SigmaDif = Loss_Function_SigmaDif + &
    (SigmaWageDifSim(s)**2-sigmaWageDifData(s)**2)**2 / VarSigma2WageDifData(s)
end do

do s = 1, NSECTORS
    do ed = 0, 1
        do t = FIRST_YR, LAST_YR
            Loss_Function_Labor = Loss_Function_Labor + 100.0*(LaborSharesSim(s,ed,t)-LaborSharesData(s,ed,t))**2
        end do
    end do
end do

do s = 1, NSECTORS
    do t = FIRST_YR, LAST_YR
        Loss_Function_Capital = Loss_Function_Capital + 100.0*(CapitalSharesSim(s,t)-CapitalSharesData(s,t))**2
    end do
end do

Loss_Function_Wage     = Loss_Function_Wage     / 10000.0
Loss_Function_Emp      = Loss_Function_Emp      / 10000.0
Loss_Function_Tr       = Loss_Function_Tr       / 10000.0
Loss_Function_Return   = Loss_Function_Return   / 10000.0
Loss_Function_Sigma    = Loss_Function_Sigma    / 10000.0
Loss_Function_SigmaDif = Loss_Function_SigmaDif / 10000.0
Loss_Function_Pers1998 = Loss_Function_Pers1998 / 10000.0
Loss_Function_Pers2000 = Loss_Function_Pers2000 / 10000.0
Loss_Function_Pers2005 = Loss_Function_Pers2005 / 10000.0
Loss_Function_Freq     = Loss_Function_Freq     / 10000.0


print*
print*, 'Loss Function Wage:          ' , Loss_Function_Wage
print*, 'Loss Function Emp:           ' , Loss_Function_Emp
print*, 'Loss Function Tr:            ' , Loss_Function_Tr
print*, 'Loss Function Return:        ' , Loss_Function_Return
print*, 'Loss Function Sigma:         ' , Loss_Function_Sigma
print*, 'Loss Function Sigma Dif:     ' , Loss_Function_SigmaDif
print*, 'Loss Function Pers 1998:     ' , Loss_Function_Pers1998
print*, 'Loss Function Pers 2000:     ' , Loss_Function_Pers2000
print*, 'Loss Function Pers 2005:     ' , Loss_Function_Pers2005
print*, 'Loss Function Freq:          ' , Loss_Function_Freq
print*, 'Loss Function Labor Shares:  ' , Loss_Function_Labor
print*, 'Loss Function K Shares:      ' , Loss_Function_Capital
print*


Loss_Function = Loss_Function_Wage + Loss_Function_Emp + Loss_Function_Tr + Loss_Function_Return + 20.0*Loss_Function_SigmaDif + &
                10.0*Loss_Function_Sigma + Loss_Function_Pers1998 + Loss_Function_Pers2000 + Loss_Function_Pers2005 + &
                Loss_Function_Freq + Loss_Function_Labor + Loss_Function_Capital


if (ISNAN(Loss_Function)) then
    Loss_Function = 999999.9999
end if


Call WriteCoefToFile(betaSim,SigmaWageDifSim,SigmaWageSim,gammaSim,phiSim,rhoSim,xsi1998Sim,xsi2000Sim,& 
                           xsi2005Sim,etaSim,Loss_Function)


end Function Loss_Function











subroutine WriteCoefToFile(betaSim,SigmaWageDifSim,SigmaWageSim,gammaSim,phiSim,rhoSim,xsi1998Sim,xsi2000Sim,& 
                           xsi2005Sim,etaSim,Loss_Function)

! ********************************************************************************
! This subroutine writes coefficients betaSim,sigmaWageSim,gammaSim,phiSim to file
! ********************************************************************************
     
USE Global_Data
                              
implicit none

real(KIND=DOUBLE), intent(in) :: betaSim(:,:), & 
                                 SigmaWageDifSim(:), &
                                 SigmaWageSim(:), &
                                 gammaSim(0:,:), & 
                                 phiSim(0:,0:,:), &
                                 rhoSim(0:,:), &
                                 xsi1998Sim(0:,:), &
                                 xsi2000Sim(0:,:), &
                                 xsi2005Sim(0:,:), &
                                 etaSim(0:,:), &
                                 Loss_Function
                                                      

integer s, s1, s2, i, gen

open(unit = 1, file = 'CompareCoefficients.csv')

do s = 1, NSECTORS
    do i = 1, NREG_Wage

        write(1,10) 'beta', s, i, betaData(s,i), betaSim(s,i), sqrt(COVbetaData(s,i,i))

    end do
end do


do s = 1, NSECTORS

    write(1,20) 'SigmaWageDif', s, (SigmaWageDifData(s))**2, (SigmaWageDifSim(s))**2, sqrt(VarSigma2WageDifData(s))

end do

do s = 1, NSECTORS

    write(1,20) 'SigmaWage', s, (SigmaWageData(s))**2, (SigmaWageSim(s))**2, sqrt(VarSigma2Data(s))

end do


do s = 0, NSECTORS
    do i = 1, NREG_Emp

        write(1,10) 'gamma', s, i, gammaData(s,i), gammaSim(s,i), sqrt(COVgammaData(s,i,i))

    end do
end do

do s1 = 0, NSECTORS
    do s2 = 0, NSECTORS
        do i = 1, NREG_Tr

            write(1,30) 'phi', s1, s2, i, phiData(s1,s2,i), phiSim(s1,s2,i), sqrt(COVphiData(s1,s2,i,i))

        end do
    end do
end do   

do s = 0, NSECTORS
    do i = 1, NREG_Return
       
        write(1,10) 'rho', s, i, rhoData(s,i), rhoSim(s,i), sqrt(COVrhoData(s,i,i))
        
    end do
end do


do s = 0, NSECTORS
    do i = 1, NREG_Pers
       
        write(1,10) 'xsi1998', s, i, xsi1998Data(s,i), xsi1998Sim(s,i), sqrt(COVxsi1998Data(s,i,i))
        
    end do
end do


do s = 0, NSECTORS
    do i = 1, NREG_Pers
       
        write(1,10) 'xsi2000', s, i, xsi2000Data(s,i), xsi2000Sim(s,i), sqrt(COVxsi2000Data(s,i,i))
        
    end do
end do


do s = 0, NSECTORS
    do i = 1, NREG_Pers
       
        write(1,10) 'xsi2005', s, i, xsi2005Data(s,i), xsi2005Sim(s,i), sqrt(COVxsi2005Data(s,i,i))
        
    end do
end do


do s = 0, NSECTORS
    do i = 1, NREG_Freq
       
        write(1,10) 'eta', s, i, etaData(s,i), etaSim(s,i), sqrt(COVetaData(s,i,i))
        
    end do
end do

close(1)



    write(101,100) PARAM_NUMBER, PARAM_VALUE, (betaSim(1,i), i = 1, NREG_Wage), Loss_Function
    write(102,100) PARAM_NUMBER, PARAM_VALUE, (betaSim(2,i), i = 1, NREG_Wage), Loss_Function
    write(103,100) PARAM_NUMBER, PARAM_VALUE, (betaSim(3,i), i = 1, NREG_Wage), Loss_Function
    write(104,100) PARAM_NUMBER, PARAM_VALUE, (betaSim(4,i), i = 1, NREG_Wage), Loss_Function
    write(105,100) PARAM_NUMBER, PARAM_VALUE, (betaSim(5,i), i = 1, NREG_Wage), Loss_Function
    write(106,100) PARAM_NUMBER, PARAM_VALUE, (betaSim(6,i), i = 1, NREG_Wage), Loss_Function
    write(107,100) PARAM_NUMBER, PARAM_VALUE, (betaSim(7,i), i = 1, NREG_Wage), Loss_Function
    
    write(200,100) PARAM_NUMBER, PARAM_VALUE, (gammaSim(0,i), i = 1, NREG_Emp), Loss_Function
    write(201,100) PARAM_NUMBER, PARAM_VALUE, (gammaSim(1,i), i = 1, NREG_Emp), Loss_Function
    write(202,100) PARAM_NUMBER, PARAM_VALUE, (gammaSim(2,i), i = 1, NREG_Emp), Loss_Function
    write(203,100) PARAM_NUMBER, PARAM_VALUE, (gammaSim(3,i), i = 1, NREG_Emp), Loss_Function
    write(204,100) PARAM_NUMBER, PARAM_VALUE, (gammaSim(4,i), i = 1, NREG_Emp), Loss_Function 
    write(205,100) PARAM_NUMBER, PARAM_VALUE, (gammaSim(5,i), i = 1, NREG_Emp), Loss_Function 
    write(206,100) PARAM_NUMBER, PARAM_VALUE, (gammaSim(6,i), i = 1, NREG_Emp), Loss_Function 
    write(207,100) PARAM_NUMBER, PARAM_VALUE, (gammaSim(7,i), i = 1, NREG_Emp), Loss_Function 
    
    write(300,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(0,0,i), i = 1, NREG_Tr), Loss_Function
    write(301,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(0,1,i), i = 1, NREG_Tr), Loss_Function
    write(302,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(0,2,i), i = 1, NREG_Tr), Loss_Function
    write(303,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(0,3,i), i = 1, NREG_Tr), Loss_Function
    write(304,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(0,4,i), i = 1, NREG_Tr), Loss_Function    
    write(305,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(0,5,i), i = 1, NREG_Tr), Loss_Function 
    write(306,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(0,6,i), i = 1, NREG_Tr), Loss_Function 
    write(307,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(0,7,i), i = 1, NREG_Tr), Loss_Function 
    
    write(310,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(1,0,i), i = 1, NREG_Tr), Loss_Function
    write(311,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(1,1,i), i = 1, NREG_Tr), Loss_Function
    write(312,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(1,2,i), i = 1, NREG_Tr), Loss_Function
    write(313,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(1,3,i), i = 1, NREG_Tr), Loss_Function
    write(314,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(1,4,i), i = 1, NREG_Tr), Loss_Function
    write(315,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(1,5,i), i = 1, NREG_Tr), Loss_Function
    write(316,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(1,6,i), i = 1, NREG_Tr), Loss_Function
    write(317,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(1,7,i), i = 1, NREG_Tr), Loss_Function
    
    write(320,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(2,0,i), i = 1, NREG_Tr), Loss_Function
    write(321,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(2,1,i), i = 1, NREG_Tr), Loss_Function
    write(322,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(2,2,i), i = 1, NREG_Tr), Loss_Function
    write(323,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(2,3,i), i = 1, NREG_Tr), Loss_Function
    write(324,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(2,4,i), i = 1, NREG_Tr), Loss_Function
    write(325,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(2,5,i), i = 1, NREG_Tr), Loss_Function
    write(326,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(2,6,i), i = 1, NREG_Tr), Loss_Function
    write(327,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(2,7,i), i = 1, NREG_Tr), Loss_Function
    
    write(330,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(3,0,i), i = 1, NREG_Tr), Loss_Function
    write(331,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(3,1,i), i = 1, NREG_Tr), Loss_Function
    write(332,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(3,2,i), i = 1, NREG_Tr), Loss_Function
    write(333,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(3,3,i), i = 1, NREG_Tr), Loss_Function
    write(334,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(3,4,i), i = 1, NREG_Tr), Loss_Function
    write(335,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(3,5,i), i = 1, NREG_Tr), Loss_Function
    write(336,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(3,6,i), i = 1, NREG_Tr), Loss_Function
    write(337,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(3,7,i), i = 1, NREG_Tr), Loss_Function
    
    write(340,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(4,0,i), i = 1, NREG_Tr), Loss_Function
    write(341,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(4,1,i), i = 1, NREG_Tr), Loss_Function
    write(342,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(4,2,i), i = 1, NREG_Tr), Loss_Function
    write(343,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(4,3,i), i = 1, NREG_Tr), Loss_Function
    write(344,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(4,4,i), i = 1, NREG_Tr), Loss_Function
    write(345,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(4,5,i), i = 1, NREG_Tr), Loss_Function
    write(346,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(4,6,i), i = 1, NREG_Tr), Loss_Function
    write(347,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(4,7,i), i = 1, NREG_Tr), Loss_Function
    
    write(350,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(5,0,i), i = 1, NREG_Tr), Loss_Function
    write(351,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(5,1,i), i = 1, NREG_Tr), Loss_Function
    write(352,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(5,2,i), i = 1, NREG_Tr), Loss_Function
    write(353,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(5,3,i), i = 1, NREG_Tr), Loss_Function
    write(354,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(5,4,i), i = 1, NREG_Tr), Loss_Function
    write(355,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(5,5,i), i = 1, NREG_Tr), Loss_Function
    write(356,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(5,6,i), i = 1, NREG_Tr), Loss_Function
    write(357,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(5,7,i), i = 1, NREG_Tr), Loss_Function
    
    write(360,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(6,0,i), i = 1, NREG_Tr), Loss_Function
    write(361,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(6,1,i), i = 1, NREG_Tr), Loss_Function
    write(362,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(6,2,i), i = 1, NREG_Tr), Loss_Function
    write(363,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(6,3,i), i = 1, NREG_Tr), Loss_Function
    write(364,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(6,4,i), i = 1, NREG_Tr), Loss_Function
    write(365,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(6,5,i), i = 1, NREG_Tr), Loss_Function
    write(366,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(6,6,i), i = 1, NREG_Tr), Loss_Function
    write(367,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(6,7,i), i = 1, NREG_Tr), Loss_Function
    
    write(370,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(7,0,i), i = 1, NREG_Tr), Loss_Function
    write(371,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(7,1,i), i = 1, NREG_Tr), Loss_Function
    write(372,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(7,2,i), i = 1, NREG_Tr), Loss_Function
    write(373,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(7,3,i), i = 1, NREG_Tr), Loss_Function
    write(374,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(7,4,i), i = 1, NREG_Tr), Loss_Function
    write(375,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(7,5,i), i = 1, NREG_Tr), Loss_Function
    write(376,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(7,6,i), i = 1, NREG_Tr), Loss_Function
    write(377,100) PARAM_NUMBER, PARAM_VALUE, (phiSim(7,7,i), i = 1, NREG_Tr), Loss_Function
    
    write(401,101) PARAM_NUMBER, PARAM_VALUE, (sigmaWageSim(s), s = 1, NSECTORS), Loss_Function
    
    write(501,101) PARAM_NUMBER, PARAM_VALUE, (SigmaWageDifSim(s), s = 1, NSECTORS), Loss_Function
    
    write(600,102) PARAM_NUMBER, PARAM_VALUE, (xsi1998Sim(0,i), i = 1, NREG_Pers), Loss_Function
    write(601,102) PARAM_NUMBER, PARAM_VALUE, (xsi1998Sim(1,i), i = 1, NREG_Pers), Loss_Function
    write(602,102) PARAM_NUMBER, PARAM_VALUE, (xsi1998Sim(2,i), i = 1, NREG_Pers), Loss_Function
    write(603,102) PARAM_NUMBER, PARAM_VALUE, (xsi1998Sim(3,i), i = 1, NREG_Pers), Loss_Function
    write(604,102) PARAM_NUMBER, PARAM_VALUE, (xsi1998Sim(4,i), i = 1, NREG_Pers), Loss_Function
    write(605,102) PARAM_NUMBER, PARAM_VALUE, (xsi1998Sim(5,i), i = 1, NREG_Pers), Loss_Function
    write(606,102) PARAM_NUMBER, PARAM_VALUE, (xsi1998Sim(6,i), i = 1, NREG_Pers), Loss_Function
    write(607,102) PARAM_NUMBER, PARAM_VALUE, (xsi1998Sim(7,i), i = 1, NREG_Pers), Loss_Function
    
    write(700,102) PARAM_NUMBER, PARAM_VALUE, (xsi2000Sim(0,i), i = 1, NREG_Pers), Loss_Function
    write(701,102) PARAM_NUMBER, PARAM_VALUE, (xsi2000Sim(1,i), i = 1, NREG_Pers), Loss_Function
    write(702,102) PARAM_NUMBER, PARAM_VALUE, (xsi2000Sim(2,i), i = 1, NREG_Pers), Loss_Function
    write(703,102) PARAM_NUMBER, PARAM_VALUE, (xsi2000Sim(3,i), i = 1, NREG_Pers), Loss_Function
    write(704,102) PARAM_NUMBER, PARAM_VALUE, (xsi2000Sim(4,i), i = 1, NREG_Pers), Loss_Function
    write(705,102) PARAM_NUMBER, PARAM_VALUE, (xsi2000Sim(5,i), i = 1, NREG_Pers), Loss_Function
    write(706,102) PARAM_NUMBER, PARAM_VALUE, (xsi2000Sim(6,i), i = 1, NREG_Pers), Loss_Function
    write(707,102) PARAM_NUMBER, PARAM_VALUE, (xsi2000Sim(7,i), i = 1, NREG_Pers), Loss_Function
    
    write(800,102) PARAM_NUMBER, PARAM_VALUE, (xsi2005Sim(0,i), i = 1, NREG_Pers), Loss_Function
    write(801,102) PARAM_NUMBER, PARAM_VALUE, (xsi2005Sim(1,i), i = 1, NREG_Pers), Loss_Function
    write(802,102) PARAM_NUMBER, PARAM_VALUE, (xsi2005Sim(2,i), i = 1, NREG_Pers), Loss_Function
    write(803,102) PARAM_NUMBER, PARAM_VALUE, (xsi2005Sim(3,i), i = 1, NREG_Pers), Loss_Function
    write(804,102) PARAM_NUMBER, PARAM_VALUE, (xsi2005Sim(4,i), i = 1, NREG_Pers), Loss_Function
    write(805,102) PARAM_NUMBER, PARAM_VALUE, (xsi2005Sim(5,i), i = 1, NREG_Pers), Loss_Function
    write(806,102) PARAM_NUMBER, PARAM_VALUE, (xsi2005Sim(6,i), i = 1, NREG_Pers), Loss_Function
    write(807,102) PARAM_NUMBER, PARAM_VALUE, (xsi2005Sim(7,i), i = 1, NREG_Pers), Loss_Function
    
    write(900,102) PARAM_NUMBER, PARAM_VALUE, (etaSim(0,i), i = 1, NREG_Freq), Loss_Function
    write(901,102) PARAM_NUMBER, PARAM_VALUE, (etaSim(1,i), i = 1, NREG_Freq), Loss_Function
    write(902,102) PARAM_NUMBER, PARAM_VALUE, (etaSim(2,i), i = 1, NREG_Freq), Loss_Function
    write(903,102) PARAM_NUMBER, PARAM_VALUE, (etaSim(3,i), i = 1, NREG_Freq), Loss_Function
    write(904,102) PARAM_NUMBER, PARAM_VALUE, (etaSim(4,i), i = 1, NREG_Freq), Loss_Function
    write(905,102) PARAM_NUMBER, PARAM_VALUE, (etaSim(5,i), i = 1, NREG_Freq), Loss_Function
    write(906,102) PARAM_NUMBER, PARAM_VALUE, (etaSim(6,i), i = 1, NREG_Freq), Loss_Function
    write(907,102) PARAM_NUMBER, PARAM_VALUE, (etaSim(7,i), i = 1, NREG_Freq), Loss_Function
    
    write(1000,100) PARAM_NUMBER, PARAM_VALUE, (rhoSim(0,i), i = 1, NREG_Return), Loss_Function
    write(1001,100) PARAM_NUMBER, PARAM_VALUE, (rhoSim(1,i), i = 1, NREG_Return), Loss_Function
    write(1002,100) PARAM_NUMBER, PARAM_VALUE, (rhoSim(2,i), i = 1, NREG_Return), Loss_Function
    write(1003,100) PARAM_NUMBER, PARAM_VALUE, (rhoSim(3,i), i = 1, NREG_Return), Loss_Function
    write(1004,100) PARAM_NUMBER, PARAM_VALUE, (rhoSim(4,i), i = 1, NREG_Return), Loss_Function
    write(1005,100) PARAM_NUMBER, PARAM_VALUE, (rhoSim(5,i), i = 1, NREG_Return), Loss_Function
    write(1006,100) PARAM_NUMBER, PARAM_VALUE, (rhoSim(6,i), i = 1, NREG_Return), Loss_Function
    write(1007,100) PARAM_NUMBER, PARAM_VALUE, (rhoSim(7,i), i = 1, NREG_Return), Loss_Function

10 format(a10, ',' , i4, ',', i4, ', ,',  3(f10.6,','))
20 format(a10, ',' , i4, ', , ,',  3(f10.6,','))
30 format(a10, ',' , i4, ',', i4, ',', i4, ',', 3(f10.6,','))

100 format(i4,',',26(f16.8,','))
101 format(i4,',',9(f16.8,','))
102 format(i4,',',23(f16.8,','))

end subroutine WriteCoefToFile




subroutine Write_Sim_Data(ChoiceSim,WageSim,ExperSim,typ,VSim,Cost1Sim,Cost2Sim,CohortWgt, &
                          CapitalSim,z,sk_sup_eq,rsk_eq)

! ****************************
! Write Generated Data to File
! ****************************

USE Global_Data

implicit none
          
integer          , intent(in) :: & 
ChoiceSim(:,FIRST_COH:,(FIRST_YR-9):), &
ExperSim(:,:,FIRST_COH:,FIRST_YR:), &
typ(:,FIRST_COH:)

real(KIND=DOUBLE), intent(in) :: WageSim(:,FIRST_COH:,FIRST_YR:), &
                                 VSim(:,FIRST_COH:,FIRST_YR:), &
                                 Cost1Sim(:,FIRST_COH:,FIRST_YR:), &
                                 Cost2Sim(:,FIRST_COH:,FIRST_YR:), &
                                 CohortWgt(FIRST_COH:,0:), &
                                 CapitalSim(:,FIRST_YR:), &
                                 z(:,FIRST_YR:), &
                                 sk_sup_eq(:,0:,FIRST_YR:), &
                                 rsk_eq(:,0:,FIRST_YR:)
                                 
integer coh, n, age, ed, s, t                                   
        


open(unit = 1, file = 'Data_set.csv')
do coh = FIRST_COH, LAST_COH
    do n = 1, NPOP
        do t = FIRST_YR, LAST_YR
            age = t - coh
            if (EducData(n,coh) <= 2) then
                ed = 0
            else if (EducData(n,coh) >= 3) then
                ed = 1
            end if
            if (age >= FIRST_AGE .and. age <= LAST_AGE) then
                write(1,1234) n, coh, t, EducData(n,coh), GenderData(n,coh), &
                 ChoiceSim(n,coh,coh+age), ChoiceSim(n,coh,coh+age-1), ChoiceSim(n,coh,coh+age-2), & 
                 ChoiceSim(n,coh,coh+age-3), ChoiceSim(n,coh,coh+age-4), ChoiceSim(n,coh,coh+age-5), & 
                 ChoiceSim(n,coh,coh+age-6), ChoiceSim(n,coh,coh+age-7), ChoiceSim(n,coh,coh+age-8), &
                 ChoiceSim(n,coh,coh+age-9), ExperSim(1,n,coh,coh+age), ExperSim(2,n,coh,coh+age), & 
                 ExperSim(3,n,coh,coh+age), ExperSim(4,n,coh,coh+age), ExperSim(5,n,coh,coh+age), &
                 ExperSim(6,n,coh,coh+age), ExperSim(7,n,coh,coh+age), WageSim(n,coh,coh+age), VSim(n,coh,coh+age), & 
                 Cost1Sim(n,coh,coh+age), Cost2Sim(n,coh,coh+age), CohortWgt(coh,ed)
            end if                              
        end do
    end do
end do    
close(1)


open(unit = 1, file = 'Init_Simulation.csv')
do coh = (LAST_YR-LAST_AGE), (LAST_YR-FIRST_AGE)
    age = LAST_YR - coh
    do n = 1, NPOP
        if (EducData(n,coh) <= 2) then
            ed = 0
        else if (EducData(n,coh) >= 3) then
            ed = 1
        end if
        write(1,4321) n, LAST_YR, coh, EducData(n,coh), GenderData(n,coh), & 
         ChoiceSim(n,coh,coh+age), ChoiceSim(n,coh,coh+age-1), ChoiceSim(n,coh,coh+age-2), &
         ChoiceSim(n,coh,coh+age-3), ChoiceSim(n,coh,coh+age-4), ChoiceSim(n,coh,coh+age-5), &
         ChoiceSim(n,coh,coh+age-6), ChoiceSim(n,coh,coh+age-7), ChoiceSim(n,coh,coh+age-8), &
         ChoiceSim(n,coh,coh+age-9), ed, typ(n,coh), CohortWgt(coh,ed)
    end do
end do
close(1)




open(unit = 1, file = 'Outcomes.csv')
write(1,1111) 'Year', 'Z1', 'Z2', 'Z3', 'Z4', 'Z5', 'Z6', 'Z7', & 
              'RSK_EQ1_0', 'RSK_EQ2_0', 'RSK_EQ3_0', 'RSK_EQ4_0', 'RSK_EQ5_0', 'RSK_EQ6_0', 'RSK_EQ7_0', & 
              'RSK_EQ1_1', 'RSK_EQ2_1', 'RSK_EQ3_1', 'RSK_EQ4_1', 'RSK_EQ5_1', 'RSK_EQ6_1', 'RSK_EQ7_1', & 
              'SK_SUP_EQ1_0', 'SK_SUP_EQ2_0', 'SK_SUP_EQ3_0', 'SK_SUP_EQ4_0', 'SK_SUP_EQ5_0', 'SK_SUP_EQ6_0', 'SK_SUP_EQ7_0', &
              'SK_SUP_EQ1_1', 'SK_SUP_EQ2_1', 'SK_SUP_EQ3_1', 'SK_SUP_EQ4_1', 'SK_SUP_EQ5_1', 'SK_SUP_EQ6_1', 'SK_SUP_EQ7_1', &
              'CAPSIM1', 'CAPSIM2', 'CAPSIM3', 'CAPSIM4', 'CAPSIM5', 'CAPSIM6', 'CAPSIM7', 'CAPDATA'
do t = FIRST_YR, LAST_YR
    write(1,1212) t, z(1,t), z(2,t), z(3,t), z(4,t), z(5,t), z(6,t), z(7,t), & 
                  rsk_eq(1,0,t), rsk_eq(2,0,t), rsk_eq(3,0,t), rsk_eq(4,0,t), rsk_eq(5,0,t), rsk_eq(6,0,t), rsk_eq(7,0,t), &
                  rsk_eq(1,1,t), rsk_eq(2,1,t), rsk_eq(3,1,t), rsk_eq(4,1,t), rsk_eq(5,1,t), rsk_eq(6,1,t), rsk_eq(7,1,t), &
                  sk_sup_eq(1,0,t), sk_sup_eq(2,0,t), sk_sup_eq(3,0,t), sk_sup_eq(4,0,t), sk_sup_eq(5,0,t), sk_sup_eq(6,0,t), sk_sup_eq(7,0,t), & 
                  sk_sup_eq(1,1,t), sk_sup_eq(2,1,t), sk_sup_eq(3,1,t), sk_sup_eq(4,1,t), sk_sup_eq(5,1,t), sk_sup_eq(6,1,t), sk_sup_eq(7,1,t), & 
                  CapitalSim(1,t), CapitalSim(2,t), CapitalSim(3,t), CapitalSim(4,t), CapitalSim(5,t), CapitalSim(6,t), CapitalSim(7,t), & 
                  CapitalData(t)
end do
close(1)


1111 format(44(a14,','))
1212 format(i5,',',43(f20.8,','))
1234 format(22(i5,','),5(f16.8,','))
4321 format(17(i5,','),f9.4,',')

end subroutine Write_Sim_Data


end module Loss_Function_MOD

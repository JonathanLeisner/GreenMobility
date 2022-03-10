program Main

! ************************************************************************
! This program estimates a dynamic labor market equilibrium with 7 sectors
! Main.f90 consists in:
! 1) Reading  coefficients and cross products matrices from the auxiliary
! models
! 2) Reading Initial Conditions, Value Added for each sector, etc...
! 3) Invoke two optimizers: NEWUOA and/or Nelder-Mead
! It contains two subroutines:
! Read_Coef
! Read_Data
! ************************************************************************


! ***********************
! Including other modules
! ***********************
USE MKL_VSL_TYPE
USE MKL_VSL
USE Global_Data
USE Loss_Function_MOD
USE Emax_MOD
USE Nelder_Mead
USE MPI


! **********************
! Variables Declarations
! **********************

implicit none


! Counting Variables
integer i, j, k, n, t, s, coh, ed, tp


! Model parameters
real(KIND=DOUBLE) theta(7), beta(NSECTORS,12), & 
                  sigma(0:NSECTORS), kappa(26), tau(0:NSECTORS), SigmaPref, &
                  omega(2:NTYPES,0:NSECTORS), lambda(NTYPES), gamma(2:NTYPES,0:NSECTORS+4), aa(NSECTORS,0:1,0:1), &
                  sigma_prod(NSECTORS), rsk_init(NSECTORS,0:1), param(NPARAM)
        

real(KIND=DOUBLE) f


! Variables for expansion factor and wage correction        
real(KIND=DOUBLE) sum_wage, correction


! Variables used in the portion that allows for the optimization
! In a subset of variables
integer NP          
real(KIND=DOUBLE), allocatable, dimension(:) :: param_in


! For NEWUOA
real(KIND=DOUBLE)  rho_end, rho_beg
integer            iprint, maxfun, NPT, SizeW
real(KIND=DOUBLE), allocatable, dimension(:) :: wrkspace


! For Nelder-Mead
real(KIND=DOUBLE) fvalue, ftol, simp
integer           MAXFCN, NLOOP, iquad, ifault
real(KIND=DOUBLE), allocatable, dimension(:) ::  stepsize, var


! Random Number Generation
TYPE (VSL_STREAM_STATE) :: stream
integer(kind=4)   errcode
integer           brng, method, method1, method2, method3, method4, seed
real(KIND=DOUBLE) rvec_eps(NSECTORS+1), rvec_eta(NSECTORS+1), scale, loc, unif(1)

! MPI Variables
integer ierr, size_mpi, rank, tag, source, count
integer status(MPI_STATUS_SIZE)

! *******************
! Program Starts Here
! *******************

Call MPI_INIT(ierr)
Call MPI_COMM_SIZE(MPI_COMM_WORLD,size_mpi,ierr)
Call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)


if (rank == 0) then

! *******************
! Program Starts Here
! *******************

    ! *****************************************
    ! Read Coefficients of the Auxiliary Models
    ! *****************************************
    
    Call Read_Coef



    ! ***************
    ! Read other Data
    ! ***************

    Call Read_Data




    !***************************************************************************************
    ! Calculating the expansion factor Using wage data for 1995 only
    ! Calibration: WageBill(1995) = 0.65*Y(1995)
    ! Expansion_Factor = 0.65*Y(1995)/Sum_wages
    ! Where Sum_wages = cohort weighted sum of wages in the initial conditions
    ! Expansion Factor expands the Wage Bill in the Data in order to get an estimate of the 
    ! Wage Bill for the population
    !***************************************************************************************

    sum_wage   = 0.0

    do t = FIRST_YR, LAST_YR
        do coh = (t-LAST_AGE), (t-FIRST_AGE)
            do n = 1, NPOP
                if (WageData_ExpFactor(n,coh,t) > 0.0) then
                    if (EducData_ExpFactor(n,coh,t) <= 2) then
                        sum_wage = sum_wage + & 
                        (real(CohortSizeData(coh,0)) / & 
                        real(CohortSizeData(FIRST_COH,0)))*WageData_ExpFactor(n,coh,t)
                    else if (EducData_ExpFactor(n,coh,t) >= 3) then
                        sum_wage = sum_wage + & 
                        (real(CohortSizeData(coh,1)) / & 
                        real(CohortSizeData(FIRST_COH,0)))*WageData_ExpFactor(n,coh,t)
                    end if
                end if
            end do
        end do
    end do

    ExpansionFactorGlobal = 0.65*sum(OutputData) / sum_wage


    do t = FIRST_YR, LAST_YR
        rKData(t) = rKData(t) / 100.0
    end do



    ! ***********************************************************
    ! Correction factor in order to have sum(wages) = 0.65*output
    ! This will be used when I compute sector specific labor
    ! shares. On the aggregate and over the sample period, the 
    ! labor share is 0.65 However, the sector specific labor 
    ! shares can change over time
    ! ***********************************************************

    correction = 0.65*sum(OutputData(:,:)) / sum(TotalHrWageData(:,:,:))



    ! *************
    ! Labor Shares
    ! *************

    do t = FIRST_YR, LAST_YR
        do s = 1, NSECTORS
            LaborSharesData(s,0,t) = correction*sum(TotalHrWageData(s,t,1:2)) / OutputData(s,t)
            LaborSharesData(s,1,t) = correction*sum(TotalHrWageData(s,t,3:4)) / OutputData(s,t)
        end do
    end do

    ! Read Human Capital (Labor) Shares
    open(unit = 1, file = 'Labor_Shares.csv')
    write(1,9876)
    do ed = 0, 1
        do t = FIRST_YR, LAST_YR
            write(1,8765) t, ed, (LaborSharesData(s,ed,t), s = 1, NSECTORS)
        end do
    end do
    close(1)

    9876 format('Year, Educ Level, Agr/Min, LT Manuf, HT Manuf, Const, Trade, Trans/Util, Service')
    8765 format(2(i5,','),7(f7.4,','))
         
    do t = FIRST_YR, LAST_YR     
        do s = 1, NSECTORS
            CapitalSharesData(s,t) = 1 - LaborSharesData(s,0,t) - LaborSharesData(s,1,t)
        end do
    end do
    

    !****************************
    ! Initializing the parameters
    !****************************

    open(unit = 1, file ='Starting_Point.csv')

        do i = 1, 7
            read(1,*) n, theta(i), scaling_global(n)
        end do

        do s = 1, NSECTORS
            do i = 1, 12
                read(1,*) n, beta(s,i), scaling_global(n)
            end do
        end do    

        do s = 0, NSECTORS
            read(1,*) n, sigma(s), scaling_global(n)
        end do

        do i = 1, 26
            read(1,*) n, kappa(i), scaling_global(n)
        end do  
    
        do i = 1, 6
            read(1,*) n, tau(i+1), scaling_global(n)
        end do
    
        read(1,*) n, SigmaPref, scaling_global(n)
        
        do tp = 2, NTYPES
            do i = 0, NSECTORS
                read(1,*) n, omega(tp,i), scaling_global(n)
            end do
        end do
        
        do i = 2, NTYPES
            read(1,*) n, lambda(i), scaling_global(n)
        end do
    
        do i = 2, NTYPES
            do s = 0, NSECTORS+4
                read(1,*) n, gamma(i,s), scaling_global(n)
            end do
        end do
        
        do ed = 0, 1
            do s = 1, NSECTORS
                do i = 0, 1
                    read(1,*) n, aa(s,ed,i), scaling_global(n)
                end do
            end do
        end do
        
        do s = 1, NSECTORS
            read(1,*), n, sigma_prod(s), scaling_global(n)
        end do

        do ed = 0, 1
            do s = 1, NSECTORS
                read(1,*), n, rsk_init(s,ed), scaling_global(n)  
            end do
        end do

    close(1)


    ! *******************************************************************************
    ! Scaling the Starting Point so that each component is approximately equal to 1.0
    ! *******************************************************************************=

    param(1:7)      = theta(1:7)       / scaling_global(1:7)
    param(8:19)     = beta(1,1:12)     / scaling_global(8:19)
    param(20:31)    = beta(2,1:12)     / scaling_global(20:31)
    param(32:43)    = beta(3,1:12)     / scaling_global(32:43)
    param(44:55)    = beta(4,1:12)     / scaling_global(44:55) 
    param(56:67)    = beta(5,1:12)     / scaling_global(56:67) 
    param(68:79)    = beta(6,1:12)     / scaling_global(68:79) 
    param(80:91)    = beta(7,1:12)     / scaling_global(80:91) 
    param(92:99)    = log(sigma(0:7))  / scaling_global(92:99)
    param(100:125)  = kappa(1:26)      / scaling_global(100:125)
    param(126:131)  = tau(2:7)         / scaling_global(126:131)
    param(132)      = log(SigmaPref)   / scaling_global(132)
    param(133:140)  = omega(2,0:7)     / scaling_global(133:140)
    param(141:148)  = omega(3,0:7)     / scaling_global(141:148)
    param(149:150)  = lambda(2:3)      / scaling_global(149:150)
    param(151:162)  = gamma(2,0:11)    / scaling_global(151:162)
    param(163:174)  = gamma(3,0:11)    / scaling_global(163:174)
    param(175)      = log(aa(1,0,0))   / scaling_global(175)
    param(176)      = aa(1,0,1)        / scaling_global(176)
    param(177)      = log(aa(2,0,0))   / scaling_global(177)
    param(178)      = aa(2,0,1)        / scaling_global(178)
    param(179)      = log(aa(3,0,0))   / scaling_global(179)
    param(180)      = aa(3,0,1)        / scaling_global(180)
    param(181)      = log(aa(4,0,0))   / scaling_global(181)
    param(182)      = aa(4,0,1)        / scaling_global(182)
    param(183)      = log(aa(5,0,0))   / scaling_global(183)
    param(184)      = aa(5,0,1)        / scaling_global(184)
    param(185)      = log(aa(6,0,0))   / scaling_global(185)
    param(186)      = aa(6,0,1)        / scaling_global(186)
    param(187)      = log(aa(7,0,0))   / scaling_global(187)
    param(188)      = aa(7,0,1)        / scaling_global(188)    
    param(189)      = log(aa(1,1,0))   / scaling_global(189)
    param(190)      = aa(1,1,1)        / scaling_global(190)
    param(191)      = log(aa(2,1,0))   / scaling_global(191)
    param(192)      = aa(2,1,1)        / scaling_global(192)
    param(193)      = log(aa(3,1,0))   / scaling_global(193)
    param(194)      = aa(3,1,1)        / scaling_global(194)
    param(195)      = log(aa(4,1,0))   / scaling_global(195)
    param(196)      = aa(4,1,1)        / scaling_global(196)
    param(197)      = log(aa(5,1,0))   / scaling_global(197)
    param(198)      = aa(5,1,1)        / scaling_global(198)
    param(199)      = log(aa(6,1,0))   / scaling_global(199)
    param(200)      = aa(6,1,1)        / scaling_global(200)
    param(201)      = log(aa(7,1,0))   / scaling_global(201)
    param(202)      = aa(7,1,1)        / scaling_global(202)
    param(203:209)  = log(sigma_prod(1:7)) / scaling_global(203:209)
    param(210:216)  = log(rsk_init(1:7,0)) / scaling_global(210:216)
    param(217:223)  = log(rsk_init(1:7,1)) / scaling_global(217:223)

           
    !*****************************************************************
    ! Generate the random draws for the simulation performed in
    ! Subroutine SMM_OBJ_FCN once and for all
    !*****************************************************************

    brng   = VSL_BRNG_MT19937
    method = VSL_METHOD_DGAUSSIAN_BOXMULLER
    seed   = 112233445

    ! ***** Initializing *****
    errcode = vslnewstream(stream, brng,seed)

    ! ***** Warming Up *****
    do j = 1, WARMUP
        errcode = vdrnggaussian(method, stream, NSECTORS+1, &
                                rvec_eps, real(0.0,DOUBLE), real(1.0,DOUBLE))
    end do 

    do coh = FIRST_COH, LAST_COH
        do n = 1, NPOP
            do t = FIRST_YR, LAST_YR
                errcode = vdrnggaussian(method, stream, NSECTORS+1, &
                                        rvec_eps, real(0.0,DOUBLE), real(1.0,DOUBLE))
                eps_global(n,coh,t,0:NSECTORS) = rvec_eps
            end do
        end do
    end do

    ! ***** Deinitialize *****
    errcode = vsldeletestream(stream)


    !*****************************************************************
    ! Normalized Preference Shocks from Gunbel distribution
    !*****************************************************************

    brng    = VSL_BRNG_MT19937
    method3 = VSL_METHOD_DGUMBEL_ICDF
    seed    = 223344556

    ! ***** Initializing *****
    errcode = vslnewstream(stream,brng,seed)

    loc   = 0.0
    scale = 1.0

    ! ***** Initializing *****
    errcode = vslnewstream(stream, brng, seed)

    ! ***** Warming Up *****
    do j = 1, WARMUP
        errcode = vdrnggumbel(method3, stream, NSECTORS+1, rvec_eta, loc, scale)
    end do 

    do coh = FIRST_COH, LAST_COH
        do n = 1, NPOP
            do t = FIRST_YR, LAST_YR
                errcode = vdrnggumbel(method3, stream, NSECTORS+1, rvec_eta, loc, scale)
                eta_global(n,coh,t,0:NSECTORS) = rvec_eta
            end do
        end do
    end do

    ! ***** Deinitialize *****
    errcode = vsldeletestream(stream)
    
    
    
    !*****************************************************************
    ! Uniform Draws for the probability of types
    !*****************************************************************

    ! ***** Initializing *****
    brng    = VSL_BRNG_MT19937
    seed    = 98765456
    method4 = VSL_METHOD_IUNIFORM_STD
    errcode = vslnewstream(stream, brng, seed)

    ! ***** Warming Up *****
    do j = 1, WARMUP
        errcode = vdrnguniform(method4, stream, 1, unif, &
                  real(0.0,DOUBLE), real(1.0,DOUBLE))
    end do 

    do coh = FIRST_COH, LAST_COH
        do n = 1, NPOP
            errcode = vdrnguniform(method4, stream, 1, unif, &
                real(0.0,DOUBLE), real(1.0,DOUBLE))    
            unif_global(n,coh,:) = unif
        end do
    end do      

    ! ***** Deinitialize *****
    errcode = vsldeletestream(stream)

end if


Call MPI_BCAST(param,NPARAM,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

count = size(EducData)
Call MPI_BCAST(EducData,count,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

count = size(GenderData)
Call MPI_BCAST(GenderData,count,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

Call MPI_BCAST(ExpansionFactorGlobal,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

count = size(CohortSizeData)
Call MPI_BCAST(CohortSizeData,count,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! Default: code does not compute Standard Errors
! This variable is reset to 1 if subroutine NumDer is called
StdErrors_FLAG = 0
    
! **********************************************************
! Set a mask if one wants the optimization to be carried out
! on a subset of parameters
! **********************************************************

param_global = param
mask_global = 0
!mask_global(13:19) = 1
!mask_global(25:31) = 1
!mask_global(37:43) = 1
!mask_global(49:55) = 1
!mask_global(61:67) = 1
!mask_global(73:79) = 1
!mask_global(85:91) = 1
!mask_global(100:132) = 1
!mask_global(1:174) = 1
mask_global(1:209) = 1

! Number of parameters in the subset of parameters
NP = sum(mask_global)

allocate(param_in(NP))
k = 1
do i = 1, NPARAM
    if(mask_global(i) == 1) then
        param_in(k) = param(i)
        k = k + 1
    end if
end do

allocate(stepsize(NP))
allocate(var(NP))

 ! *************************
 ! Optimization: Nelder Mead
 ! *************************

ALGORITHM     = 2
BATCH         = 1
FUNCTION_ITER = 0
                 
ftol     = 0.001
stepsize = 0.25
MAXFCN   = 3300
iprint   = -1
NLOOP    = 2*NP
iquad    = 1
simp     = 0.000001                 
              
Call SMM_OBJ_FCN(param_in,F)
   
Call minim(param_in, stepsize, NP, fvalue, MAXFCN, iprint, ftol, NLOOP, iquad,  &
                 simp, var, SMM_OBJ_FCN, ifault) 



! ********************
! Optimization: NEWUOA
! ********************
!
!ALGORITHM     = 1
!BATCH         = 2
!FUNCTION_ITER = 0
!
!iprint  = 0
!rho_end = 0.001
!rho_beg = 0.25
!maxfun  = 3000
!NPT     = 2*NP+1
!SizeW   = (NPT+13)*(NPT+NP)+3*NP*(NP+3)
!
!allocate(wrkspace(SizeW))
!
!wrkspace       = 0.0
!
!Call SMM_OBJ_FCN(param_in,F)
!
!call newuoa(NP,npt,param_in,rho_beg,rho_end,iprint,maxfun,wrkspace)
              
Call SMM_OBJ_FCN(param_in,F)



! ***************************
! Scaling Back the parameters
! ***************************

if (rank == 0) then

k = 1
do i = 1, NPARAM
    if(mask_global(i) == 1) then
        param(i) = param_in(k)
        k = k + 1
    end if
end do

param(1:91)     = param(1:91)         * scaling_global(1:91)
param(92:99)    = exp(param(92:99)    * scaling_global(92:99))
param(100:131)  = param(100:131)      * scaling_global(100:131)
param(132)      = exp(param(132)      * scaling_global(132))
param(133:174)  = param(133:174)      * scaling_global(133:174)
param(175)      = exp(param(175)      * scaling_global(175))
param(176)      = param(176)          * scaling_global(176)
param(177)      = exp(param(177)      * scaling_global(177))
param(178)      = param(178)          * scaling_global(178)
param(179)      = exp(param(179)      * scaling_global(179))
param(180)      = param(180)          * scaling_global(180)
param(181)      = exp(param(181)      * scaling_global(181))
param(182)      = param(182)          * scaling_global(182)
param(183)      = exp(param(183)      * scaling_global(183))
param(184)      = param(184)          * scaling_global(184)
param(185)      = exp(param(185)      * scaling_global(185))
param(186)      = param(186)          * scaling_global(186)
param(187)      = exp(param(187)      * scaling_global(187))
param(188)      = param(188)          * scaling_global(188)
param(189)      = exp(param(189)      * scaling_global(189))
param(190)      = param(190)          * scaling_global(190)
param(191)      = exp(param(191)      * scaling_global(191))
param(192)      = param(192)          * scaling_global(192)
param(193)      = exp(param(193)      * scaling_global(193))
param(194)      = param(194)          * scaling_global(194)
param(195)      = exp(param(195)      * scaling_global(195))
param(196)      = param(196)          * scaling_global(196)
param(197)      = exp(param(197)      * scaling_global(197))
param(198)      = param(198)          * scaling_global(198)
param(199)      = exp(param(199)      * scaling_global(199))
param(200)      = param(200)          * scaling_global(200)
param(201)      = exp(param(201)      * scaling_global(201))
param(202)      = param(202)          * scaling_global(202)
param(203:223)  = exp(param(203:223)  * scaling_global(203:223))


! Save the optimal parameter vector
open(unit = 1, file = 'Opt_Param.csv')
    do i = 1, NPARAM
        write(1,*) i, ',', param(i)
    end do
close(1)

end if

deallocate(param_in,stepsize,var)

Call MPI_FINALIZE(ierr)



!************
! SUBROUTINES
!************

Contains     

Subroutine Read_Data

! ***********************************************************
! Read individual wage data, cohort size, initial conditions,
! value added, returns to capital, total wages
! ***********************************************************

implicit none
    
! Counting variables
integer i, j, n, t, k, l, ls, s, s1, s2, e, a, g, it, &
        coh, year, x, x1, x2, sector
    
! IOSTAT variable for reading files
integer status



! ******************
! Initial Conditions
! ******************

WageData = -999

open(unit = 1 , file = 'initial_conditions.csv')
read(1,*)
do 
    read(1,*,IOSTAT=status) n, coh, FirstYearData(n,coh), AgeData(n,coh), &
                            EducData(n,coh), GenderData(n,coh), &
                            ChoiceData(n,coh,coh+AgeData(n,coh)), &
                            ChoiceData(n,coh,coh+AgeData(n,coh)-1), & 
                            ChoiceData(n,coh,coh+AgeData(n,coh)-2), &
                            ChoiceData(n,coh,coh+AgeData(n,coh)-3), &
                            ChoiceData(n,coh,coh+AgeData(n,coh)-4), &
                            ChoiceData(n,coh,coh+AgeData(n,coh)-5), &
                            ChoiceData(n,coh,coh+AgeData(n,coh)-6), &
                            ChoiceData(n,coh,coh+AgeData(n,coh)-7), &
                            ChoiceData(n,coh,coh+AgeData(n,coh)-8), &
                            ChoiceData(n,coh,coh+AgeData(n,coh)-9), &
                            WageData(n,coh,coh+AgeData(n,coh))                           
    if (status /= 0 ) exit
end do
close(1)


WageData_ExpFactor = -999

open(unit = 2 , file = 'wagedata_expfactor.csv')
read(2,*)
do 
    read(2,*,IOSTAT=status) n, coh, year, s, EducData_ExpFactor(n,coh,year), WageData_ExpFactor(n,coh,year)
    if (status /= 0 ) exit
end do
close(2)





!*************
! Cohort Sizes
!*************

open(unit = 3 , file = 'cohort_sizes.csv')
read(3,*)
do 
    read(3,*,IOSTAT=status) year, ed, CohortSizeData(year,ed)
    if (status /= 0 ) exit
end do
close(3)



! ******************
! Value Added Series
! ******************

open(unit = 4 , file = 'va_series.csv')
read(4,*)
do 
    read(4,*,IOSTAT=status) year, OutputData(1,year), OutputData(2,year), & 
                            OutputData(3,year), OutputData(4,year), OutputData(5,year), &
                            OutputData(6,year), OutputData(7,year)
    if (status /= 0 ) exit
end do
close(4)



! *****************************
! Capital and Return to Capital
! *****************************

open(unit = 5 , file = 'capital and return to capital.csv')
!read(5,*)
do
    read(5,*,IOSTAT=status) year, CapitalData(year), rKData(year)
    if (status /= 0 ) exit
end do
close(5)



! **********
! Total Wage
! **********

open(unit = 6 , file = 'total_hr_wage.csv')
    read(6,*)
    do 
        read(6,*,IOSTAT=status) year, sector, ed, TotalHrWageData(sector,year,ed), &
                                TotalWageData(sector,year,ed)
        if (status /= 0 ) exit
    end do
close(6)


end subroutine Read_Data



Subroutine Read_Coef

! ***********************************************************************
! Read the coefficients obtained by estimating the auxiliary models
! to the data. Also reads cross product matrices and covariance matrices.
! The covariance matrices are computed assuming homokedasticity.
! ***********************************************************************

implicit none

integer i, s, s1, s2, k, info
real(KIND=DOUBLE) blah, A(NREG_Wage,NREG_Wage), B(NREG_Emp,NREG_Emp), &
                  C(NREG_Tr,NREG_Tr), D(NREG_Pers,NREG_Pers), E(NREG_Freq,NREG_Freq), &
                  UW(NREG_Wage,NREG_Wage), UE(NREG_Emp,NREG_Emp), UT(NREG_Tr,NREG_Tr), &
                  UPers(NREG_Pers,NREG_Pers), UFreq(NREG_Freq,NREG_Freq)
character char, char2
character*25 filename, blah2



! ********************
! Log-Wage Regressions
! ********************

do s = 1, NSECTORS
    write(char,1000) s
    filename = 'wagesector'
    filename = filename(1:10) // char // '.csv'
    open(unit = 1, file = filename)
        read(1,*)
        read(1,*) sigmaWageData(s), (betaData(s,i), i = 1, NREG_Wage)
        do i = 1, NREG_Wage
            read(1,*) blah, (COVbetaData(s,i,j), j = 1, NREG_Wage)
        end do
    close(1)
    
    filename = 'cp' // filename
    open(unit = 2, file = filename)
        read(2,*)
        do i = 1, NREG_Wage
            read(2,*) NobsWageData(s), blah2, (CPWageData(s,i,j), j = 1, NREG_Wage)
        end do
    close(2)
    
    ! Check if CPWageData is positive definite
    
    A = CPWageData(s,:,:)
    
    call dpotrf('L', NREG_Wage, A, NREG_Wage, info)
    if (info .ne. 0) then
        print*, 'CPWageData is not symmetric positive definite.'
        print*, 'Aborting program.'
        pause
        stop
    end if
    
    UW = COVbetaData(s,:,:)
    
    call dpotrf('U',NREG_Wage,UW,NREG_Wage,info)    
    
    invCOVbetaData(s,:,:) = UW
    
    call dpotri('U',NREG_Wage,invCOVbetaData(s,:,:),NREG_Wage,info)
    
    do i = 1, NREG_Wage
        do j = 1, NREG_Wage
            if (i > j) then
                invCOVbetaData(s,i,j) = invCOVbetaData(s,j,i)
            end if
        end do
    end do
    
end do


open(unit = 1, file = 'var_residuals.csv')

    read(1,*)
    do s = 1, NSECTORS
        read(1,*) i, VarSigma2Data(i), n
        VarSigma2Data(i) = VarSigma2Data(i) / n
    end do       

close(1)

! ********************************************************************
! ********************************************************************

! **************************
! Log-Wage First Differences
! **************************

do s = 1, NSECTORS
    write(char,1000) s
    filename = 'wagedif'
    filename = filename(1:7) // char // '.csv'
    open(unit = 1, file = filename)
        read(1,*)
        read(1,*) sigmaWageDifData(s)
    close(1)
end do    

open(unit = 1, file = 'var_residuals_wagedif.csv')

    read(1,*)
    do s = 1, NSECTORS
        read(1,*) i, VarSigma2WageDifData(i), n
        VarSigma2WageDifData(i) = VarSigma2WageDifData(i) / n
    end do       

close(1)

! ********************************************************************
! ********************************************************************

! ***************************
! Sectoral Choice Regressions
! ***************************

do s = 0, NSECTORS
    write(char,1000) s
    filename = 'emp'
    filename = filename(1:3) // char // '.csv'
    open(unit = 1, file = filename)
        read(1,*)
        read(1,*) sigmaEmpData(s), (gammaData(s,i), i = 1, NREG_Emp)
        do i = 1, NREG_Emp
            read(1,*) blah, (COVgammaData(s,i,j), j = 1, NREG_Emp)
        end do
    close(1)
    
    filename = 'cp' // filename
    open(unit = 2, file = filename)
        read(2,*)
        do i = 1, NREG_Emp
            read(2,*) NobsEmpData(s), blah2, (CPEmpData(s,i,j), j = 1, NREG_Emp)
        end do
    close(2)    
    
    ! Check if CPEmpData is positive definite
    
    B = CPEmpData(s,:,:)
    
    call dpotrf('L', NREG_Emp, B, NREG_Emp, info)
    if (info .ne. 0) then
        print*, 'CPEmpData is not symmetric positive definite.'
        pause
        print*, 'Aborting program.'
        stop
    end if
    
    UE = COVgammaData(s,:,:)
    
    call dpotrf('U',NREG_Emp,UE,NREG_Emp,info)    
    
    invCOVgammaData(s,:,:) = UE
    
    call dpotri('U',NREG_Emp,invCOVgammaData(s,:,:),NREG_Emp,info)
    
    do i = 1, NREG_Emp
        do j = 1, NREG_Emp
            if (i > j) then
                invCOVgammaData(s,i,j) = invCOVgammaData(s,j,i)
            end if
        end do
    end do
    
end do

! ***************************************************************
! ***************************************************************


! *******************
! Return Regressions
! *******************

do s = 0, NSECTORS
    write(char,1000) s
    filename = 'Return'
    filename = filename(1:6) // char // '.csv'
    open(unit = 1, file = filename)
        read(1,*)
        read(1,*) sigmaReturnData(s), (rhoData(s,i), i = 1, NREG_Return)
        do i = 1, NREG_Return
            read(1,*) blah, (COVrhoData(s,i,j), j = 1, NREG_Return)
        end do
    close(1)
    
    filename = 'cp' // filename
    open(unit = 2, file = filename)
        read(2,*)
        do i = 1, NREG_Return
            read(2,*) NobsReturnData(s), blah2, (CPReturnData(s,i,j), j = 1, NREG_Emp)
        end do
    close(2)    
    
    ! Check if CPreturnData is positive definite
    
    B = CPReturnData(s,:,:)
    
    call dpotrf('L', NREG_Return, B, NREG_Return, info)
    if (info .ne. 0) then
        print*, 'CPReturnData is not symmetric positive definite.'
        pause
        print*, 'Aborting program.'
        stop
    end if
    
    UE = COVrhoData(s,:,:)
    
    call dpotrf('U',NREG_Return,UE,NREG_Return,info)    
    
    invCOVrhoData(s,:,:) = UE
    
    call dpotri('U',NREG_Return,invCOVrhoData(s,:,:),NREG_Return,info)
    
    do i = 1, NREG_Return
        do j = 1, NREG_Return
            if (i > j) then
                invCOVrhoData(s,i,j) = invCOVrhoData(s,j,i)
            end if
        end do
    end do
    
end do

! ***************************************************************
! ***************************************************************


! **********************
! Transition Regressions
! **********************

do s1 = 0, NSECTORS
    do s2 = 0, NSECTORS
        write(char,1000) s1
        write(char2,1000) s2
        filename = 'tr'
        filename = filename(1:2) // char // char2 // '.csv'
        open(unit = 1, file = filename)
            read(1,*)
            read(1,*) sigmaTrData(s1,s2), (phiData(s1,s2,i), i = 1, NREG_Tr)
            do i = 1, NREG_Tr
                read(1,*) blah, (COVphiData(s1,s2,i,j), j = 1, NREG_Tr)
            end do
        close(1)
        
        filename = 'cp' // filename
        open(unit = 2, file = filename)
            read(2,*)
            do i = 1, NREG_Tr
                read(2,*) NobsTrData(s1,s2), blah2, (CPTrData(s1,s2,i,j), j = 1, NREG_Tr)
            end do
        close(2)
        
        ! Check if CPTrData is positive definite
        
        C = CPTrData(s1,s2,:,:)
        
        call dpotrf('L', NREG_Tr, C, NREG_Tr, info)
        if (info .ne. 0) then
            print*, 'CPTrData is not symmetric positive definite.'
            print*, 'Aborting program.'
            pause
            stop
        end if
        
        UT = COVphiData(s1,s2,:,:)
    
        call dpotrf('U',NREG_Tr,UT,NREG_Tr,info)    
    
        invCOVphiData(s1,s2,:,:) = UT
    
        call dpotri('U',NREG_Emp,invCOVphiData(s1,s2,:,:),NREG_Tr,info)
    
        do i = 1, NREG_Tr
            do j = 1, NREG_Tr
                if (i > j) then
                    invCOVphiData(s1,s2,i,j) = invCOVphiData(s1,s2,j,i)
                end if
            end do
        end do
        
    end do
end do

! ***************************************************************
! ***************************************************************


! ***********************
! Persistence Regressions
! ***********************

do s = 0, NSECTORS
    write(char,1000) s
    filename = 'Persistence'
    filename = filename(1:11) // char // '_1998.csv'
    open(unit = 1, file = filename)
        read(1,*)
        read(1,*) sigmaPers1998Data(s), (xsi1998Data(s,i), i = 1, NREG_Pers)
        do i = 1, NREG_Pers
            read(1,*) blah, (COVxsi1998Data(s,i,j), j = 1, NREG_Pers)
        end do
    close(1)
    
    filename = 'cp' // filename
    open(unit = 2, file = filename)
        read(2,*)
        do i = 1, NREG_Pers
            read(2,*) NobsPers1998Data(s), blah2, (CPPers1998Data(s,i,j), j = 1, NREG_Pers)
        end do
    close(2)
    
    
    ! Check if CPPers1998Data is positive definite
        
    D = CPPers1998Data(s,:,:)
        
    call dpotrf('L', NREG_Pers, D, NREG_Pers, info)
    if (info .ne. 0) then
        print*, 'CPPers1998Data is not symmetric positive definite.'
        print*, 'Aborting program.'
        pause
        stop
    end if
        
    UPers = COVxsi1998Data(s,:,:)
    
    call dpotrf('U',NREG_Pers,UPers,NREG_Pers,info)    
    
    invCOVxsi1998Data(s,:,:) = UPers
    
    call dpotri('U',NREG_Pers,invCOVxsi1998Data(s,:,:),NREG_Pers,info)
    
    do i = 1, NREG_Pers
        do j = 1, NREG_Pers
            if (i > j) then
                invCOVxsi1998Data(s,i,j) = invCOVxsi1998Data(s,j,i)
            end if
        end do
    end do
    
end do    



do s = 0, NSECTORS
    write(char,1000) s
    filename = 'Persistence'
    filename = filename(1:11) // char // '_2000.csv'
    open(unit = 1, file = filename)
        read(1,*)
        read(1,*) sigmaPers2000Data(s), (xsi2000Data(s,i), i = 1, NREG_Pers)
        do i = 1, NREG_Pers
            read(1,*) blah, (COVxsi2000Data(s,i,j), j = 1, NREG_Pers)
        end do
    close(1)
    
    filename = 'cp' // filename
    open(unit = 2, file = filename)
        read(2,*)
        do i = 1, NREG_Pers
            read(2,*) NobsPers2000Data(s), blah2, (CPPers2000Data(s,i,j), j = 1, NREG_Pers)
        end do
    close(2)
    
    
    ! Check if CPPers2000Data is positive definite
        
    D = CPPers2000Data(s,:,:)
        
    call dpotrf('L', NREG_Pers, D, NREG_Pers, info)
    if (info .ne. 0) then
        print*, 'CPPers2000Data is not symmetric positive definite.'
        print*, 'Aborting program.'
        pause
        stop
    end if
        
    UPers = COVxsi2000Data(s,:,:)
    
    call dpotrf('U',NREG_Pers,UPers,NREG_Pers,info)    
    
    invCOVxsi2000Data(s,:,:) = UPers
    
    call dpotri('U',NREG_Pers,invCOVxsi2000Data(s,:,:),NREG_Pers,info)
    
    do i = 1, NREG_Pers
        do j = 1, NREG_Pers
            if (i > j) then
                invCOVxsi2000Data(s,i,j) = invCOVxsi2000Data(s,j,i)
            end if
        end do
    end do
    
end do 



do s = 0, NSECTORS
    write(char,1000) s
    filename = 'Persistence'
    filename = filename(1:11) // char // '_2005.csv'
    open(unit = 1, file = filename)
        read(1,*)
        read(1,*) sigmaPers2005Data(s), (xsi2005Data(s,i), i = 1, NREG_Pers)
        do i = 1, NREG_Pers
            read(1,*) blah, (COVxsi2005Data(s,i,j), j = 1, NREG_Pers)
        end do
    close(1)
    
    filename = 'cp' // filename
    open(unit = 2, file = filename)
        read(2,*)
        do i = 1, NREG_Pers
            read(2,*) NobsPers2005Data(s), blah2, (CPPers2005Data(s,i,j), j = 1, NREG_Pers)
        end do
    close(2)
    
    
    ! Check if CPPers2005Data is positive definite
        
    D = CPPers2005Data(s,:,:)
        
    call dpotrf('L', NREG_Pers, D, NREG_Pers, info)
    if (info .ne. 0) then
        print*, 'CPPers2005Data is not symmetric positive definite.'
        print*, 'Aborting program.'
        pause
        stop
    end if
        
    UPers = COVxsi2005Data(s,:,:)
    
    call dpotrf('U',NREG_Pers,UPers,NREG_Pers,info)    
    
    invCOVxsi2005Data(s,:,:) = UPers
    
    call dpotri('U',NREG_Pers,invCOVxsi2005Data(s,:,:),NREG_Pers,info)
    
    do i = 1, NREG_Pers
        do j = 1, NREG_Pers
            if (i > j) then
                invCOVxsi2005Data(s,i,j) = invCOVxsi2005Data(s,j,i)
            end if
        end do
    end do
    
end do 

! ***************************************************************
! ***************************************************************


! *********************
! Frequency Regressions
! *********************

do s = 0, NSECTORS
    write(char,1000) s
    filename = 'freq'
    filename = filename(1:4) // char // '.csv'
    open(unit = 1, file = filename)
        read(1,*)
        read(1,*) sigmaFreqData(s), (etaData(s,i), i = 1, NREG_Freq)
        do i = 1, NREG_Freq
            read(1,*) blah, (COVetaData(s,i,j), j = 1, NREG_Freq)
        end do
    close(1)
    
    filename = 'cp' // filename
    open(unit = 2, file = filename)
        read(2,*)
        do i = 1, NREG_Freq
            read(2,*) NobsFreqData(s), blah2, (CPFreqData(s,i,j), j = 1, NREG_Freq)
        end do
    close(2)
    
    
    ! Check if CPPers2005Data is positive definite
        
    E = CPFreqData(s,:,:)
        
    call dpotrf('L', NREG_Freq, E, NREG_Freq, info)
    if (info .ne. 0) then
        print*, 'CPFreqData is not symmetric positive definite.'
        print*, 'Aborting program.'
        pause
        stop
    end if
        
    UFreq = COVetaData(s,:,:)
    
    call dpotrf('U',NREG_Freq,UFreq,NREG_Freq,info)    
    
    invCOVetaData(s,:,:) = UFreq
    
    call dpotri('U',NREG_Freq,invCOVetaData(s,:,:),NREG_Freq,info)
    
    do i = 1, NREG_Freq
        do j = 1, NREG_Freq
            if (i > j) then
                invCOVetaData(s,i,j) = invCOVetaData(s,j,i)
            end if
        end do
    end do
    
end do    

! ***************************************************************
! ***************************************************************


1000 format(i1)

end subroutine Read_Coef

end program Main





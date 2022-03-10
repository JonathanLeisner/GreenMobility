program Main

! ***************************************************************
! This program simulates a trade shock in High-Tech Manufacturing
! Different assumptions on the mobility of physical capital are
! imposed: perfect mobility (Simulation1), no capital mobility
! (Simulation2) and imperfect capital mobility (Simulation3)
! ***************************************************************



! ***********************
! Including other modules
! ***********************
USE MKL_VSL_TYPE
USE MKL_VSL
USE Global_Data
USE Emax_MOD
USE Simulation_MOD1
!USE Simulation_MOD2
!USE Simulation_MOD3
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


! Random Number Generation
TYPE (VSL_STREAM_STATE) :: stream
integer(kind=4)   errcode
integer           brng, method, method1, method2, method3, method4, seed
real(KIND=DOUBLE) rvec_eps(NSECTORS+1), rvec_eta(NSECTORS+1)

! MPI Variables
integer ierr, size, rank, tag, source, count
integer status(MPI_STATUS_SIZE)


! *******************
! Program Starts Here
! *******************

Call MPI_INIT(ierr)
Call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
Call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

! *******************
! Program Starts Here
! *******************

if (rank == 0) then

    ! ************************
    ! Read Data for Simulation
    ! ************************

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
            aaGlobal(s,0,t) = correction*sum(TotalHrWageData(s,t,1:2)) / OutputData(s,t)
            aaGlobal(s,1,t) = correction*sum(TotalHrWageData(s,t,3:4)) / OutputData(s,t)
        end do
    end do

    open(unit = 1, file = 'Labor_Shares.csv')
    write(1,9876)
    do ed = 0, 1
        do t = FIRST_YR, LAST_YR
           write(1,8765) t, ed, (aaGlobal(s,ed,t), s = 1, NSECTORS)
        end do
    end do
    close(1)

    9876 format('Year, Educ Level, Agr/Min, LT Manuf, HT Manuf, Const, Trade, Trans, Service')
    8765 format(2(i5,','),7(f7.4,','))


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


    param(1:7)      = theta(1:7)       
    param(8:19)     = beta(1,1:12)     
    param(20:31)    = beta(2,1:12)     
    param(32:43)    = beta(3,1:12)     
    param(44:55)    = beta(4,1:12)     
    param(56:67)    = beta(5,1:12)     
    param(68:79)    = beta(6,1:12)     
    param(80:91)    = beta(7,1:12)     
    param(92:99)    = sigma(0:7)
    param(100:125)  = kappa(1:26)      
    param(126:131)  = tau(2:7)         
    param(132)      = SigmaPref
    param(133:140)  = omega(2,0:7)     
    param(141:148)  = omega(3,0:7)     
    param(149:150)  = lambda(2:3)      
    param(151:162)  = gamma(2,0:11)    
    param(163:174)  = gamma(3,0:11)    
    param(175)      = aa(1,0,0)
    param(176)      = aa(1,0,1)        
    param(177)      = aa(2,0,0)
    param(178)      = aa(2,0,1)        
    param(179)      = aa(3,0,0)
    param(180)      = aa(3,0,1)        
    param(181)      = aa(4,0,0)
    param(182)      = aa(4,0,1)        
    param(183)      = aa(5,0,0)
    param(184)      = aa(5,0,1)        
    param(185)      = aa(6,0,0)
    param(186)      = aa(6,0,1)        
    param(187)      = aa(7,0,0)
    param(188)      = aa(7,0,1)        
    param(189)      = aa(1,1,0)
    param(190)      = aa(1,1,1)        
    param(191)      = aa(2,1,0)
    param(192)      = aa(2,1,1)        
    param(193)      = aa(3,1,0)
    param(194)      = aa(3,1,1)        
    param(195)      = aa(4,1,0)
    param(196)      = aa(4,1,1)        
    param(197)      = aa(5,1,0)
    param(198)      = aa(5,1,1)        
    param(199)      = aa(6,1,0)
    param(200)      = aa(6,1,1)        
    param(201)      = aa(7,1,0)
    param(202)      = aa(7,1,1)        
    param(203:209)  = sigma_prod(1:7)
    param(210:216)  = rsk_init(1:7,0)
    param(217:223)  = rsk_init(1:7,1)
    
end if

Call MPI_BCAST(param,NPARAM,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)    

! Perfect Capital Mobility
Call Simulation1(param)

! No Capital Mobility
! Call Simulation2(param)

! Imperfect Capital Mobility
! Call Simulation3(param)

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

end program Main





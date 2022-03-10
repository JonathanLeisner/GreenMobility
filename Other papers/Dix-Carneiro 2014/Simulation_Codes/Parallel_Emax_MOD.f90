include 'mkl_vsl.fi'
include 'lapack.f90'


MODULE Emax_MOD

! **********************************************************
! This module contains the subroutine that approximates
! The Emax function
! it contains:
subroutine EmaxCoef_Myopic(param,cut,gender,educ,tp,PI_Myopic,R_sq_Myopic)

! **********************************************************

USE Global_Data
USE MKL_VSL_TYPE
USE MKL_VSL
USE LinReg_MOD


Contains


! Computes Emax imposing static expectations (first iteration)
subroutine EmaxCoef_Myopic(param,cut,gender,educ,tp,PI_Myopic,R_sq_Myopic)


implicit none

! *********************************************************************************
! *********************************************************************************

! ************************************
! Inputs and Outputs of the Subroutine
! ************************************

real(KIND=DOUBLE), intent(in)   :: param(:), cut(:,0:)
integer, intent(in)             :: gender, educ, tp
real(KIND=DOUBLE), intent(out)  :: PI_Myopic(FIRST_AGE+1:,0:,:), &
                                   R_sq_Myopic(FIRST_AGE+1:,0:)

real(KIND=DOUBLE) PI_COEF(NREG)

! ******************
! Counting Variables
! ******************

integer   i, ii, j, t, k, kk, l, n, s, s1, s2


! ***********************
! Parameters of the model
! ***********************

real(KIND=DOUBLE) theta(7), beta(NSECTORS,12), & 
                  sigma(0:NSECTORS), kappa(26), tau(0:NSECTORS), SigmaPref, &
                  omega(2:NTYPES,0:NSECTORS), lambda(NTYPES)


! ***********************
! Wage variables and Emax
! ***********************

real(KIND=DOUBLE) w(0:NSECTORS), hc(NSECTORS), V(0:NSECTORS), VMAX, VM, & 
                  Emax(0:NSECTORS), cost(0:NSECTORS,0:NSECTORS), dummy(NSECTORS), &
                  SumExpV, Emax0(0:NSECTORS)
                  
integer           exper(NSECTORS), lag, age, &
                  educ2, educ3, educ4, ed, ExperTomorrow(NSECTORS)


! *******************************     
! Variables for the approximation
! *******************************
                 
real(KIND=DOUBLE) COEF(NREG,1), SSE, SST, MeanEmax
real(KIND=DOUBLE) rvec_eps(0:NSECTORS), rvec_rsk(NSECTORS)
real(KIND=DOUBLE), allocatable, dimension(:,:,:) :: eps
real(KIND=DOUBLE) rsk(NSECTORS)
integer           ch(9)

real(KIND=DOUBLE), allocatable, dimension(:,:)   :: Y
real(KIND=DOUBLE), allocatable, dimension(:,:)   :: X


! ************************
! Random Number Generation
! ************************

TYPE (VSL_STREAM_STATE) :: stream1, stream2, stream3
integer(KIND=4) errcode1, errcode2, errcode3
integer brng1, brng2, brng3, method1, method2, method3, seed1, seed2, seed3
integer info 


! *********************************************************************************
! *********************************************************************************


allocate(eps(FIRST_AGE+1:LAST_AGE,NDRAWS,0:NSECTORS))
allocate(X(INTP,NREG))
allocate(Y(INTP,1))

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



! *******************************************************************************
! Generating Random Draws for the Monte Carlo Integration
! We generate NDRAWS for each age, in order to save computational time
! Just as a reminder, in the implementation of Keane and Wolpin (1994),
! Drawing NDRAWS just once and for all the states and all ages did not perform
! very well for approximating Emax. Drawing the same points for all states by age
! performed much better
! *******************************************************************************

! *******************************
! Random Generator Initialization 
! for HC Shocks
! *******************************

brng1    = VSL_BRNG_MT19937
method1  = VSL_METHOD_DGAUSSIAN_BOXMULLER
seed1    = 12345678
errcode1 = vslnewstream(stream1,brng1,seed1)

! ***** Warming Up *****
do j = 1, WARMUP
    errcode1 = vdrnggaussian(method1, stream1, NSECTORS+1, &
               rvec_eps, real(0.0,DOUBLE), real(1.0,DOUBLE))
end do                        

do age = FIRST_AGE+1, LAST_AGE
    do j = 1, NDRAWS
        errcode1 = vdrnggaussian(method1, stream1, NSECTORS+1, &
                   rvec_eps, real(0.0,DOUBLE), real(1.0,DOUBLE))
        do i = 0, NSECTORS
            eps(age,j,i) = sigma(i)*rvec_eps(i)
        end do
    end do
end do

! ***********************************
! Deinitialize Stream for Wage Shocks
! ***********************************

errcode1 = vsldeletestream(stream1)




! ********************************
! Random Generator Initialization 
! for Returns to Skill And Choices
! ********************************

! ***** Initializing *****
! HC rpices
brng2    = VSL_BRNG_MT19937
seed2    = 98765456
method2 = VSL_METHOD_IUNIFORM_STD
errcode2 = vslnewstream(stream2,brng2,seed2)
! Choices
brng3    = VSL_BRNG_MT19937
seed3    = 25109825
method3 = VSL_METHOD_IBERNOULLI_ICDF
errcode3 = vslnewstream(stream3,brng3,seed3)



! ***** Warming Up *****
do j = 1, WARMUP
    errcode2 = vdrnguniform(method2, stream2, NSECTORS, rsk, &
              real(0.0,DOUBLE), real(1.0,DOUBLE))
    errcode3 = virnguniform(method3, stream3, 9, ch, 0, NSECTORS+1)
end do  

! Skilled / Unskilled Worker
if (educ <= 2) then
    ed = 0
else if (educ >= 3) then
    ed = 1
end if


PI_Myopic   = 0.0
R_sq_Myopic = 0.0


age_loop: do age = LAST_AGE, (FIRST_AGE+1), -1

    lag_loop: do lag = 0, NSECTORS

        INTP_loop: do kk = 1, INTP
        
    
            ! *****************************************************************
            ! Random Number Generation - returns to skill and past choices
            ! *****************************************************************
           
            errcode2 = vdrnguniform(method2, stream2, NSECTORS, rvec_rsk, &
                       real(0.0,DOUBLE), real(1.0,DOUBLE))
        
            rsk = LOW_RSK(:,ed) + (UP_RSK(:,ed) - LOW_RSK(:,ed))*rvec_rsk
          
            errcode3 = virnguniform(method3, stream3, 8, ch(2:9), 0, NSECTORS+1) 
        
            ! *****************************************************************
            ! *****************************************************************
        
                
            ch(1) = lag
            
            
            ! experience accumulated in the last 8 years
            exper = 0
            do s = 1, NSECTORS
                do i = 1, 8
                    if (ch(i) == s) then
                        exper(s) = exper(s) + 1
                    end if
                end do
            end do
            
            dummy = 0.0
            do s = 1, NSECTORS
                if (rsk(s) >= cut(s,ed)) then
                    dummy(s) = 1
                end if
            end do
            
            if (age < LAST_AGE) then
               
                do s = 0, NSECTORS
                    
                    ExperTomorrow = exper
                    if (s >= 1) then
                        ExperTomorrow(s) = exper(s) + 1
                    end if
                    
                    PI_COEF = PI_Myopic(age+1,s,:)
                    Emax(s) = Emax_hat(PI_COEF, rsk, ExperTomorrow, cut(:,ed))
                    
                end do
                    
            else
                Emax = 0.0            
            end if
            
            
            ! experience accumulated in the last 9 years
            do s = 1, NSECTORS
                if (ch(9) == s) then
                    exper(s) = exper(s) + 1
                end if
            end do
        
        
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
                 
            if (tp == 1) then
                w(0) = exp(theta(1) + theta(2)*(gender-1) + theta(3)*educ2 + &
                theta(4)*educ3 + theta(5)*educ4 + theta(6)*(age-25) + & 
                theta(7)*(age-25)**2)     
            else if (tp == 2 .or. tp == 3) then
                w(0) = exp(omega(tp,0) + theta(1) + theta(2)*(gender-1) + theta(3)*educ2 + &
                theta(4)*educ3 + theta(5)*educ4 + theta(6)*(age-25) + & 
                theta(7)*(age-25)**2) 
            end if
            
            do s = 1, NSECTORS
                if (tp == 1) then
                    hc(s) = exp(beta(s,1)*(gender-1) + beta(s,2)*educ2 + & 
                    beta(s,3)*educ4 + beta(s,4)*(age-25) + beta(s,5)*(age-25)**2 + &
                    beta(s,6)*exper(1) + beta(s,7)*exper(2) + beta(s,8)*exper(3) + & 
                    beta(s,9)*exper(4) + beta(s,10)*exper(5) + beta(s,11)*exper(6) + beta(s,12)*exper(7))               
                else if (tp == 2 .or. tp == 3) then
                    hc(s) = exp(omega(tp,s) + beta(s,1)*(gender-1) + beta(s,2)*educ2 + & 
                    beta(s,3)*educ4 + beta(s,4)*(age-25) + beta(s,5)*(age-25)**2 + &
                    beta(s,6)*exper(1) + beta(s,7)*exper(2) + beta(s,8)*exper(3) + & 
                    beta(s,9)*exper(4) + beta(s,10)*exper(5) + beta(s,11)*exper(6) + beta(s,12)*exper(7))
                end if
            end do
            
            cost = 0.0
            
            do s1 = 0, NSECTORS
                do s2 = 1, NSECTORS
                    if (s1 == 0) then
                        cost(s1,s2) = exp(kappa(s2) + kappa(21)*(gender-1) + kappa(22)*educ2 + & 
                        kappa(23)*educ3 + kappa(24)*educ4 + kappa(25)*(age-25) + & 
                        kappa(26)*(age-25)**2 + lambda(tp))
                    else if (s1 == 1 .and. (s1 /= s2)) then
                        cost(s1,s2) = exp(kappa(13+s2) + kappa(21)*(gender-1) + kappa(22)*educ2 + & 
                        kappa(23)*educ3 + kappa(24)*educ4 + kappa(25)*(age-25) + & 
                        kappa(26)*(age-25)**2 + lambda(tp))
                    else if (s1 /= 0 .and. s1 /= 1 .and. (s1 /= s2)) then
                        cost(s1,s2) = exp(kappa(6+s1) + kappa(13+s2) + kappa(21)*(gender-1) + kappa(22)*educ2 + & 
                        kappa(23)*educ3 + kappa(24)*educ4 + kappa(25)*(age-25) + & 
                        kappa(26)*(age-25)**2 + lambda(tp))
                    end if                 
                end do
            end do            
             
                
            MeanEmax = 0.0
            do n = 1, NDRAWS
                
                V(0) = tau(0) + w(0) + eps(age,n,0) + rho*Emax(0)
                do s = 1, NSECTORS
                    V(s) = tau(s) + rsk(s)*hc(s)*exp(eps(age,n,s)) - cost(lag,s) + rho*Emax(s)
                end do

                VM = maxval(V)
                SumExpV = 0.0
                do s = 0, NSECTORS
                    SumExpV = SumExpV + exp((V(s)-VM)/sigmaPref)
                end do
                VMAX = VM + sigmaPref*log(SumExpV)

                MeanEmax = MeanEmax + VMAX

            end do
            MeanEmax = MeanEmax / NDRAWS
           
            
           
            ! We approximate Emax with a 2nd order polynomial
            ! In the state variables, in each of the regions
            ! Defined by the vector cut
            
           
            Y(kk,1)       = MeanEmax
        
            X(kk,1)       = 1.0
            
            X(kk,2:8)     = rsk(1:7)
            
            X(kk,9:15)    = dummy(1:7)*(rsk(1:7)-cut(1:7,ed))
            
            X(kk,16:22)   = exper(1:7)
                        
            X(kk,23:29)   = rsk(1:7)**2
            
            X(kk,30:36)   = dummy(1:7)*(rsk(1:7)-cut(1:7,ed))**2
            
            X(kk,37:43)   = exper(1:7)**2
            
            X(kk,44:49)   = rsk(1)*rsk(2:7)
           
            X(kk,50:54)   = rsk(2)*rsk(3:7)
          
            X(kk,55:58)   = rsk(3)*rsk(4:7)
            
            X(kk,59:61) = rsk(4)*rsk(5:7)
            
            X(kk,62:63) = rsk(5)*rsk(6:7)
            
            X(kk,64)     = rsk(6)*rsk(7)
            
            X(kk,65:71) = rsk(1)*exper(1:7)
                        
            X(kk,72:78) = rsk(2)*exper(1:7)
            
            X(kk,79:85) = rsk(3)*exper(1:7)
            
            X(kk,86:92) = rsk(4)*exper(1:7)
            
            X(kk,93:99) = rsk(5)*exper(1:7)
            
            X(kk,100:106) = rsk(6)*exper(1:7)
            
            X(kk,107:113) = rsk(7)*exper(1:7)
            
            X(kk,114:119) = exper(1)*exper(2:7)
        
            X(kk,120:124) = exper(2)*exper(3:7)
        
            X(kk,125:128) = exper(3)*exper(4:7)
            
            X(kk,129:131) = exper(4)*exper(5:7)
            
            X(kk,132:133) = exper(5)*exper(6:7)
            
            X(kk,134)     = exper(6)*exper(7)

        end do INTP_loop              
    
        ! ***** Least Squares Fitting ****
        ! ********************************
    
         Call LinReg(Y,X,INTP,NREG,COEF,SST,SSE,info)
             
        ! ********************************
        ! ********************************
    
        PI_Myopic(age,lag,:) = COEF(:,1)
        R_sq_Myopic(age,lag) = 1 - (SSE/SST)    
        
    end do lag_loop

end do age_loop

364 format(6(i6,','),2(f20.8,','))

! *********************************
! Deinitialize Stream for 
! returns to skill and past choices
! *********************************

errcode2 = vsldeletestream(stream2)
errcode3 = vsldeletestream(stream3)

deallocate(Y,X,eps)

end subroutine EmaxCoef_Myopic




! *********************************************************************************************************************************
! *********************************************************************************************************************************
! *********************************************************************************************************************************
! *********************************************************************************************************************************
! *********************************************************************************************************************************
! *********************************************************************************************************************************
! *********************************************************************************************************************************
! *********************************************************************************************************************************
! *********************************************************************************************************************************
! *********************************************************************************************************************************
! *********************************************************************************************************************************
! *********************************************************************************************************************************


! Computes Emax under Rational Expectations
! Notice that there is no aggregate uncertainty in the simulation.
! Therefore, the Rational Expectations equilibrium is a perfect foresight equilibrium
! For the algorithm used in computing the perfect foresight equilibrium, refer to:
! "An Estimable Dynamic General Equilibrium Model of Work, Schooling and Occupational Choice," by Donghoon Lee
! International Economic Review, pg. 46, February, 2005
subroutine EmaxCoef_RE(param,avec,cut,gender,educ,tp,PI_RE,R_sq_RE,PI_Myopic)

implicit none


real(KIND=DOUBLE), intent(in)   :: param(:), cut(:,0:), avec(:,0:,FIRST_YR_SIM:)
integer, intent(in)             :: gender, educ, tp
real(KIND=DOUBLE), intent(out)  :: PI_RE(YR_ANNMNT+1:,FIRST_AGE+1:,0:,:), &
                                   R_sq_RE(YR_ANNMNT+1:,FIRST_AGE+1:,0:)
real(KIND=DOUBLE), intent(in)  ::  PI_Myopic(FIRST_AGE+1:,0:,:,:,:,:)

real(KIND=DOUBLE) PI_COEF(NREG), PI_COEF1(NREG), PI_COEF2(NREG)

! ******************
! Counting Variables
! ******************

integer   i, ii, j, t, k, kk, l, n, s, s1, s2, year


! ***********************
! Parameters of the model
! ***********************

real(KIND=DOUBLE) theta(7), beta(NSECTORS,12), & 
                  sigma(0:NSECTORS), kappa(26), tau(0:NSECTORS), SigmaPref, &
                  omega(2:NTYPES,0:NSECTORS), lambda(NTYPES)


! ***********************
! Wage variables and Emax
! ***********************

real(KIND=DOUBLE) w(0:NSECTORS), hc(NSECTORS), V(0:NSECTORS), VMAX, VM, & 
                  Emax(0:NSECTORS), cost(0:NSECTORS,1:NSECTORS), dummy(NSECTORS), &
                  SumExpV
                  
integer           exper(NSECTORS), lag, age, &
                  educ2, educ3, educ4, ed, ExperTomorrow(NSECTORS)


! *******************************     
! Variables for the approximation
! *******************************
                 
real(KIND=DOUBLE) COEF(NREG,1), SSE, SST, MeanEmax
real(KIND=DOUBLE) rvec_eps(0:NSECTORS), rvec_rsk(NSECTORS)
real(KIND=DOUBLE), allocatable, dimension(:,:,:) :: eps
real(KIND=DOUBLE) rsk(NSECTORS), rsk_tomorrow(NSECTORS)
integer           ch(9)

real(KIND=DOUBLE), allocatable, dimension(:,:)   :: Y
real(KIND=DOUBLE), allocatable, dimension(:,:)   :: X


TYPE (VSL_STREAM_STATE) :: stream1, stream2, stream3
integer(KIND=4) errcode1, errcode2, errcode3
integer brng1, brng2, brng3, method1, method2, method3, seed1, seed2, seed3
integer info 


! *********************************************************************************
! *********************************************************************************


allocate(eps(FIRST_AGE+1:LAST_AGE,NDRAWS,0:NSECTORS))
allocate(X(INTP,NREG))
allocate(Y(INTP,1))

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



PI_RE   = 0.0
R_sq_RE = 0.0


year_loop: do year = (LAST_YR_SIM-1), (YR_ANNMNT+1), -1

    ! *******************************************************************************
    ! Generating Random Draws for the Monte Carlo Integration
    ! We generate NDRAWS for each age, in order to save computational time
    ! Just as a reminder, in the implementation of Keane and Wolpin (1994),
    ! Drawing NDRAWS just once and for all the states and all ages did not perform
    ! very well for approximating Emax. Drawing the same points for all states by age
    ! performed much better
    ! *******************************************************************************

    ! *******************************
    ! Random Generator Initialization 
    ! for HC Shocks
    ! *******************************

    brng1    = VSL_BRNG_MT19937
    method1  = VSL_METHOD_DGAUSSIAN_BOXMULLER
    seed1    = 12345678
    errcode1 = vslnewstream(stream1,brng1,seed1)

    ! ***** Warming Up *****
    do j = 1, WARMUP
        errcode1 = vdrnggaussian(method1, stream1, NSECTORS+1, &
                   rvec_eps, real(0.0,DOUBLE), real(1.0,DOUBLE))
    end do                        

    do age = FIRST_AGE+1, LAST_AGE
        do j = 1, NDRAWS
            errcode1 = vdrnggaussian(method1, stream1, NSECTORS+1, &
                       rvec_eps, real(0.0,DOUBLE), real(1.0,DOUBLE))
            do i = 0, NSECTORS
                eps(age,j,i) = sigma(i)*rvec_eps(i)
            end do
        end do
    end do

    ! ***********************************
    ! Deinitialize Stream for Wage Shocks
    ! ***********************************

    errcode1 = vsldeletestream(stream1)




    ! ********************************
    ! Random Generator Initialization 
    ! for Returns to Skill And Choices
    ! ********************************

    ! ***** Initializing *****
    ! HC rpices
    brng2    = VSL_BRNG_MT19937
    seed2    = 98765456
    method2 = VSL_METHOD_IUNIFORM_STD
    errcode2 = vslnewstream(stream2,brng2,seed2)
    ! Choices
    brng3    = VSL_BRNG_MT19937
    seed3    = 25109825
    method3 = VSL_METHOD_IBERNOULLI_ICDF
    errcode3 = vslnewstream(stream3,brng3,seed3)
    
    ! ***** Warming Up *****
    do j = 1, WARMUP
        errcode2 = vdrnguniform(method2, stream2, NSECTORS, rsk, &
                  real(0.0,DOUBLE), real(1.0,DOUBLE))
        errcode3 = virnguniform(method3, stream3, 9, ch, 0, NSECTORS+1)
    end do  

    
    
    age_loop: do age = LAST_AGE, (FIRST_AGE+1), -1
    
        lag_loop: do lag = 0, NSECTORS

            INTP_loop: do kk = 1, INTP
        
                if (educ <= 2) then
                    ed = 0
                else if (educ >= 3) then
                    ed = 1
                end if
        
    
                ! *****************************************************************
                ! Random Number Generation - returns to skill and past choices
                ! *****************************************************************
           
                errcode2 = vdrnguniform(method2, stream2, NSECTORS, rvec_rsk, &
                           real(0.0,DOUBLE), real(1.0,DOUBLE))
        
                rsk = LOW_RSK(:,ed) + (UP_RSK(:,ed) - LOW_RSK(:,ed))*rvec_rsk
          
                errcode3 = virnguniform(method3, stream3, 8, ch(2:9), 0, NSECTORS+1) 
        
                ! *****************************************************************
                ! *****************************************************************
        
                
                ch(1) = lag
            
            
                ! experience accumulated in the last 8 years
                exper = 0
                do s = 1, NSECTORS
                    do i = 1, 8
                        if (ch(i) == s) then
                            exper(s) = exper(s) + 1
                        end if
                    end do
                end do
               
                           
                if (age < LAST_AGE) then
                
                    rsk_tomorrow = avec(:,ed,year+1)*rsk

		            do s = 0, NSECTORS
                        
                        ExperTomorrow = exper
                        if (s >= 1) then
                            ExperTomorrow(s) = exper(s) + 1
                        end if
                        
                        if (year <= (LAST_YR_SIM-2)) then  
                            PI_COEF = PI_RE(year+1,age+1,s,:)
                            PI_COEF1 = PI_Myopic(age+1,s,gender,educ,tp,:)
                            PI_COEF2 = PI_RE(year+1,age+1,s,:)
                        else
                            PI_COEF = PI_Myopic(age+1,s,gender,educ,tp,:)
                        end if
                        
                        Emax(s) = Emax_hat(PI_COEF, rsk_tomorrow, ExperTomorrow, cut(:,ed))
                    
                    end do                    
                    
                    
                else
                    Emax = 0.0            
                end if
            
            
               ! experience accumulated in the last 9 years
                do s = 1, NSECTORS
                    if (ch(9) == s) then
                        exper(s) = exper(s) + 1
                    end if
                end do        
        
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
                 
                if (tp == 1) then
                    w(0) = exp(theta(1) + theta(2)*(gender-1) + theta(3)*educ2 + &
                    theta(4)*educ3 + theta(5)*educ4 + theta(6)*(age-25) + & 
                    theta(7)*(age-25)**2)     
                else if (tp == 2 .or. tp == 3) then
                    w(0) = exp(omega(tp,0) + theta(1) + theta(2)*(gender-1) + theta(3)*educ2 + &
                    theta(4)*educ3 + theta(5)*educ4 + theta(6)*(age-25) + & 
                    theta(7)*(age-25)**2) 
                end if
            
                do s = 1, NSECTORS
                    if (tp == 1) then
                        hc(s) = exp(beta(s,1)*(gender-1) + beta(s,2)*educ2 + & 
                        beta(s,3)*educ4 + beta(s,4)*(age-25) + beta(s,5)*(age-25)**2 + &
                        beta(s,6)*exper(1) + beta(s,7)*exper(2) + beta(s,8)*exper(3) + & 
                        beta(s,9)*exper(4) + beta(s,10)*exper(5) + beta(s,11)*exper(6) + beta(s,12)*exper(7))               
                    else if (tp == 2 .or. tp == 3) then
                        hc(s) = exp(omega(tp,s) + beta(s,1)*(gender-1) + beta(s,2)*educ2 + & 
                        beta(s,3)*educ4 + beta(s,4)*(age-25) + beta(s,5)*(age-25)**2 + &
                        beta(s,6)*exper(1) + beta(s,7)*exper(2) + beta(s,8)*exper(3) + & 
                        beta(s,9)*exper(4) + beta(s,10)*exper(5) + beta(s,11)*exper(6) + beta(s,12)*exper(7))
                    end if
                end do
            
                cost = 0.0
            
                do s1 = 0, NSECTORS
                    do s2 = 1, NSECTORS
                        if (s1 == 0) then
                            cost(s1,s2) = exp(kappa(s2) + kappa(21)*(gender-1) + kappa(22)*educ2 + & 
                            kappa(23)*educ3 + kappa(24)*educ4 + kappa(25)*(age-25) + & 
                            kappa(26)*(age-25)**2 + lambda(tp))
                        else if (s1 == 1 .and. (s1 /= s2)) then
                            cost(s1,s2) = exp(kappa(13+s2) + kappa(21)*(gender-1) + kappa(22)*educ2 + & 
                            kappa(23)*educ3 + kappa(24)*educ4 + kappa(25)*(age-25) + & 
                            kappa(26)*(age-25)**2 + lambda(tp))
                        else if (s1 /= 0 .and. s1 /= 1 .and. (s1 /= s2)) then
                            cost(s1,s2) = exp(kappa(6+s1) + kappa(13+s2) + kappa(21)*(gender-1) + kappa(22)*educ2 + & 
                            kappa(23)*educ3 + kappa(24)*educ4 + kappa(25)*(age-25) + & 
                            kappa(26)*(age-25)**2 + lambda(tp))
                        end if                 
                    end do
                end do            
             
                
                MeanEmax = 0.0
                do n = 1, NDRAWS
                
                    V(0) = tau(0) + w(0) + eps(age,n,0) + rho*Emax(0)
                    do s = 1, NSECTORS
                        V(s) = tau(s) + rsk(s)*hc(s)*exp(eps(age,n,s)) - cost(lag,s) + rho*Emax(s)
                    end do

                    VM = maxval(V)
                    SumExpV = 0.0
                    do s = 0, NSECTORS
                        SumExpV = SumExpV + exp((V(s)-VM)/sigmaPref)
                    end do
                    VMAX = VM + sigmaPref*log(SumExpV)

                    MeanEmax = MeanEmax + VMAX

                end do
                MeanEmax = MeanEmax / NDRAWS
           
            
                dummy = 0.0
                do s = 1, NSECTORS
                    if (rsk(s) >= cut(s,ed)) then
                        dummy(s) = 1
                    end if
                end do
            
                ! We approximate Emax with a 2nd order polynomial
                ! In the state variables, in each of the regions
                ! Defined by the vector cut
            
           
                Y(kk,1)       = MeanEmax
        
                X(kk,1)       = 1.0
            
                X(kk,2:8)     = rsk(1:7)
            
                X(kk,9:15)    = dummy(1:7)*(rsk(1:7)-cut(1:7,ed))
            
                X(kk,16:22)   = exper(1:7)
                        
                X(kk,23:29)   = rsk(1:7)**2
            
                X(kk,30:36)   = dummy(1:7)*(rsk(1:7)-cut(1:7,ed))**2
            
                X(kk,37:43)   = exper(1:7)**2
            
                X(kk,44:49)   = rsk(1)*rsk(2:7)
           
                X(kk,50:54)   = rsk(2)*rsk(3:7)
          
                X(kk,55:58)   = rsk(3)*rsk(4:7)
            
                X(kk,59:61) = rsk(4)*rsk(5:7)
            
                X(kk,62:63) = rsk(5)*rsk(6:7)
            
                X(kk,64)     = rsk(6)*rsk(7)
            
                X(kk,65:71) = rsk(1)*exper(1:7)
                        
                X(kk,72:78) = rsk(2)*exper(1:7)
            
                X(kk,79:85) = rsk(3)*exper(1:7)
            
                X(kk,86:92) = rsk(4)*exper(1:7)
            
                X(kk,93:99) = rsk(5)*exper(1:7)
            
                X(kk,100:106) = rsk(6)*exper(1:7)
            
                X(kk,107:113) = rsk(7)*exper(1:7)
            
                X(kk,114:119) = exper(1)*exper(2:7)
        
                X(kk,120:124) = exper(2)*exper(3:7)
        
                X(kk,125:128) = exper(3)*exper(4:7)
            
                X(kk,129:131) = exper(4)*exper(5:7)
            
                X(kk,132:133) = exper(5)*exper(6:7)
            
                X(kk,134)     = exper(6)*exper(7)

            end do INTP_loop
            
            ! ***** Least Squares Fitting ****
            ! ********************************
    
            Call LinReg(Y,X,INTP,NREG,COEF,SST,SSE,info)
             
            ! ********************************
            ! ********************************
    
            PI_RE(year,age,lag,:) = COEF(:,1)
            R_sq_RE(year,age,lag) = 1 - (SSE/SST)
            
        end do lag_loop

    end do age_loop
end do year_loop    

end subroutine EmaxCoef_RE




function Emax_hat(PI, rsk, exper, cut)


!*******************************************************
! Calculates the approximate value for Emax at the state 
! described by rsk, exper, cut and given the
! regression coefficients PI
!*******************************************************

implicit none

integer          , intent(in) :: exper(:)
real(KIND=DOUBLE), intent(in) :: PI(:), rsk(:), cut(:)

integer s, i

real(KIND=DOUBLE) XR(NREG), Emax_hat

real(KIND=DOUBLE) dummy(NSECTORS)


dummy = 0.0
do s = 1, NSECTORS
    if (rsk(s) >= cut(s)) then
        dummy(s) = 1
    end if
end do
        
XR(1)       = 1.0
            
XR(2:8)     = rsk(1:7)
            
XR(9:15)    = dummy(1:7)*(rsk(1:7)-cut(1:7))
            
XR(16:22)   = exper(1:7)
                        
XR(23:29)   = rsk(1:7)**2
            
XR(30:36)   = dummy(1:7)*(rsk(1:7)-cut(1:7))**2
            
XR(37:43)   = exper(1:7)**2
            
XR(44:49)   = rsk(1)*rsk(2:7)
           
XR(50:54)   = rsk(2)*rsk(3:7)
          
XR(55:58)   = rsk(3)*rsk(4:7)
            
XR(59:61) = rsk(4)*rsk(5:7)
            
XR(62:63) = rsk(5)*rsk(6:7)
            
XR(64)     = rsk(6)*rsk(7)
            
XR(65:71) = rsk(1)*exper(1:7)
                        
XR(72:78) = rsk(2)*exper(1:7)
            
XR(79:85) = rsk(3)*exper(1:7)
            
XR(86:92) = rsk(4)*exper(1:7)
            
XR(93:99) = rsk(5)*exper(1:7)
            
XR(100:106) = rsk(6)*exper(1:7)
            
XR(107:113) = rsk(7)*exper(1:7)
            
XR(114:119) = exper(1)*exper(2:7)
        
XR(120:124) = exper(2)*exper(3:7)
        
XR(125:128) = exper(3)*exper(4:7)
            
XR(129:131) = exper(4)*exper(5:7)
            
XR(132:133) = exper(5)*exper(6:7)
            
XR(134)     = exper(6)*exper(7)

Emax_hat = sum(PI*XR)

end function Emax_hat


end module Emax_MOD

include 'mkl_vsl.fi'
include 'lapack.f90'


MODULE Emax_MOD

! **********************************************************
! This module contains the subroutine that approximates
! The Emax function
! it contains:
! subroutine Calc_Emax(param, cut, PI, R_sq)
! subroutine EmaxCoef(param,cut,gender,educ,PI,R_sq)
! function   Emax_hat(PI, rsk, exper, cut)
! **********************************************************

USE Global_Data
USE MKL_VSL_TYPE
USE MKL_VSL
USE LinReg_MOD


Contains



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
                    
                    Emax(s) = &
                    PI_Myopic(age+1,s,1)*1.0 + &
                    PI_Myopic(age+1,s,2)*rsk(1) + &
                    PI_Myopic(age+1,s,3)*rsk(2) + &
                    PI_Myopic(age+1,s,4)*rsk(3) + &
                    PI_Myopic(age+1,s,5)*rsk(4) + &
                    PI_Myopic(age+1,s,6)*rsk(5) + &
                    PI_Myopic(age+1,s,7)*rsk(6) + &
                    PI_Myopic(age+1,s,8)*rsk(7) + &
                    
                    PI_Myopic(age+1,s,9)*dummy(1)*(rsk(1)-cut(1,ed)) + &
                    PI_Myopic(age+1,s,10)*dummy(2)*(rsk(2)-cut(2,ed)) + &
                    PI_Myopic(age+1,s,11)*dummy(3)*(rsk(3)-cut(3,ed)) + &
                    PI_Myopic(age+1,s,12)*dummy(4)*(rsk(4)-cut(4,ed)) + &
                    PI_Myopic(age+1,s,13)*dummy(5)*(rsk(5)-cut(5,ed)) + &
                    PI_Myopic(age+1,s,14)*dummy(6)*(rsk(6)-cut(6,ed)) + &
                    PI_Myopic(age+1,s,15)*dummy(7)*(rsk(7)-cut(7,ed)) + &
                    
                    PI_Myopic(age+1,s,16)*ExperTomorrow(1) + &
                    PI_Myopic(age+1,s,17)*ExperTomorrow(2) + &
                    PI_Myopic(age+1,s,18)*ExperTomorrow(3) + &
                    PI_Myopic(age+1,s,19)*ExperTomorrow(4) + &
                    PI_Myopic(age+1,s,20)*ExperTomorrow(5) + &
                    PI_Myopic(age+1,s,21)*ExperTomorrow(6) + &
                    PI_Myopic(age+1,s,22)*ExperTomorrow(7) + &
                    
                    PI_Myopic(age+1,s,23)*rsk(1)**2 + &
                    PI_Myopic(age+1,s,24)*rsk(2)**2 + &
                    PI_Myopic(age+1,s,25)*rsk(3)**2 + &
                    PI_Myopic(age+1,s,26)*rsk(4)**2 + &
                    PI_Myopic(age+1,s,27)*rsk(5)**2 + &
                    PI_Myopic(age+1,s,28)*rsk(6)**2 + &
                    PI_Myopic(age+1,s,29)*rsk(7)**2 + &
                    
                    PI_Myopic(age+1,s,30)*dummy(1)*(rsk(1)-cut(1,ed))**2 + &
                    PI_Myopic(age+1,s,31)*dummy(2)*(rsk(2)-cut(2,ed))**2 + &
                    PI_Myopic(age+1,s,32)*dummy(3)*(rsk(3)-cut(3,ed))**2 + &
                    PI_Myopic(age+1,s,33)*dummy(4)*(rsk(4)-cut(4,ed))**2 + &
                    PI_Myopic(age+1,s,34)*dummy(5)*(rsk(5)-cut(5,ed))**2 + &
                    PI_Myopic(age+1,s,35)*dummy(6)*(rsk(6)-cut(6,ed))**2 + &
                    PI_Myopic(age+1,s,36)*dummy(7)*(rsk(7)-cut(7,ed))**2 + &
                
                    PI_Myopic(age+1,s,37)*ExperTomorrow(1)**2 + &
                    PI_Myopic(age+1,s,38)*ExperTomorrow(2)**2 + &
                    PI_Myopic(age+1,s,39)*ExperTomorrow(3)**2 + &
                    PI_Myopic(age+1,s,40)*ExperTomorrow(4)**2 + &
                    PI_Myopic(age+1,s,41)*ExperTomorrow(5)**2 + &
                    PI_Myopic(age+1,s,42)*ExperTomorrow(6)**2 + &
                    PI_Myopic(age+1,s,43)*ExperTomorrow(7)**2 + &
                    
                    PI_Myopic(age+1,s,44)*rsk(1)*rsk(2) + &
                    PI_Myopic(age+1,s,45)*rsk(1)*rsk(3) + &
                    PI_Myopic(age+1,s,46)*rsk(1)*rsk(4) + &
                    PI_Myopic(age+1,s,47)*rsk(1)*rsk(5) + &
                    PI_Myopic(age+1,s,48)*rsk(1)*rsk(6) + &
                    PI_Myopic(age+1,s,49)*rsk(1)*rsk(7) + &
                    
                    PI_Myopic(age+1,s,50)*rsk(2)*rsk(3) + &
                    PI_Myopic(age+1,s,51)*rsk(2)*rsk(4) + &
                    PI_Myopic(age+1,s,52)*rsk(2)*rsk(5) + &
                    PI_Myopic(age+1,s,53)*rsk(2)*rsk(6) + &
                    PI_Myopic(age+1,s,54)*rsk(2)*rsk(7) + &
                    
                    PI_Myopic(age+1,s,55)*rsk(3)*rsk(4) + &
                    PI_Myopic(age+1,s,56)*rsk(3)*rsk(5) + &
                    PI_Myopic(age+1,s,57)*rsk(3)*rsk(6) + &
                    PI_Myopic(age+1,s,58)*rsk(3)*rsk(7) + &
                    
                    PI_Myopic(age+1,s,59)*rsk(4)*rsk(5) + &
                    PI_Myopic(age+1,s,60)*rsk(4)*rsk(6) + &
                    PI_Myopic(age+1,s,61)*rsk(4)*rsk(7) + &
                    
                    PI_Myopic(age+1,s,62)*rsk(5)*rsk(6) + &
                    PI_Myopic(age+1,s,63)*rsk(5)*rsk(7) + &
                    
                    PI_Myopic(age+1,s,64)*rsk(6)*rsk(7) + &
                    
                    PI_Myopic(age+1,s,65)*rsk(1)*ExperTomorrow(1) + &
                    PI_Myopic(age+1,s,66)*rsk(1)*ExperTomorrow(2) + &
                    PI_Myopic(age+1,s,67)*rsk(1)*ExperTomorrow(3) + &
                    PI_Myopic(age+1,s,68)*rsk(1)*ExperTomorrow(4) + &
                    PI_Myopic(age+1,s,69)*rsk(1)*ExperTomorrow(5) + &
                    PI_Myopic(age+1,s,70)*rsk(1)*ExperTomorrow(6) + &
                    PI_Myopic(age+1,s,71)*rsk(1)*ExperTomorrow(7) + &
                    
                    PI_Myopic(age+1,s,72)*rsk(2)*ExperTomorrow(1) + &
                    PI_Myopic(age+1,s,73)*rsk(2)*ExperTomorrow(2) + &
                    PI_Myopic(age+1,s,74)*rsk(2)*ExperTomorrow(3) + &
                    PI_Myopic(age+1,s,75)*rsk(2)*ExperTomorrow(4) + &
                    PI_Myopic(age+1,s,76)*rsk(2)*ExperTomorrow(5) + &
                    PI_Myopic(age+1,s,77)*rsk(2)*ExperTomorrow(6) + &
                    PI_Myopic(age+1,s,78)*rsk(2)*ExperTomorrow(7) + &
                    
                    PI_Myopic(age+1,s,79)*rsk(3)*ExperTomorrow(1) + &
                    PI_Myopic(age+1,s,80)*rsk(3)*ExperTomorrow(2) + &
                    PI_Myopic(age+1,s,81)*rsk(3)*ExperTomorrow(3) + &
                    PI_Myopic(age+1,s,82)*rsk(3)*ExperTomorrow(4) + &
                    PI_Myopic(age+1,s,83)*rsk(3)*ExperTomorrow(5) + &
                    PI_Myopic(age+1,s,84)*rsk(3)*ExperTomorrow(6) + &
                    PI_Myopic(age+1,s,85)*rsk(3)*ExperTomorrow(7) + &
                    
                    PI_Myopic(age+1,s,86)*rsk(4)*ExperTomorrow(1) + &
                    PI_Myopic(age+1,s,87)*rsk(4)*ExperTomorrow(2) + &
                    PI_Myopic(age+1,s,88)*rsk(4)*ExperTomorrow(3) + &
                    PI_Myopic(age+1,s,89)*rsk(4)*ExperTomorrow(4) + &
                    PI_Myopic(age+1,s,90)*rsk(4)*ExperTomorrow(5) + &
                    PI_Myopic(age+1,s,91)*rsk(4)*ExperTomorrow(6) + &
                    PI_Myopic(age+1,s,92)*rsk(4)*ExperTomorrow(7) + &
                    
                    PI_Myopic(age+1,s,93)*rsk(5)*ExperTomorrow(1) + &
                    PI_Myopic(age+1,s,94)*rsk(5)*ExperTomorrow(2) + &
                    PI_Myopic(age+1,s,95)*rsk(5)*ExperTomorrow(3) + &
                    PI_Myopic(age+1,s,96)*rsk(5)*ExperTomorrow(4) + &
                    PI_Myopic(age+1,s,97)*rsk(5)*ExperTomorrow(5) + &
                    PI_Myopic(age+1,s,98)*rsk(5)*ExperTomorrow(6) + &
                    PI_Myopic(age+1,s,99)*rsk(5)*ExperTomorrow(7) + &
                    
                    PI_Myopic(age+1,s,100)*rsk(6)*ExperTomorrow(1) + &
                    PI_Myopic(age+1,s,101)*rsk(6)*ExperTomorrow(2) + &
                    PI_Myopic(age+1,s,102)*rsk(6)*ExperTomorrow(3) + &
                    PI_Myopic(age+1,s,103)*rsk(6)*ExperTomorrow(4) + &
                    PI_Myopic(age+1,s,104)*rsk(6)*ExperTomorrow(5) + &
                    PI_Myopic(age+1,s,105)*rsk(6)*ExperTomorrow(6) + &
                    PI_Myopic(age+1,s,106)*rsk(6)*ExperTomorrow(7) + &
                    
                    PI_Myopic(age+1,s,107)*rsk(7)*ExperTomorrow(1) + &
                    PI_Myopic(age+1,s,108)*rsk(7)*ExperTomorrow(2) + &
                    PI_Myopic(age+1,s,109)*rsk(7)*ExperTomorrow(3) + &
                    PI_Myopic(age+1,s,110)*rsk(7)*ExperTomorrow(4) + &
                    PI_Myopic(age+1,s,111)*rsk(7)*ExperTomorrow(5) + &
                    PI_Myopic(age+1,s,112)*rsk(7)*ExperTomorrow(6) + &
                    PI_Myopic(age+1,s,113)*rsk(7)*ExperTomorrow(7) + &
                    
                    PI_Myopic(age+1,s,114)*ExperTomorrow(1)*ExperTomorrow(2) + &
                    PI_Myopic(age+1,s,115)*ExperTomorrow(1)*ExperTomorrow(3) + &
                    PI_Myopic(age+1,s,116)*ExperTomorrow(1)*ExperTomorrow(4) + &
                    PI_Myopic(age+1,s,117)*ExperTomorrow(1)*ExperTomorrow(5) + &
                    PI_Myopic(age+1,s,118)*ExperTomorrow(1)*ExperTomorrow(6) + &
                    PI_Myopic(age+1,s,119)*ExperTomorrow(1)*ExperTomorrow(7) + &
                    PI_Myopic(age+1,s,120)*ExperTomorrow(2)*ExperTomorrow(3) + &
                    PI_Myopic(age+1,s,121)*ExperTomorrow(2)*ExperTomorrow(4) + &
                    PI_Myopic(age+1,s,122)*ExperTomorrow(2)*ExperTomorrow(5) + &
                    PI_Myopic(age+1,s,123)*ExperTomorrow(2)*ExperTomorrow(6) + &
                    PI_Myopic(age+1,s,124)*ExperTomorrow(2)*ExperTomorrow(7) + &
                    PI_Myopic(age+1,s,125)*ExperTomorrow(3)*ExperTomorrow(4) + &
                    PI_Myopic(age+1,s,126)*ExperTomorrow(3)*ExperTomorrow(5) + &
                    PI_Myopic(age+1,s,127)*ExperTomorrow(3)*ExperTomorrow(6) + &
                    PI_Myopic(age+1,s,128)*ExperTomorrow(3)*ExperTomorrow(7) + &
                    PI_Myopic(age+1,s,129)*ExperTomorrow(4)*ExperTomorrow(5) + &
                    PI_Myopic(age+1,s,130)*ExperTomorrow(4)*ExperTomorrow(6) + &
                    PI_Myopic(age+1,s,131)*ExperTomorrow(4)*ExperTomorrow(7) + &
                    PI_Myopic(age+1,s,132)*ExperTomorrow(5)*ExperTomorrow(6) + &
                    PI_Myopic(age+1,s,133)*ExperTomorrow(5)*ExperTomorrow(7) + &
                    PI_Myopic(age+1,s,134)*ExperTomorrow(6)*ExperTomorrow(7)
                
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
        
        !do kk = 1, INTP
        !    rsk   = X(kk,2:5)
        !    exper = X(kk,10:13)
        !    write(1029,364) gender, educ, tp, age, lag, kk, Y(kk,1), Emax_hat(PI_Myopic(age,lag,:), rsk , exper, cut(:,ed))
        !end do
        
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

    end subroutine





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


    
subroutine EmaxCoef_RE(param,COEF_EXP,covar_EXP,cut,gender,educ,tp,PI_RE,R_sq_RE,FLAG_EMAX)

implicit none

real(KIND=DOUBLE), intent(in)   :: param(:), cut(:,0:), COEF_EXP(:,0:,:,:), covar_EXP(:,:,0:)
integer, intent(in)             :: gender, educ, tp
real(KIND=DOUBLE), intent(out)  :: PI_RE(FIRST_AGE+1:,0:,:), &
                                   R_sq_RE(FIRST_AGE+1:,0:)
integer, intent(out)            :: FLAG_EMAX

real(KIND=DOUBLE) PI_COEF(NREG), CHOL_LOWER(NSECTORS,NSECTORS), CHOL(NSECTORS,NSECTORS)

! ******************
! Counting Variables
! ******************

integer   i, ii, j, t, k, kk, l, n, s, s1, s2, year


! ***********************
! Parameters of the model
! ***********************

real(KIND=DOUBLE) theta(7), beta(NSECTORS,12), & 
                  sigma(0:NSECTORS), kappa(26), tau(0:NSECTORS), SigmaPref, &
                  omega(2:NTYPES,0:NSECTORS), lambda(NTYPES), &
                  alpha_EXP(NSECTORS), beta_EXP(NSECTORS), sigma2_EXP(NSECTORS)


! ***********************
! Wage variables and Emax
! ***********************

real(KIND=DOUBLE) w(0:NSECTORS), hc(NSECTORS), V(0:NSECTORS), VMAX, VM, & 
                  Emax(0:NSECTORS), Emax_init(0:NSECTORS), dummy(NSECTORS), &
                  cost(0:NSECTORS,1:NSECTORS), &
                  SumExpV, Emax0(0:NSECTORS)
                  
integer           exper(NSECTORS), exper8(NSECTORS), lag, age, &
                  educ2, educ3, educ4, ed, ExperTomorrow(NSECTORS)

real(KIND=DOUBLE) term1(NSECTORS), term2(NSECTORS)


! *******************************     
! Variables for the approximation
! *******************************
                 
real(KIND=DOUBLE) COEF(NREG,1), SSE, SST, MeanEmax
real(KIND=DOUBLE) rvec_eps(0:NSECTORS), rvec_rsk(NSECTORS), rvec_xsi(NSECTORS)
real(KIND=DOUBLE), allocatable, dimension(:,:,:) :: eps, xsi
real(KIND=DOUBLE) rsk_draw(NSECTORS), rsk_lag(NSECTORS), rsk_lag2(NSECTORS), new_rsk(NSECTORS), new_rsk_tom(NSECTORS), &
                  log_rsk_draw(NSECTORS), log_new_rsk(NSECTORS), log_new_rsk_tom(NSECTORS), &
                  rsk1(INTP,NSECTORS), rsk2(INTP,NSECTORS)
integer           ch(9)

real(KIND=DOUBLE), allocatable, dimension(:,:)   :: Y
real(KIND=DOUBLE), allocatable, dimension(:,:)   :: X


! ************************
! Random Number Generation
! ************************

TYPE (VSL_STREAM_STATE) :: stream1, stream2, stream3, stream4, stream5
integer(KIND=4) errcode1, errcode2, errcode3, errcode4, errcode5
integer brng1, brng2, brng3, brng4, brng5, &
        method1, method2, method3, method4, method5, &
        seed1, seed2, seed3, seed4, seed5
integer info 


! *********************************************************************************
! *********************************************************************************


allocate(eps(FIRST_AGE+1:LAST_AGE,NDRAWS,0:NSECTORS))
allocate(xsi(FIRST_AGE+1:LAST_AGE,NDRAWS,NSECTORS))
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


! Skilled / Unskilled Worker
if (educ <= 2) then
    ed = 0
else if (educ >= 3) then
    ed = 1
end if

! Parameters of Expectations
alpha_EXP = COEF_EXP(:,ed,1,1)
beta_EXP  = COEF_EXP(:,ed,2,1)

CHOL = covar_EXP(:,:,ed)
call dpotrf('L', NSECTORS, CHOL, NSECTORS, info)
do s1 = 1, NSECTORS
    do s2 = 1, s1
        CHOL_LOWER(s1,s2) = CHOL(s1,s2)
    end do
end do    

do s1 = 1, NSECTORS
    do s2 = s1+1, NSECTORS
        CHOL_LOWER(s1,s2) = 0.0
    end do
end do

do s = 1, NSECTORS
    sigma2_EXP(s) = covar_EXP(s,s,ed)
end do    

! Myopic Expectations
!alpha_EXP = 0.0
!beta_EXP = 0.0
!sigma2_EXP = 0.0

PI_RE   = 0.0
R_sq_RE = 0.0


! *******************************************************************************
! Generating Random Draws for the Monte Carlo Integration
! We generate NDRAWS for each age, in order to save computational time
! Just as a reminder, in the implementation of Keane and Wolpin (1994),
! Drawing NDRAWS just once and for all the states and all ages did not perform
! very well for approximating Emax. Drawing the same points for all states by age
! performed much better
! *******************************************************************************

! *********************************
! Random Generator Initialization 
! for HC Shocks and HC price shocks
! *********************************

! epsilon shocks (HC)
brng1    = VSL_BRNG_MT19937
method1  = VSL_METHOD_DGAUSSIAN_BOXMULLER
seed1    = 12345678
errcode1 = vslnewstream(stream1,brng1,seed1)

! xsi shocks (HC prices)
brng2     = VSL_BRNG_MT19937
method2   = VSL_METHOD_DGAUSSIAN_BOXMULLER
seed2     = 98767241
errcode2  = vslnewstream(stream2,brng2,seed2)

! ***** Warming Up *****
do j = 1, WARMUP
errcode1  = vdrnggaussian(method1, stream1, NSECTORS+1, &
            rvec_eps, real(0.0,DOUBLE), real(1.0,DOUBLE))
errcode2  = vdrnggaussian(method2, stream2, NSECTORS, &
            rvec_xsi, real(0.0,DOUBLE), real(1.0,DOUBLE))
end do                        

do age = FIRST_AGE+1, LAST_AGE
    do j = 1, NDRAWS
        errcode1  = vdrnggaussian(method1, stream1, NSECTORS+1, &
                    rvec_eps, real(0.0,DOUBLE), real(1.0,DOUBLE))
        errcode2  = vdrnggaussian(method2, stream2, NSECTORS, &
                    rvec_xsi, real(0.0,DOUBLE), real(1.0,DOUBLE))
        do i = 0, NSECTORS
            eps(age,j,i) = sigma(i)*rvec_eps(i)
        end do
        do s1 = 1, NSECTORS
            xsi(age,j,s1) = 0.0
            do s2 = 1, s1
                xsi(age,j,s1) = xsi(age,j,s1) + CHOL_LOWER(s1,s2)*rvec_xsi(s2)
            end do
        end do
    end do
end do


    
! **********************************************
! Deinitialize Stream for HC and HC price Shocks
! **********************************************

errcode1 = vsldeletestream(stream1)
errcode2 = vsldeletestream(stream2)



! ********************************
! Random Generator Initialization 
! for Returns to Skill And Choices
! ********************************

! ***** Initializing *****
brng3    = VSL_BRNG_MT19937
seed3    = 98765456
method3  = VSL_METHOD_IUNIFORM_STD
errcode3 = vslnewstream(stream3, brng3, seed3)

brng4    = VSL_BRNG_MT19937
seed4    = 65745423
method4  = VSL_METHOD_IUNIFORM_STD
errcode4 = vslnewstream(stream4, brng4, seed4)

brng5    = VSL_BRNG_MT19937
seed5    = 25109825
method5  = VSL_METHOD_IBERNOULLI_ICDF
errcode5 = vslnewstream(stream5, brng5, seed5)




! ***** Warming Up *****
do j = 1, WARMUP
    errcode3 = vdrnguniform(method3, stream3, NSECTORS, rsk_lag, &
               real(0.0,DOUBLE), real(1.0,DOUBLE))
    errcode4 = vdrnguniform(method4, stream4, NSECTORS, rsk_lag2, &
               real(0.0,DOUBLE), real(1.0,DOUBLE))
    errcode5 = virnguniform(method5, stream5, 9, ch, 0, NSECTORS+1)    
end do 

    
age_loop: do age = LAST_AGE, (FIRST_AGE+1), -1
    
    if (FLAG_EMAX == 1) then
        exit
    end if
    
    lag_loop: do lag = 0, NSECTORS
        
        if (FLAG_EMAX == 1) then
            exit
        end if

        INTP_loop: do kk = 1, INTP
        
    
            ! *****************************************************************
            ! Random Number Generation - returns to skill and past choices
            ! *****************************************************************
           
            errcode3 = vdrnguniform(method3, stream3, NSECTORS, rvec_rsk, &
                       real(0.0,DOUBLE), real(1.0,DOUBLE))
        
            rsk_lag = LOW_RSK(:,ed) + (UP_RSK(:,ed) - LOW_RSK(:,ed))*rvec_rsk
                
            errcode4 = vdrnguniform(method4, stream4, NSECTORS, rvec_rsk, &
                       real(0.0,DOUBLE), real(1.0,DOUBLE))
                
            rsk_lag2 = 0.85*rsk_lag + 0.3*rsk_lag*rvec_rsk 
            
            ! rsk_lag2 = LOW_RSK(:,ed) + (UP_RSK(:,ed) - LOW_RSK(:,ed))*rvec_rsk
          
            errcode5 = virnguniform(method5, stream5, 8, ch(2:9), 0, NSECTORS+1) 
        
            ! *****************************************************************
            ! *****************************************************************
                       
            ch(1) = lag
            
            ! experience accumulated in the last 9 years
            exper = 0
            do s = 1, NSECTORS
                do i = 1, 9
                    if (ch(i) == s) then
                        exper(s) = exper(s) + 1
                    end if
                end do
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
            
            ! human capital net of shocks
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
            
            ! experience accumulated in the last 8 years
            exper8 = 0
            do s = 1, NSECTORS
                do i = 1, 8
                    if (ch(i) == s) then
                        exper8(s) = exper8(s) + 1
                    end if
                end do
            end do
            
            
            if (age < LAST_AGE) then    
                
                do s = 0, NSECTORS
                    
                    ExperTomorrow = exper8
                    if (s >= 1) then
                        ExperTomorrow(s) = exper8(s) + 1
                    end if
                            
                    Emax_init(s) = &
                        
                    PI_RE(age+1,s,1)*1.0 + &
                    
                    PI_RE(age+1,s,16)*ExperTomorrow(1) + &
                    PI_RE(age+1,s,17)*ExperTomorrow(2) + &
                    PI_RE(age+1,s,18)*ExperTomorrow(3) + &
                    PI_RE(age+1,s,19)*ExperTomorrow(4) + &
                    PI_RE(age+1,s,20)*ExperTomorrow(5) + &
                    PI_RE(age+1,s,21)*ExperTomorrow(6) + &
                    PI_RE(age+1,s,22)*ExperTomorrow(7) + &
                    
                    PI_RE(age+1,s,37)*ExperTomorrow(1)**2 + &
                    PI_RE(age+1,s,38)*ExperTomorrow(2)**2 + &
                    PI_RE(age+1,s,39)*ExperTomorrow(3)**2 + &
                    PI_RE(age+1,s,40)*ExperTomorrow(4)**2 + &
                    PI_RE(age+1,s,41)*ExperTomorrow(5)**2 + &
                    PI_RE(age+1,s,42)*ExperTomorrow(6)**2 + &
                    PI_RE(age+1,s,43)*ExperTomorrow(7)**2 + &
                    
                    PI_RE(age+1,s,114)*ExperTomorrow(1)*ExperTomorrow(2) + &
                    PI_RE(age+1,s,115)*ExperTomorrow(1)*ExperTomorrow(3) + &
                    PI_RE(age+1,s,116)*ExperTomorrow(1)*ExperTomorrow(4) + &
                    PI_RE(age+1,s,117)*ExperTomorrow(1)*ExperTomorrow(5) + &
                    PI_RE(age+1,s,118)*ExperTomorrow(1)*ExperTomorrow(6) + &
                    PI_RE(age+1,s,119)*ExperTomorrow(1)*ExperTomorrow(7) + &
                    PI_RE(age+1,s,120)*ExperTomorrow(2)*ExperTomorrow(3) + &
                    PI_RE(age+1,s,121)*ExperTomorrow(2)*ExperTomorrow(4) + &
                    PI_RE(age+1,s,122)*ExperTomorrow(2)*ExperTomorrow(5) + &
                    PI_RE(age+1,s,123)*ExperTomorrow(2)*ExperTomorrow(6) + &
                    PI_RE(age+1,s,124)*ExperTomorrow(2)*ExperTomorrow(7) + &
                    PI_RE(age+1,s,125)*ExperTomorrow(3)*ExperTomorrow(4) + &
                    PI_RE(age+1,s,126)*ExperTomorrow(3)*ExperTomorrow(5) + &
                    PI_RE(age+1,s,127)*ExperTomorrow(3)*ExperTomorrow(6) + &
                    PI_RE(age+1,s,128)*ExperTomorrow(3)*ExperTomorrow(7) + &
                    PI_RE(age+1,s,129)*ExperTomorrow(4)*ExperTomorrow(5) + &
                    PI_RE(age+1,s,130)*ExperTomorrow(4)*ExperTomorrow(6) + &
                    PI_RE(age+1,s,131)*ExperTomorrow(4)*ExperTomorrow(7) + &
                    PI_RE(age+1,s,132)*ExperTomorrow(5)*ExperTomorrow(6) + &
                    PI_RE(age+1,s,133)*ExperTomorrow(5)*ExperTomorrow(7) + &
                    PI_RE(age+1,s,134)*ExperTomorrow(6)*ExperTomorrow(7)
                        
                end do
                                    
                MeanEmax = 0.0
                do n = 1, NDRAWS
                    
                    ! contemporaneous rsk draw, conditional on rsk_lag and rsk_lag2
                    do s = 1, NSECTORS
                        rsk_draw(s) = rsk_lag(s) * ((rsk_lag(s)/rsk_lag2(s))**beta_EXP(s)) * exp(alpha_EXP(s) + xsi(age,n,s))
                    end do
                    
                    do s = 1 , NSECTORS
                        new_rsk_tom(s) = rsk_draw(s) * ((rsk_draw(s)/rsk_lag(s))**beta_EXP(s)) * exp(alpha_EXP(s) + sigma2_EXP(s)/2)
                    end do
                        
                    dummy = 0.0
                    do s = 1, NSECTORS
                        if (new_rsk_tom(s) >= cut(s,ed)) then
                            dummy(s) = 1
                        end if
                    end do
                    
                    do s = 0, NSECTORS
                    
                        ExperTomorrow = exper8
                        if (s >= 1) then
                            ExperTomorrow(s) = exper8(s) + 1
                        end if
                            
                        Emax(s) = Emax_init(s) + &
                        
                        !PI_RE(age+1,s,1)*1.0 + &
                        
                        PI_RE(age+1,s,2)*new_rsk_tom(1) + &
                        PI_RE(age+1,s,3)*new_rsk_tom(2) + &
                        PI_RE(age+1,s,4)*new_rsk_tom(3) + &
                        PI_RE(age+1,s,5)*new_rsk_tom(4) + &
                        PI_RE(age+1,s,6)*new_rsk_tom(5) + &
                        PI_RE(age+1,s,7)*new_rsk_tom(6) + &
                        PI_RE(age+1,s,8)*new_rsk_tom(7) + &
                    
                        PI_RE(age+1,s,9)*dummy(1)*(new_rsk_tom(1)-cut(1,ed)) + &
                        PI_RE(age+1,s,10)*dummy(2)*(new_rsk_tom(2)-cut(2,ed)) + &
                        PI_RE(age+1,s,11)*dummy(3)*(new_rsk_tom(3)-cut(3,ed)) + &
                        PI_RE(age+1,s,12)*dummy(4)*(new_rsk_tom(4)-cut(4,ed)) + &
                        PI_RE(age+1,s,13)*dummy(5)*(new_rsk_tom(5)-cut(5,ed)) + &
                        PI_RE(age+1,s,14)*dummy(6)*(new_rsk_tom(6)-cut(6,ed)) + &
                        PI_RE(age+1,s,15)*dummy(7)*(new_rsk_tom(7)-cut(7,ed)) + &
                    
                        !PI_RE(age+1,s,16)*ExperTomorrow(1) + &
                        !PI_RE(age+1,s,17)*ExperTomorrow(2) + &
                        !PI_RE(age+1,s,18)*ExperTomorrow(3) + &
                        !PI_RE(age+1,s,19)*ExperTomorrow(4) + &
                        !PI_RE(age+1,s,20)*ExperTomorrow(5) + &
                        !PI_RE(age+1,s,21)*ExperTomorrow(6) + &
                        !PI_RE(age+1,s,22)*ExperTomorrow(7) + &
                    
                        PI_RE(age+1,s,23)*new_rsk_tom(1)**2 + &
                        PI_RE(age+1,s,24)*new_rsk_tom(2)**2 + &
                        PI_RE(age+1,s,25)*new_rsk_tom(3)**2 + &
                        PI_RE(age+1,s,26)*new_rsk_tom(4)**2 + &
                        PI_RE(age+1,s,27)*new_rsk_tom(5)**2 + &
                        PI_RE(age+1,s,28)*new_rsk_tom(6)**2 + &
                        PI_RE(age+1,s,29)*new_rsk_tom(7)**2 + &
                    
                        PI_RE(age+1,s,30)*dummy(1)*(new_rsk_tom(1)-cut(1,ed))**2 + &
                        PI_RE(age+1,s,31)*dummy(2)*(new_rsk_tom(2)-cut(2,ed))**2 + &
                        PI_RE(age+1,s,32)*dummy(3)*(new_rsk_tom(3)-cut(3,ed))**2 + &
                        PI_RE(age+1,s,33)*dummy(4)*(new_rsk_tom(4)-cut(4,ed))**2 + &
                        PI_RE(age+1,s,34)*dummy(5)*(new_rsk_tom(5)-cut(5,ed))**2 + &
                        PI_RE(age+1,s,35)*dummy(6)*(new_rsk_tom(6)-cut(6,ed))**2 + &
                        PI_RE(age+1,s,36)*dummy(7)*(new_rsk_tom(7)-cut(7,ed))**2 + &
                
                        !PI_RE(age+1,s,37)*ExperTomorrow(1)**2 + &
                        !PI_RE(age+1,s,38)*ExperTomorrow(2)**2 + &
                        !PI_RE(age+1,s,39)*ExperTomorrow(3)**2 + &
                        !PI_RE(age+1,s,40)*ExperTomorrow(4)**2 + &
                        !PI_RE(age+1,s,41)*ExperTomorrow(5)**2 + &
                        !PI_RE(age+1,s,42)*ExperTomorrow(6)**2 + &
                        !PI_RE(age+1,s,43)*ExperTomorrow(7)**2 + &
                    
                        PI_RE(age+1,s,44)*new_rsk_tom(1)*new_rsk_tom(2) + &
                        PI_RE(age+1,s,45)*new_rsk_tom(1)*new_rsk_tom(3) + &
                        PI_RE(age+1,s,46)*new_rsk_tom(1)*new_rsk_tom(4) + &
                        PI_RE(age+1,s,47)*new_rsk_tom(1)*new_rsk_tom(5) + &
                        PI_RE(age+1,s,48)*new_rsk_tom(1)*new_rsk_tom(6) + &
                        PI_RE(age+1,s,49)*new_rsk_tom(1)*new_rsk_tom(7) + &
                    
                        PI_RE(age+1,s,50)*new_rsk_tom(2)*new_rsk_tom(3) + &
                        PI_RE(age+1,s,51)*new_rsk_tom(2)*new_rsk_tom(4) + &
                        PI_RE(age+1,s,52)*new_rsk_tom(2)*new_rsk_tom(5) + &
                        PI_RE(age+1,s,53)*new_rsk_tom(2)*new_rsk_tom(6) + &
                        PI_RE(age+1,s,54)*new_rsk_tom(2)*new_rsk_tom(7) + &
                    
                        PI_RE(age+1,s,55)*new_rsk_tom(3)*new_rsk_tom(4) + &
                        PI_RE(age+1,s,56)*new_rsk_tom(3)*new_rsk_tom(5) + &
                        PI_RE(age+1,s,57)*new_rsk_tom(3)*new_rsk_tom(6) + &
                        PI_RE(age+1,s,58)*new_rsk_tom(3)*new_rsk_tom(7) + &
                    
                        PI_RE(age+1,s,59)*new_rsk_tom(4)*new_rsk_tom(5) + &
                        PI_RE(age+1,s,60)*new_rsk_tom(4)*new_rsk_tom(6) + &
                        PI_RE(age+1,s,61)*new_rsk_tom(4)*new_rsk_tom(7) + &
                    
                        PI_RE(age+1,s,62)*new_rsk_tom(5)*new_rsk_tom(6) + &
                        PI_RE(age+1,s,63)*new_rsk_tom(5)*new_rsk_tom(7) + &
                    
                        PI_RE(age+1,s,64)*new_rsk_tom(6)*new_rsk_tom(7) + &
                    
                        PI_RE(age+1,s,65)*new_rsk_tom(1)*ExperTomorrow(1) + &
                        PI_RE(age+1,s,66)*new_rsk_tom(1)*ExperTomorrow(2) + &
                        PI_RE(age+1,s,67)*new_rsk_tom(1)*ExperTomorrow(3) + &
                        PI_RE(age+1,s,68)*new_rsk_tom(1)*ExperTomorrow(4) + &
                        PI_RE(age+1,s,69)*new_rsk_tom(1)*ExperTomorrow(5) + &
                        PI_RE(age+1,s,70)*new_rsk_tom(1)*ExperTomorrow(6) + &
                        PI_RE(age+1,s,71)*new_rsk_tom(1)*ExperTomorrow(7) + &
                    
                        PI_RE(age+1,s,72)*new_rsk_tom(2)*ExperTomorrow(1) + &
                        PI_RE(age+1,s,73)*new_rsk_tom(2)*ExperTomorrow(2) + &
                        PI_RE(age+1,s,74)*new_rsk_tom(2)*ExperTomorrow(3) + &
                        PI_RE(age+1,s,75)*new_rsk_tom(2)*ExperTomorrow(4) + &
                        PI_RE(age+1,s,76)*new_rsk_tom(2)*ExperTomorrow(5) + &
                        PI_RE(age+1,s,77)*new_rsk_tom(2)*ExperTomorrow(6) + &
                        PI_RE(age+1,s,78)*new_rsk_tom(2)*ExperTomorrow(7) + &
                    
                        PI_RE(age+1,s,79)*new_rsk_tom(3)*ExperTomorrow(1) + &
                        PI_RE(age+1,s,80)*new_rsk_tom(3)*ExperTomorrow(2) + &
                        PI_RE(age+1,s,81)*new_rsk_tom(3)*ExperTomorrow(3) + &
                        PI_RE(age+1,s,82)*new_rsk_tom(3)*ExperTomorrow(4) + &
                        PI_RE(age+1,s,83)*new_rsk_tom(3)*ExperTomorrow(5) + &
                        PI_RE(age+1,s,84)*new_rsk_tom(3)*ExperTomorrow(6) + &
                        PI_RE(age+1,s,85)*new_rsk_tom(3)*ExperTomorrow(7) + &
                    
                        PI_RE(age+1,s,86)*new_rsk_tom(4)*ExperTomorrow(1) + &
                        PI_RE(age+1,s,87)*new_rsk_tom(4)*ExperTomorrow(2) + &
                        PI_RE(age+1,s,88)*new_rsk_tom(4)*ExperTomorrow(3) + &
                        PI_RE(age+1,s,89)*new_rsk_tom(4)*ExperTomorrow(4) + &
                        PI_RE(age+1,s,90)*new_rsk_tom(4)*ExperTomorrow(5) + &
                        PI_RE(age+1,s,91)*new_rsk_tom(4)*ExperTomorrow(6) + &
                        PI_RE(age+1,s,92)*new_rsk_tom(4)*ExperTomorrow(7) + &
                    
                        PI_RE(age+1,s,93)*new_rsk_tom(5)*ExperTomorrow(1) + &
                        PI_RE(age+1,s,94)*new_rsk_tom(5)*ExperTomorrow(2) + &
                        PI_RE(age+1,s,95)*new_rsk_tom(5)*ExperTomorrow(3) + &
                        PI_RE(age+1,s,96)*new_rsk_tom(5)*ExperTomorrow(4) + &
                        PI_RE(age+1,s,97)*new_rsk_tom(5)*ExperTomorrow(5) + &
                        PI_RE(age+1,s,98)*new_rsk_tom(5)*ExperTomorrow(6) + &
                        PI_RE(age+1,s,99)*new_rsk_tom(5)*ExperTomorrow(7) + &
                    
                        PI_RE(age+1,s,100)*new_rsk_tom(6)*ExperTomorrow(1) + &
                        PI_RE(age+1,s,101)*new_rsk_tom(6)*ExperTomorrow(2) + &
                        PI_RE(age+1,s,102)*new_rsk_tom(6)*ExperTomorrow(3) + &
                        PI_RE(age+1,s,103)*new_rsk_tom(6)*ExperTomorrow(4) + &
                        PI_RE(age+1,s,104)*new_rsk_tom(6)*ExperTomorrow(5) + &
                        PI_RE(age+1,s,105)*new_rsk_tom(6)*ExperTomorrow(6) + &
                        PI_RE(age+1,s,106)*new_rsk_tom(6)*ExperTomorrow(7) + &
                    
                        PI_RE(age+1,s,107)*new_rsk_tom(7)*ExperTomorrow(1) + &
                        PI_RE(age+1,s,108)*new_rsk_tom(7)*ExperTomorrow(2) + &
                        PI_RE(age+1,s,109)*new_rsk_tom(7)*ExperTomorrow(3) + &
                        PI_RE(age+1,s,110)*new_rsk_tom(7)*ExperTomorrow(4) + &
                        PI_RE(age+1,s,111)*new_rsk_tom(7)*ExperTomorrow(5) + &
                        PI_RE(age+1,s,112)*new_rsk_tom(7)*ExperTomorrow(6) + &
                        PI_RE(age+1,s,113)*new_rsk_tom(7)*ExperTomorrow(7)
                    
                        !PI_RE(age+1,s,114)*ExperTomorrow(1)*ExperTomorrow(2) + &
                        !PI_RE(age+1,s,115)*ExperTomorrow(1)*ExperTomorrow(3) + &
                        !PI_RE(age+1,s,116)*ExperTomorrow(1)*ExperTomorrow(4) + &
                        !PI_RE(age+1,s,117)*ExperTomorrow(1)*ExperTomorrow(5) + &
                        !PI_RE(age+1,s,118)*ExperTomorrow(1)*ExperTomorrow(6) + &
                        !PI_RE(age+1,s,119)*ExperTomorrow(1)*ExperTomorrow(7) + &
                        !PI_RE(age+1,s,120)*ExperTomorrow(2)*ExperTomorrow(3) + &
                        !PI_RE(age+1,s,121)*ExperTomorrow(2)*ExperTomorrow(4) + &
                        !PI_RE(age+1,s,122)*ExperTomorrow(2)*ExperTomorrow(5) + &
                        !PI_RE(age+1,s,123)*ExperTomorrow(2)*ExperTomorrow(6) + &
                        !PI_RE(age+1,s,124)*ExperTomorrow(2)*ExperTomorrow(7) + &
                        !PI_RE(age+1,s,125)*ExperTomorrow(3)*ExperTomorrow(4) + &
                        !PI_RE(age+1,s,126)*ExperTomorrow(3)*ExperTomorrow(5) + &
                        !PI_RE(age+1,s,127)*ExperTomorrow(3)*ExperTomorrow(6) + &
                        !PI_RE(age+1,s,128)*ExperTomorrow(3)*ExperTomorrow(7) + &
                        !PI_RE(age+1,s,129)*ExperTomorrow(4)*ExperTomorrow(5) + &
                        !PI_RE(age+1,s,130)*ExperTomorrow(4)*ExperTomorrow(6) + &
                        !PI_RE(age+1,s,131)*ExperTomorrow(4)*ExperTomorrow(7) + &
                        !PI_RE(age+1,s,132)*ExperTomorrow(5)*ExperTomorrow(6) + &
                        !PI_RE(age+1,s,133)*ExperTomorrow(5)*ExperTomorrow(7) + &
                        !PI_RE(age+1,s,134)*ExperTomorrow(6)*ExperTomorrow(7)
                                                             
                    end do
                    
                    V(0) = tau(0) + w(0) + eps(age,n,0) + rho*Emax(0)
                    do s = 1, NSECTORS
                        V(s) = tau(s) + rsk_draw(s)*hc(s)*exp(eps(age,n,s)) - cost(lag,s) + rho*Emax(s)
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
                
            else
                    
                MeanEmax = 0.0
                do n = 1, NDRAWS
                    
                    ! contemporaneous rsk draw, conditional on rsk_lag and rsk_lag2
                    do s = 1, NSECTORS
                        rsk_draw(s) = rsk_lag(s) * ((rsk_lag(s)/rsk_lag2(s))**beta_EXP(s)) * exp(alpha_EXP(s) + xsi(age,n,s))
                    end do                        
                        
                    V(0) = tau(0) + w(0) + eps(age,n,0)
                    do s = 1, NSECTORS
                        V(s) = tau(s) + rsk_draw(s)*hc(s)*exp(eps(age,n,s)) - cost(lag,s)
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
            
            end if
                
            ! We approximate Emax with a 2nd order polynomial
            ! In the state variables, in each of the regions
            ! Defined by the vector cut
            
            do s = 1 , NSECTORS
                 new_rsk(s) = rsk_lag(s) * ((rsk_lag(s)/rsk_lag2(s))**beta_EXP(s)) * exp(alpha_EXP(s) + sigma2_EXP(s)/2)
            end do

            dummy = 0.0
            do s = 1, NSECTORS
                if (new_rsk(s) >= cut(s,ed)) then
                    dummy(s) = 1
                end if
            end do

            
            Y(kk,1)       = MeanEmax
        
            X(kk,1)       = 1.0
            
            X(kk,2:8)     = new_rsk(1:7)
            
            X(kk,9:15)    = dummy(1:7)*(new_rsk(1:7)-cut(1:7,ed))
            
            X(kk,16:22)   = exper(1:7)
                        
            X(kk,23:29)   = new_rsk(1:7)**2
            
            X(kk,30:36)   = dummy(1:7)*(new_rsk(1:7)-cut(1:7,ed))**2
            
            X(kk,37:43)   = exper(1:7)**2
            
            X(kk,44:49)   = new_rsk(1)*new_rsk(2:7)
           
            X(kk,50:54)   = new_rsk(2)*new_rsk(3:7)
          
            X(kk,55:58)   = new_rsk(3)*new_rsk(4:7)
            
            X(kk,59:61) = new_rsk(4)*new_rsk(5:7)
            
            X(kk,62:63) = new_rsk(5)*new_rsk(6:7)
            
            X(kk,64)     = new_rsk(6)*new_rsk(7)
            
            X(kk,65:71) = new_rsk(1)*exper(1:7)
                        
            X(kk,72:78) = new_rsk(2)*exper(1:7)
            
            X(kk,79:85) = new_rsk(3)*exper(1:7)
            
            X(kk,86:92) = new_rsk(4)*exper(1:7)
            
            X(kk,93:99) = new_rsk(5)*exper(1:7)
            
            X(kk,100:106) = new_rsk(6)*exper(1:7)
            
            X(kk,107:113) = new_rsk(7)*exper(1:7)
            
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
        
        if (info /= 0) then
            print*, 'Problem running the regression'
            FLAG_EMAX = 1
            exit
        end if            
             
        ! ********************************
        ! ********************************
    
        PI_RE(age,lag,:) = COEF(:,1)
        R_sq_RE(age,lag) = 1 - (SSE/SST)
            
    end do lag_loop

end do age_loop   

! *********************************
! Deinitialize Stream for 
! returns to skill and past choices
! *********************************

errcode3 = vsldeletestream(stream3)
errcode4 = vsldeletestream(stream4)
errcode5 = vsldeletestream(stream5)

deallocate(Y,X,eps,xsi)

364 format(6(i6,','),2(f20.8,','))

end subroutine EmaxCoef_RE    


end module Emax_MOD

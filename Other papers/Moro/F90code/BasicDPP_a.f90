! CODER: Meghan
! DATE OF THIS CODE: February 19, 2015
! CODE NAME: BasicDPP_1.f90

!*****************************************************************************************************
! NOTES: This is a dynamic programming code that simulates a 5 period 
! 	model where agents decide whether to work or not.
!	For those who worked in the prior period, offers are available with certainty.  
!	For those who did not work, an offer occurs with some probability.
!	In the terminal period, we assume no one gets offers, and thus no one works.
!*****************************************************************************************************

module initialize_data
  implicit none

! DOUBLE PRECISION

  integer, parameter :: prec = 8
  real(prec), parameter :: zero=0._prec, one=1._prec

! INDICES

  integer :: a, lw, ed, xper, mar, i_g, t, n, i, d

! DIMENSIONS 

  integer,parameter :: dt		= 5000		! Number of draws for expected values (i.e. for Monte Carlo integration).
  integer,parameter :: nobs		= 1500	! Total number of simulated individuals.
  integer,parameter :: nt 		= 5	! Number of time periods in the dynamic programming problem.
  integer,parameter :: ntdata 		= 5		! Number of time periods to be simulated.
								! Note, in practice, in many cases nt and ntdata won't be the same.
								! For example, you might want to solve the model past the time period you observe in your data (i.e. nt>ntdata).
  integer,parameter :: nsh_u		= 1		! Number of utility shocks.
  integer,parameter :: maxexp	= 4		! Highest experience allowed.  
  integer,parameter :: numst		= 60		! Number of combinations of state space elements in the model. Must change when maxexp changes.
								! Or when additional states or dimensions are added to the stock index.
								! Currently 5x3x2x2=60.  See below for variables in the state space.
  integer,parameter :: sim_data 	= 11		! Unit for data file to be written.
  integer,parameter :: init_data	= 12		! Unit for initial conditions to be read in from.
  integer,parameter :: nvs      		= 5		! Number of variables in dataset of initial conditions.

! In what follows, when a dimension is marked with a ':' that denotes we specify the exact dimension
! in the subroutine that follows called allocate_equations.

! CHOICES

   integer, dimension(:,:), allocatable :: work	! Decisions are to work or not.  0 denotes non-work and 1 denotes work.

! STATE VARIABLES

  integer, dimension(:,:), allocatable :: exper, lastwork, married
				
				! exper is years of experience.
				! lastwork is an indicator for whether the individual worked last period or not.
				! married is an indicator for whether the individual is married or not.

  integer, dimension(:), allocatable :: educ

				! Education categories, which are time-invariant, hence only 1 dimension.

  integer, dimension(:,:,:,:), allocatable :: stock

				! Tracks experience, education categories, whether worked last period or not, whether married or not.
				! We do this because the maximum rank (i.e. number of dimensions) of an array in F90 is 7.
				! We can get around this by putting an array (i.e. stock) within an array.

! ERRORS

! Includes two sets of errors, one to be used in value function calculations indexed by d.
! One to be used in the simulations indexed by n.

	integer, parameter:: LSTATE=16
 	double precision, dimension(:,:), allocatable:: ush, ush_sim		! Utility shock to working. 
												! Note that if you had more than one shock, you would likely add a dimension here.
	double precision, dimension(:,:), allocatable:: wsh,wsh_sim		! Wage shock.

 	integer :: IRANK, ISEED, STATE(LSTATE),INFO
 	double precision :: cov_u(nsh_u,nsh_u), var_w				! Covariance matrix for utility shocks and variance of wage shock.
												! Even though there is only 1 utility shock right now, I code it this way so it 
												! is easy to extend the model to include more choices and shocks.
 double precision :: XMU_U(nsh_u), XMU_W					! Variables for the mean of the utility shock and the wage shock.
 double precision :: shock_u_val(dt,nsh_u), shock_w_val(dt)		! Vector that will hold the draws of utility and wage shocks for calculating value functions.
 double precision :: shock_u_sim(nobs,nsh_u), shock_w_sim(nobs)	! Vector that holds the draws of utility and wage shocks for the simulations.

! COVARIANCE-MATRIX AND VARIANCE TERMS

! Note there is a covariance-variance matrix for utility shocks and a wage shock variance.
! Below are the parameters that govern those variances.
! This will become clearer in the subroutine that follows.

real(prec):: alpha_cov_u_1_1, alpha_var_wsh

! UTILITY, WAGE, NON-LABOR INCOME, VALUE FUNCTIONS

! Capital E in front means expected value.
! sim in front means used in simulations when we simulate individuals' choices,
! and not used in the dynamic programming problem solution.

  real(prec), dimension(:,:,:), allocatable     		:: utfunc
  real(prec), dimension(:,:,:), allocatable			:: sim_utfunc
  real(prec), dimension(:,:), allocatable			:: wage
  real(prec), dimension(:,:), allocatable			:: sim_wage
  real(prec), dimension(:,:,:), allocatable			:: nlincome
  real(prec), dimension(:,:,:), allocatable			:: sim_nlincome
  real(prec), dimension(:,:,:), allocatable			:: v
  real(prec), dimension(:,:), allocatable 			:: v_off_max, sumv_off_max, Ev_off_max
  real(prec), dimension(:,:), allocatable 			:: v_no_max, sumv_no_max, Ev_no_max
  real(prec), dimension(:,:,:), allocatable			:: sim_v

! RANDOM NUMBER TO DETERMINE IF OFFER RECEIVED

  real(prec), dimension(:,:), allocatable			:: joffrn

! UNIFORM DISTRIBUTION RANDOM NUMBER SEED

! Necessary for job offer random numbers and for replicability.

  integer:: sizer
  integer, allocatable :: joffseed(:)	

! MODEL PARAMETERS

real(prec) :: offerprob(0:2,0:1)		! Job offer probability.  First index is education, second is married or not.
							! For those who worked last period they get an offer with certainty this period.
							! For those who did not work last period they get an offer with probability offerprob this period.

real(prec) :: offerprob_intercept		! Intercept in job offer probability.
real(prec) :: offerprob_educ(0:2)		! Coefficient on education in job offer probability.
real(prec) :: offerprob_married		! Coefficient on married in job offer probability.

real(prec) :: w_intercept			! Intercept in earnings equation.
real(prec) :: w_xper, w_xpersq		! Coefficients on experience and experience squared in earnings equation.
real(prec) :: w_educ(0:2)			! Coefficient on education indicators in earnings equation.
real(prec) :: w_lastnw				! Coefficient on not working last period in earnings equation.


real(prec) :: u_work				! Disutility from work.

real(prec) :: nli_intercept			! Intercept in non-labor income equation.
real(prec) :: nli_work				! Coefficient on working in non-labor income equation.
real(prec) :: nli_educ(0:2)			! Coefficient on education indicators in non-labor income equation.
real(prec) :: nli_married			! Coefficient on being married in non-labor income equation.

real(prec) :: beta					! Discount rate.

! DATAFILE OF INITIAL CONDITIONS

  real(prec),dimension(:,:), allocatable :: initdata_file

! FOR CALCULATING MODEL STATISTICS OR OTHER NEEDED VARIABLES

  real(prec), dimension(:), allocatable :: totwork, pctwork
  integer, dimension(:,:), allocatable  :: receivedoffer			! This is an indicator used in the simulation section
											! that tracks whether an individual received a job offer
											! regardless of whether he accepted it or not.		

end module initialize_data

!******************************************************************************
! Here is where all the external subroutines start.
!******************************************************************************

!******************************************************************************
! Allocate space for the variables.
!******************************************************************************

subroutine allocate_equations
  use initialize_data
  implicit none

! CHOICES

allocate (work(1:nobs,1:ntdata))				! The dimensions of the simulated work choice are n and t (the individual and time).

! STOCHASTIC ELEMENTS

allocate (ush(1:dt,1:nt))  					! Utility shock.  First index is for draw, second is for time.  Used in DPP solution.
allocate (ush_sim(1:nobs,1:ntdata))			! Utility shock.  Used in simulations.
allocate (wsh(1:dt,1:nt))					! Shock to wage.  Used in DPP solution.
allocate (wsh_sim(1:nobs,1:ntdata))			! Shock to wage.  Used in simulations.
allocate (joffrn(1:nobs,1:ntdata))				! Random number to determine whether an individual got a job offer.
									! Used in simulations.

! STATE VARIABLES

allocate (exper(1:nobs,1:ntdata))				
allocate (lastwork(1:nobs,1:ntdata))			
allocate (educ(1:nobs))					! Assume education is time-invariant.
allocate (married(1:nobs,1:ntdata))
allocate (stock(0:maxexp,0:2,0:1,0:1))			! Experience (0:maxexp); education category-less than HS, HS, some college (0:2); worked or not last period (0:1); 
									! married or not (0:1).

! UTILITY FUNCTIONS, WAGES, NON-LABOR INCOME, VALUE FUNCTIONS

allocate (utfunc(0:1,1:numst,1:nt))				! Utility-work or not (0:1), stock array (1:numst), time (1:nt).
allocate (wage(1:numst,1:nt))				! Wage-stock array (1:numst), time (1:nt).
allocate (nlincome(0:1,1:numst,1:nt))			! NL Income-work or not (0:1), stock array (1:numst), time (1:nt).
allocate (v(0:1,1:numst,1:nt))					! Value function-work or not (0:1), stock array (1:numst), time (1:nt).
allocate (v_off_max(1:numst,1:nt))				! Maximum value for those who do get an offer.
allocate (v_no_max(1:numst,1:nt))				! Maximum value for those who do not get an offer.
allocate (sumv_off_max(1:numst,1:nt))			! Sum over the maximum values for those who do get an offer.
allocate (sumv_no_max(1:numst,1:nt))			! Sum over the maximum values for those who do not get an offer.
allocate (Ev_off_max(1:numst,1:nt))			! Expected maximum value for those who do get an offer.
allocate (Ev_no_max(1:numst,1:nt))			! Expected maximum value for those who do not get an offer.

! SIMULATED UTILITY FUNCTIONS, WAGES, NON-LABOR INCOME, VALUE FUNCTIONS
! Used after the DPP is solved when we simulate individuals' choices.

allocate (sim_utfunc(0:1,1:nobs,1:ntdata))			! Utility  
allocate (sim_wage(1:nobs,1:ntdata))				! Wage
allocate (sim_nlincome(0:1,1:nobs,1:ntdata))			! NL income
allocate (sim_v(0:1,1:nobs,1:ntdata))				! Value function

! DATASET OF INITIAL CONDITIONS

  allocate (initdata_file(1:nvs,1:nobs))				! The dimensions of this dataset are the number of initial conditions (nvs=5)
										! and the number of unique individuals (nobs=1500).

! FOR CALCULATING MODEL STATISTICS OR OTHER NEEDED VARIABLES

  allocate (totwork(1:ntdata))				! In the model statistics part of the code this is a count of how many
								! people work in a given period.
  allocate (pctwork(1:ntdata))				! In the model statistics part of the code this is the percent of 
								! people who work in a given period.
  allocate (receivedoffer(1:nobs,1:ntdata))	! Indicator variable used in the simulations to denote whether
								! the person received a job offer regardless of whether they accept it.

! SET RANDOM NUMBER SEED FOR UNIFORM DISTRIBUTION
! This guarantees that every time we re-run the code that the same
! set of random numbers from the uniform distribution are called.
! We use these random numbers to determine the realization of
! whether an individual got a job offer (if he did not work in the prior period).

	CALL RANDOM_SEED(SIZE=sizer)
	allocate (joffseed(sizer))
	joffseed=1336187
	CALL RANDOM_SEED(PUT=joffseed)

! SET VARIABLES TO ZERO
! This guarantees all variables are set to a value.
! We do not do this for model parameters.
! We address those later.

work				= zero
ush				= zero
ush_sim			= zero
wsh				= zero
wsh_sim			= zero
exper			= zero
lastwork			= zero
educ				= zero
married			= zero
stock				= zero
utfunc			= zero
wage 			= zero
nlincome			= zero
v				= zero
v_off_max			= zero
v_no_max			= zero
sumv_off_max		= zero
sumv_no_max		= zero
Ev_off_max		= zero
Ev_no_max		= zero
sim_v			= zero
sim_utfunc			= zero
sim_wage			= zero
sim_nlincome		= zero
joffrn				= zero
totwork			= zero
pctwork			= zero
receivedoffer		= zero

end subroutine allocate_equations

!******************************************************************************
! This subroutine calculates the value functions.
!******************************************************************************

subroutine calculate_values
  use initialize_data
  implicit none 
    
    call allocate_equations				! Call the allocate equation subroutine.

!******************************************************************************
! Set the values for all of the model parameters.
!******************************************************************************

! EARNINGS OFFER PARAMETERS

w_intercept= 	10.1_prec				! Intercept.
w_xper= 		0.05_prec				! Coefficient on experience.
w_xpersq=		-0.0006_prec			! Coefficient on experience squared.
w_educ(0)= 	zero					! Coefficient on less than HS grad (omitted category).
w_educ(1)= 	0.26_prec				! Coefficient on HS grad.
w_educ(2)=	0.65_prec				! Coefficient on some college.
w_lastnw=		-0.13_prec				! Coefficient on not working last period.

! UTILITY PARAMETERS

u_work=		-2.65_prec				! Disutility from work.

! JOB OFFER PROBABILITY PARAMETERS
! Logit formulation assumed for job offer probability.
! Use the parameters above to generate the actual job offer probabilities for all feasible state space combinations.

do ed=0,2
	do mar=0,1
		offerprob(ed,mar)=(exp(offerprob_intercept + offerprob_educ(ed) + (offerprob_married*mar))) / &
					   (one + exp(offerprob_intercept + offerprob_educ(ed) + (offerprob_married*mar))) 
	end do
end do

! NON-LABOR INCOME PARAMETRS

nli_intercept=	7.50_prec
nli_work=		-0.85_prec
nli_educ(0)=	zero
nli_educ(1)=	-0.28_prec
nli_educ(2)=	-0.60_prec
nli_married=	0.25_prec

! DISCOUNT FACTOR

beta=0.95_prec

! VARIANCE AND COVARIANCE PARAMETERS

alpha_cov_u_1_1=	0.08_prec			! Parameter that governs variance of shock to working.
								! Note, this is not the value of the variance itself (see below).

alpha_var_wsh=		0.40_prec			! Parameter that governs variance of wage shock.
								! Note, this is not the value of the variance itself (see below).					

! CREATE AN INDEX FOR THE COMBINATIONS OF STATE VARIABLES
! Different combinations of state variables will be represented by a unique integer given below.
! Think of 'stock' as taking a combination of state variables and outputting an integer unique
! to that combination of state variables.
	
	i_g=1		! Specify the stock to start at 1.
	
	do xper=0,maxexp
		do ed=0,2
			do lw=0,1
				do mar=0,1
					stock(xper,ed,lw,mar)= i_g
					i_g=i_g+1			! Once we assigned a combination of states to a unique integer,
									! tick up the i_g index by 1 so that the next combination
									! will be assigned to a unique integer.
				end do
			end do
		end do
	end do
			
!******************************************************************
! Last period value functions first. 
! This is period nt.
!******************************************************************

! Generates the errors for this period.
! Note wage shock comes from log normal distribution.
	
 	 ISEED=123457							! This seed is different from the one we specified in allocate_equations.
										! This one governs the random number draws from the multivariate normal
										! and log normal distributions.  We could have set it earlier but it is fine to do so here.
										! The other one governs the random number draws from the uniform distribution
										! used to determine whether an individual got a job offer.

	 COV_U(1,1) = exp(alpha_cov_u_1_1)			! This is the variance of the shock to the utility from work.
										! Note that if there was another shock, the additional variance
										! and covariance terms would need to be set.		

	 do i=1,nsh_u,1
  	 	XMU_U(i)=0.0_prec					! Utility shock has mean zero.
  	 end do

	 XMU_W=0.0_prec						! Wage shock has mean zero.
    	 VAR_W=exp(alpha_var_wsh)				! Note VAR_W is the variance not standard error.
	 
  	 !CALL DRANDINITIALIZE(1,1,ISEED,1,STATE,LSTATE,INFO)
  	 !CALL DRANDMULTINORMAL(DT,NSH_U,XMU_U,COV_U,NSH_U,STATE,shock_u_val,DT,INFO)  
	 !CALL DRANDLOGNORMAL(DT,XMU_W,VAR_W,STATE,shock_w_val,INFO)   	 
     call r8vec_normal_ab(DT,XMU_U,sqrt(COV_U(NSH_U,NSH_U)),iseed,shock_u_val)
     call r8vec_normal_ab(DT,XMU_W,sqrt(VAR_W),iseed,shock_w_val)
     shock_w_val=exp(shock_w_val)

! Assigns each draw from the CALL commands above to its appropriate error. 

	do d=1,dt,1
		ush(d,nt)=shock_u_val(d,1)
		
		wsh(d,nt)=shock_w_val(d)
	end do

! ************************************
! Non-Labor Income
! ************************************

! Assigned here since it does not depend on the error draws and there is no need to rewrite it for each draw.

do a=0,1					! This loops over work or not decision.
	do xper=0,maxexp		! Even though experience is not in the NL income equation, we need to loop over it
						! because it's in the stock array.
		do ed=0,2	
			do lw=0,1
				do mar=0,1
						nlincome(a,stock(xper,ed,lw,mar),nt)= &
						exp(nli_intercept + (a*nli_work) + nli_educ(ed) + (mar*nli_married))
				end do
			end do
		end do
	end do
end do

do d=1,dt,1		! Looping over the error draws.
    
! Wages, utility functions, and value functions get over-written each draw.

! *************
! Wages
! *************

! We assume no one works in the terminal period so all wages are set to zero.
! This is a bit redundant and unnecessary since wage is set to zero in allocate_equations.

do xper=0,maxexp
	do ed=0,2
		do lw=0,1
			do mar=0,1
				wage(stock(xper,ed,lw,mar),nt)=zero
			end do
		end do
	end do
end do

! ********************
! Utility Function
! ********************

! For people who do not work.
! Since we assume all individuals do not work in the terminal period,
! we do not need to specify the utility function for those who work.
! Recall it was set to zero in allocate_equations.

! Utility is just log consumption in this case and consumption is just NL income.

do xper=0,maxexp
	do ed=0,2
		do lw=0,1
			do mar=0,1
				utfunc(0,stock(xper,ed,lw,mar),nt)=log(nlincome(0,stock(xper,ed,lw,mar),nt))
			end do
		end do
	end do
end do

! ***************************************
! Terminal Period Value Functions
! ***************************************

! a indexes work choice.
! We really do not need to solve the value function in this period
! for individuals who work.  But, there is no harm in doing so in this case.
! We'll end up never using it.

do a=0,1
	do xper=0,maxexp
		do ed=0,2
			do lw=0,1
				do mar=0,1
					v(a,stock(xper,ed,lw,mar),nt)= &
					utfunc(a,stock(xper,ed,lw,mar),nt)
				end do
			end do	
		end do
	end do
end do

! ********************************************
! Max Value for Different Choice Sets
! ********************************************

! The max value function associated with having an offer should be set to zero 
! since this cannot occur in the terminal period.
! Recall max value functions already set to zero in allocate_equations.

! Do not receive an offer. 
! No choice available, just does not work.

do xper=0,maxexp
	do ed=0,2
		do lw=0,1
			do mar=0,1
				v_no_max(stock(xper,ed,lw,mar),nt)=v(0,stock(xper,ed,lw,mar),nt)
			end do
		end do
	end do
end do

! Now sum up the above max value functions to take expected values over the draws.
! Note that this is still within the draw loop.

! The sum max value function associated with having an offer should be set to zero 
! since this cannot occur in the terminal period.
! Recall sum max value functions already set to zero in allocate_equations.

do xper=0,maxexp
	do ed=0,2
		do lw=0,1
			do mar=0,1	
				sumv_no_max(stock(xper,ed,lw,mar),nt)= &
				sumv_no_max(stock(xper,ed,lw,mar),nt) + &
				v_no_max(stock(xper,ed,lw,mar),nt)
			end do
		end do
	end do
end do

end do					! Ends the do loop over the draws.

! Now get the expected values for these sums.

! The Emax function associated with having an offer should be set to zero 
! since this cannot occur in the terminal period.
! Recall Emax functions set to zero in allocate_equations.

do xper=0,maxexp
	do ed=0,2
		do lw=0,1
			do mar=0,1
				Ev_no_max(stock(xper,ed,lw,mar),nt)=sumv_no_max(stock(xper,ed,lw,mar),nt) / dt
			end do
		end do
	end do
end do

! Resets the error holding arrays to zero.

	shock_u_val=zero
	shock_w_val=zero

!******************************************************************************
! Now, solve the earlier periods by backwards recursion.
!******************************************************************************
!******************************************************************************
! First do the period right before nt (i.e. nt-1).
! We solve this period separately from the remaining periods
! since we need to make sure the value function reflects 
! that offers are not received in the terminal period.
!******************************************************************************

! Generates the errors for this period.

  	 !CALL DRANDMULTINORMAL(DT,NSH_U,XMU_U,COV_U,NSH_U,STATE,shock_u_val,DT,INFO)  
	 !CALL DRANDLOGNORMAL(DT,XMU_W,VAR_W,STATE,shock_w_val,INFO)   	 
     call r8vec_normal_ab(DT,XMU_U,sqrt(COV_U(NSH_U,NSH_U)),iseed,shock_u_val)
     call r8vec_normal_ab(DT,XMU_W,sqrt(VAR_W),iseed,shock_w_val)
     shock_w_val=exp(shock_w_val)

 ! Assigns each draw to its appropriate error. 

	do d=1,dt,1
		ush(d,nt-1)= shock_u_val(d,1)
		
		wsh(d,nt-1)= shock_w_val(d)
	end do

! ************************************
! Non-Labor Income
! ************************************

! Assigned here since it does not depend on the error draws and there is no need to rewrite it for each draw.

do a=0,1					! This loops over work or not decision.
	do xper=0,maxexp		! Even though experience is not in the NL income equation, we need to loop over it
						! because it's in the stock ardr13ray.
		do ed=0,2	
			do lw=0,1
				do mar=0,1
						nlincome(a,stock(xper,ed,lw,mar),nt-1)= &
						exp(nli_intercept + (a*nli_work) + nli_educ(ed) + (mar*nli_married))
				end do
			end do
		end do
	end do
end do

do d=1,dt,1		! Looping over the error draws.
    
! Wages, utility functions, and value functions get over-written each draw.

! *************
! Wages
! *************

do xper=0,maxexp
	do ed=0,2
		do lw=0,1
			do mar=0,1
				wage(stock(xper,ed,lw,mar),nt-1)=exp(w_intercept + (w_xper*xper) + (w_xpersq*(xper**2)) + &
										w_educ(ed) + (w_lastnw*(1-lw)))*wsh(d,nt-1)	
			end do
		end do
	end do
end do

! ********************
! Utility Function
! ********************

! For people who do not work.
! Utility is just log consumption in this case and consumption is just NL income.

do xper=0,maxexp
	do ed=0,2
		do lw=0,1
			do mar=0,1
				utfunc(0,stock(xper,ed,lw,mar),nt-1)=log(nlincome(0,stock(xper,ed,lw,mar),nt-1))
			end do
		end do
	end do
end do

! For people who work.

do xper=0,maxexp
	do ed=0,2
		do lw=0,1
			do mar=0,1
				utfunc(1,stock(xper,ed,lw,mar),nt-1)= &
				log(nlincome(1,stock(xper,ed,lw,mar),nt-1) + wage(stock(xper,ed,lw,mar),nt-1)) + &
				u_work + ush(d,nt-1)
			end do
		end do
	end do
end do

! ***************************************
! Value Functions
! ***************************************

! Regardless of whether someone works or not in this period,
! they will not get an offer in the following period
! because we assume no one gets an offer or works in the terminal period.

! Note how experience accumulates in nt if the individual works in nt-1.

do a=0,1
	do xper=0,maxexp
		do ed=0,2
			do lw=0,1
				do mar=0,1
					v(a,stock(xper,ed,lw,mar),nt-1)= &
					utfunc(a,stock(xper,ed,lw,mar),nt-1) + &
					beta*Ev_no_max(stock(min(maxexp,xper+a),ed,a,mar),nt)
				end do
			end do	
		end do
	end do
end do

! ********************************************
! Max Value for Different Choice Sets
! ********************************************

! Do not receive an offer. 
! No choice available, just does not work.
! This should only occur for those who did not work last period.
! So don't loop over lw, just replace it with 0.

do xper=0,maxexp
	do ed=0,2
		do mar=0,1
			v_no_max(stock(xper,ed,0,mar),nt-1)=v(0,stock(xper,ed,0,mar),nt-1)
		end do
	end do
end do

! Receives an offer. 
! Chooses to work or not.
! Could occur for both those who did and did not work last period.
! So loop fully over lw.

do xper=0,maxexp
	do ed=0,2
		do lw=0,1
			do mar=0,1
				v_off_max(stock(xper,ed,lw,mar),nt-1)=max(v(0,stock(xper,ed,lw,mar),nt-1),v(1,stock(xper,ed,lw,mar),nt-1))
			end do
		end do
	end do
end do

! Now sum up the above max value functions to take expected values over the draws.
! Note that this is still within the draw loop.
! Note that except in the terminal period, it should never be the case
! that an individual worked last period and does not get an offer.
! We loop over these possibilities anyway for convenience.

do xper=0,maxexp
	do ed=0,2
		do lw=0,1
			do mar=0,1	
				sumv_no_max(stock(xper,ed,lw,mar),nt-1)= &
				sumv_no_max(stock(xper,ed,lw,mar),nt-1) + &
				v_no_max(stock(xper,ed,lw,mar),nt-1)

				sumv_off_max(stock(xper,ed,lw,mar),nt-1)= &
				sumv_off_max(stock(xper,ed,lw,mar),nt-1) + &
				v_off_max(stock(xper,ed,lw,mar),nt-1)
			end do
		end do
	end do
end do

end do					! Ends the do loop over the draws.

! Now get the expected values for these sums.
! Note that except in the terminal period, it should never be the case
! that an individual worked last period and does not get an offer.
! We loop over these possibilities anyway for convenience.
! We'll end up never using those specific Emax's.

do xper=0,maxexp
	do ed=0,2
		do lw=0,1
			do mar=0,1
				Ev_no_max(stock(xper,ed,lw,mar),nt-1)=sumv_no_max(stock(xper,ed,lw,mar),nt-1) / dt

				Ev_off_max(stock(xper,ed,lw,mar),nt-1)=sumv_off_max(stock(xper,ed,lw,mar),nt-1) / dt
			end do
		end do
	end do
end do




! Resets the error holding arrays to zero.

	shock_u_val=zero
	shock_w_val=zero

!******************************************************************************
! Now solve all remaining periods by backward recursion.
!******************************************************************************

do t=nt-2,1,-1				! The -1 says to do this do loop in reverse from time nt-2 to time 1.

	!CALL DRANDMULTINORMAL(DT,NSH_U,XMU_U,COV_U,NSH_U,STATE,shock_u_val,DT,INFO)  
	!CALL DRANDLOGNORMAL(DT,XMU_W,VAR_W,STATE,shock_w_val,INFO)   	 
    call r8vec_normal_ab(DT,XMU_U,sqrt(COV_U(NSH_U,NSH_U)),iseed,shock_u_val)
    call r8vec_normal_ab(DT,XMU_W,sqrt(VAR_W),iseed,shock_w_val)
    shock_w_val=exp(shock_w_val)

 ! Assigns each draw to its appropriate error. 

	do d=1,dt,1
		ush(d,t)=shock_u_val(d,1)
		
		wsh(d,t)=shock_w_val(d)
	end do

! ************************************
! Non-Labor Income
! ************************************

! Assigned here since it does not depend on the error draws and there is no need to rewrite it for each draw.

do a=0,1					! This loops over work or not decision.
	do xper=0,maxexp		! Even though experience is not in the NL income equation, we need to loop over it
						! because it's in the stock array.
		do ed=0,2	
			do lw=0,1
				do mar=0,1
						nlincome(a,stock(xper,ed,lw,mar),t)= &
						exp(nli_intercept + (a*nli_work) + nli_educ(ed) + (mar*nli_married))
				end do
			end do
		end do
	end do
end do

do d=1,dt,1		! Looping over the error draws.
    
! Wages, utility functions, and value functions get over-written each draw.

! *************
! Wages
! *************

do xper=0,maxexp
	do ed=0,2
		do lw=0,1
			do mar=0,1
				wage(stock(xper,ed,lw,mar),t)=exp(w_intercept + (w_xper*xper) + (w_xpersq*(xper**2)) + &
										w_educ(ed) + (w_lastnw*(1-lw)))*wsh(d,t)	
			end do
		end do
	end do
end do

! ********************
! Utility Function
! ********************

! For people who do not work.
! Utility is just log consumption in this case and consumption is just NL income.

do xper=0,maxexp
	do ed=0,2
		do lw=0,1
			do mar=0,1
				utfunc(0,stock(xper,ed,lw,mar),t)=log(nlincome(0,stock(xper,ed,lw,mar),t))
			end do
		end do
	end do
end do

! For people who work.

do xper=0,maxexp
	do ed=0,2
		do lw=0,1
			do mar=0,1
				utfunc(1,stock(xper,ed,lw,mar),t)= &
				log(nlincome(1,stock(xper,ed,lw,mar),t) + wage(stock(xper,ed,lw,mar),t)) + &
				u_work + ush(d,t)
			end do
		end do
	end do
end do

! ***************************************
! Value Functions
! ***************************************

! If the individual does not work this period,
! they face a probability, offerprob, of getting an offer,
! and a probability, 1-offerprob, of not getting an offer in t+1.

! If the individual does work this period,
! he gets an offer with certainty next period.
! Note how experience accumulates entering period t+1.

do xper=0,maxexp
	do ed=0,2
		do lw=0,1
			do mar=0,1
				v(0,stock(xper,ed,lw,mar),t)= &
				utfunc(0,stock(xper,ed,lw,mar),t) + &
				beta*((offerprob(ed,mar)*Ev_off_max(stock(xper,ed,0,mar),t+1)) + &
				((one-offerprob(ed,mar))*Ev_no_max(stock(xper,ed,0,mar),t+1)))

				v(1,stock(xper,ed,lw,mar),t)= &
				utfunc(1,stock(xper,ed,lw,mar),t) + &
				beta*Ev_off_max(stock(min(maxexp,xper+1),ed,1,mar),t+1)
			end do	
		end do
	end do
end do

! ********************************************
! Max Value for Different Choice Sets
! ********************************************

! Do not receive an offer. 
! No choice available, just does not work.
! This should only occur for those who did not work last period.
! So don't loop over lw, just replace it with 0.

do xper=0,maxexp
	do ed=0,2
		do mar=0,1
			v_no_max(stock(xper,ed,0,mar),t)=v(0,stock(xper,ed,0,mar),t)
		end do
	end do
end do

! Receives an offer. 
! Chooses to work or not.
! Could occur for both those who did and did not work last period.
! So loop fully over lw.

do xper=0,maxexp
	do ed=0,2
		do lw=0,1
			do mar=0,1
				v_off_max(stock(xper,ed,lw,mar),t)=max(v(0,stock(xper,ed,lw,mar),t),v(1,stock(xper,ed,lw,mar),t))
			end do
		end do
	end do
end do

! Now sum up the above max value functions to take expected values over the draws.
! Note that this is still within the draw loop.
! Note that except in the terminal period, it should never be the case
! that an individual worked last period and does not get an offer.
! We loop over these possibilities anyway for convenience.

do xper=0,maxexp
	do ed=0,2
		do lw=0,1
			do mar=0,1	
				sumv_no_max(stock(xper,ed,lw,mar),t)= &
				sumv_no_max(stock(xper,ed,lw,mar),t) + &
				v_no_max(stock(xper,ed,lw,mar),t)

				sumv_off_max(stock(xper,ed,lw,mar),t)= &
				sumv_off_max(stock(xper,ed,lw,mar),t) + &
				v_off_max(stock(xper,ed,lw,mar),t)
			end do
		end do
	end do
end do

! Resets the error holding array to zero.

	shock_u_val=zero
	shock_w_val=zero

end do					! Ends the do loop over the draws.

! Now get the expected values for these sums.
! Note that except in the terminal period, it should never be the case
! that an individual worked last period and does not get an offer.
! We loop over these possibilities anyway for convenience.
! We'll end up never using those specific Emax's.

do xper=0,maxexp
	do ed=0,2
		do lw=0,1
			do mar=0,1
				Ev_no_max(stock(xper,ed,lw,mar),t)=sumv_no_max(stock(xper,ed,lw,mar),t) / dt

				Ev_off_max(stock(xper,ed,lw,mar),t)=sumv_off_max(stock(xper,ed,lw,mar),t) / dt
			end do
		end do
	end do
end do

end do	  ! Ends loop over time		


!print*,"t,edu,exp,lw,mar,EV_no,EV_of"
!do t = 1,nt
!do ed=0,2
!	do xper=0,maxexp
!		do lw=0,1
!			do mar=0,1
!				if (xper<= t-1) print*,t,ed,xper,lw,mar,Ev_no_max(stock(xper,ed,lw,mar),t),Ev_off_max(stock(xper,ed,lw,mar),t)
!			end do
!		end do
!	end do
!end do
!end do

end subroutine calculate_values

!*****************************************************************************
! We now have the Emax's needed to simulate individuals'
! decisions given realized draws of shocks.
! Below we simulate data using the Emax's calculated above.
! One might do such an exercise for model fit, estimation testing, or
! debugging to make sure the model is working correctly.
!*****************************************************************************

subroutine simulations
  use initialize_data
  implicit none

call calculate_values				! This calls the calculate_values subroutine which needs to be done first
							! since we need the Emax's for the simulations.

! All stochastic elements generated before the start of the simulations.
! Note the shocks we are drawing now are not for the DPP solution.
! We have already done that.
! These shocks correspond to an individual's actual realization of 
! utility and wage shocks.

do t=1,ntdata,1

!	 CALL DRANDMULTINORMAL(nobs,NSH_U,XMU_U,COV_U,NSH_U,STATE,shock_u_sim,nobs,INFO) 
!	 CALL DRANDLOGNORMAL(nobs,XMU_W,VAR_W,STATE,shock_w_sim,INFO) 
     call r8vec_normal_ab(nobs,XMU_U,sqrt(COV_U(NSH_U,NSH_U)),iseed,shock_u_sim)
     call r8vec_normal_ab(nobs,XMU_W,sqrt(VAR_W),iseed,shock_w_sim)
     shock_w_sim=exp(shock_w_sim)

! Assigns the draws to the appropriate person.

 	do n=1,nobs,1
	
	ush_sim(n,t)=shock_u_sim(n,1)

	wsh_sim(n,t)=shock_w_sim(n)

	end do 			! Ends loop over people within a period.
end do 				! Ends loop over time.

! Call the job offer shocks.
! We'll use these to determine if a person got an offer or not.

do n=1,nobs,1
	do t=1,ntdata,1
			CALL RANDOM_NUMBER(joffrn(n,t))				! This calls a random number from a uniform distribution.
	end do		
end do						! All stochastic elements for simulations drawn at this time.

! Here is where the data of initial conditions needs to be read in.
! This is a comma delimited file with nvs columns and nobs rows.
! Note you will need to change the directory to match the one you use.
! You can read about the next 3 lines from the links I sent you.
! Basically these are things Fortran requires you to specify when
! reading in a dataset.  Recall we set init_data earlier in the code.

open(unit=init_data, file="InitData_BasicDPP_1.txt", status="old", form="formatted", access="sequential")
read(init_data,*)initdata_file
close(unit=init_data)

do i=1,nobs,1				! Begins loop over people.
						! Here we assign elements from the initial conditions dataset to their appropriate variables.

n=initdata_file(1,i)
exper(n,1)=initdata_file(2,i)
educ(n)=initdata_file(3,i)
lastwork(n,1)=initdata_file(4,i)
married(n,1)=initdata_file(5,i)

end do					! Ends loop over people.

! At this point, every individual's initital conditions have been read in.
! Now we want to simulate their choices over time and record their
! state variables over time.

do n=1,nobs				! We loop over individuals
	do t=1,ntdata,1			! and over time periods, starting at t=1.
						! Note, we move forward in time because backward recursion is only needed
						! for solving the DPP, which we did.
						! We're now seeing what individuals actually choose to do in each period 
						! as time moves forward.

	! Here we specify how state variables evolve across periods.
		
	if (t>1) then							! If this is not the initial period...
		married(n,t)=married(n,t-1)			! We are assuming marital status does not change over time.

		lastwork(n,t)=work(n,t-1)			! Last work is simply the work decision from the prior period.

		exper(n,t)=min(maxexp,exper(n,t-1)+work(n,t-1))	! Experience is last period's experience plus 1 if the individual worked,
											! with the cap on experience imposed.	
	end if	


! Assign wages, non-labor income, utility functions, and value functions based on each individual's
! realized state variables entering the period.

! *************
! Wages
! *************

sim_wage(n,t)=exp(w_intercept + (w_xper*exper(n,t)) + (w_xpersq*(exper(n,t)**2)) + w_educ(educ(n)) + &
					(w_lastnw*(1-lastwork(n,t))))*wsh_sim(n,t)
				
! **********************
! Non-Labor Income
! **********************

do a=0,1
	sim_nlincome(a,n,t)=exp(nli_intercept + (a*nli_work) + nli_educ(educ(n)) + (nli_married*married(n,t))) 
end do

! ********************
! Utility Function
! ********************

! For those who do not work.

sim_utfunc(0,n,t)= log(sim_nlincome(0,n,t))

! For those who work.

sim_utfunc(1,n,t)=log(sim_nlincome(1,n,t) + sim_wage(n,t)) + u_work + ush_sim(n,t)

! ***********************************
! Value Functions
! ***********************************

! Terminal period value functions.

if (t==nt) then

	! No one works in the terminal period but we solve the simulated value function
	! for that choice anyway.
	! Later, we will impose that working is not part of the choice set in the terminal period.

do a=0,1
	sim_v(a,n,t)=sim_utfunc(a,n,t)
end do

end if			! Ends if over terminal period.

! One year before terminal period.

! Recall, regardless of whether the individual works or not,
! he will not get an offer in the terminal period.

if (t==nt-1) then

	sim_v(0,n,t)=sim_utfunc(0,n,t) + &
	beta*Ev_no_max(stock(exper(n,t),educ(n),0,married(n,t)),t+1)

	sim_v(1,n,t)=sim_utfunc(1,n,t) + &
	beta*Ev_no_max(stock(min(maxexp,exper(n,t)+1),educ(n),1,married(n,t)),t+1)

end if			! Ends if over one year before terminal period.

! Now all remaining periods.

if (t<nt-1) then

	sim_v(0,n,t)=sim_utfunc(0,n,t) + &
	beta*((offerprob(educ(n),married(n,t))*Ev_off_max(stock(exper(n,t),educ(n),0,married(n,t)),t+1)) + &
	((one-offerprob(educ(n),married(n,t)))*Ev_no_max(stock(exper(n,t),educ(n),0,married(n,t)),t+1)))

	sim_v(1,n,t)=sim_utfunc(1,n,t) + &
	beta*(Ev_off_max(stock(min(maxexp,exper(n,t)+1),educ(n),1,married(n,t)),t+1)) 

end if			! Ends if over all periods before nt-1.

!*********************************************************************
! This part of the program has simulated individuals choose the option 
! available to them that is associated with the highest value function.     
! We also track whether an individual receives an offer or not, 
! keeping in mind receivedoffer was set to zero in allocate_equations. 
!*********************************************************************

! If terminal period, make sure no one works.

if (t==nt) then
	work(n,t)=0
end if

if (t<nt) then

	if (lastwork(n,t)==1) then			! If the person worked last period,
								! they have an offer this period and choose
								! between working or not.

		do a=0,1
			if (sim_v(a,n,t)==max(sim_v(0,n,t),sim_v(1,n,t))) then
				work(n,t)=a
				receivedoffer(n,t)=1
			end if
		end do

	end if

	if (lastwork(n,t)==0) then			! If the person did not work last period,
								! their choice set this period depends on whether they
								! get an offer or not.  This is where those random numbers
								! we drew earlier come into play.
								! If the random number joffrn is less than or equal to the probability of 
								! receiving an offer, the individual receives an offer.

		if (joffrn(n,t)<=offerprob(educ(n),married(n,t))) then			! Gets an offer.
			do a=0,1
				if (sim_v(a,n,t)==max(sim_v(0,n,t),sim_v(1,n,t))) then
					work(n,t)=a
					receivedoffer(n,t)=1
				end if
			end do
		end if

		if (joffrn(n,t)>offerprob(educ(n),married(n,t))) then			! Does not get an offer.
			work(n,t)=0
		end if

	end if

end if						! Ends if over non-terminal periods.

	end do				! Ends do over time.
end do 					! Ends do over people.


!*********************************************************************
! Writes the data generated from this program to a
! comma-delimited .txt file.
!**********************************************************************

! The i's and F's below correspond to the type of number contained in each variable.
! For example, lastwork(n,t) is a 1-digit integer.  
! sim_wage(n,t) is a float where we record 10 digits including 2 decimal points.
! The A1's correspond to the commas.

open(unit=sim_data,file="SimData_BasicDPP_1.txt",form="formatted", status="unknown")
do n=1,nobs,1
do t=1,ntdata,1
	write(sim_data,'(i5,A1,i2,A1,i1,A1,i1,A1,i1,A1,i1,A1,i1,A1,i1,A1,F10.2,A1,F10.2)')n &
				,',',t,',',work(n,t),',',exper(n,t),',',educ(n),&
				',',lastwork(n,t),',',married(n,t),',',receivedoffer(n,t)&
				,',',sim_wage(n,t),',',sim_nlincome(work(n,t),n,t)
end do
end do
close(unit=sim_data)

end subroutine simulations

!*******************************************************************
! This part of the program tabulates some basic model statistics.
!*******************************************************************

subroutine model_statistics
  use initialize_data
  implicit none

! Tally up the number of individuals that work per period.

do n=1,nobs
	do t=1,ntdata

		if (work(n,t)==1) then
			totwork(t)=totwork(t)+one
		end if

	end do
end do

! Now create a variable for the percent of individuals who work per period.

do t=1,ntdata
	pctwork(t)=totwork(t)/nobs
end do

! Write the percent working per period to the screen.

do t=1,ntdata
	write(*,15)'Pct working in year ',t,' ',pctwork(t)
end do

15 format(1x,a60,i1,a3,f5.3)

end subroutine model_statistics

!****************************
! Main Program
!****************************

program basicdpp

  use initialize_data
  implicit none

external allocate_equations, calculate_values, simulations, model_statistics

! Calculate the value functions, do the simulations, and compute model statistics.
! Recall that the simulations subroutine calls the calculate_values subroutine
! which calls the allocate_equations subroutine.  So we only need to call 
! the simulations and the model_statistics subroutines here.

call simulations
call model_statistics

end program basicdpp
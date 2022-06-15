//  The model specific functions for restidental and location model
//  Is included in rwmodel.c
//  Arguments of the majority of the functions are defined in declSTATES. When calling these functions use STATES. 
//  all values are accessible through p[0] as descibed above. 
	
//*********************************************************************
//Declarations of all the of the functions below
//*********************************************************************
double f_inc(const param* p, declSTATES, declOUTCOMES);
double f_aftertax(const param* p, const double inc, const int irh);
double f_get_inc(const param* p, declSTATES, declOUTCOMES);
void   f_precomp_inc(param* p);
double f_kappa(const param* p, declSTATES, declOUTCOMES, const double inc);
double f_phi(const param* p, declSTATES, declOUTCOMES, const double inc);
double f_sqm(const param* p, declSTATES, declOUTCOMES);
double f_psqm(const param* p, declOUTCOMES);
double f_umvcost_rl(const param* p, declSTATES, declOUTCOMES);
double f_umvcost_wl(const param* p, declSTATES, declOUTCOMES);
double f_mswcost(const param* p, declSTATES, declOUTCOMES);
double f_ttimecost(const param* p, declOUTCOMES); 
double f_regtaste(const param* p, declSTATES, declOUTCOMES); 
double f_u(const param* p, declSTATES, declOUTCOMES);

// probability of succeeding in getting a new job. 
double f_pn(const param* p, declSTATES, const int irw){
	// irw is the intention of getting a job (choice)
	double pn, xb, age; 
	
	if (irw ==(p->nrwork-1)) return 1; // Always possible to choose unemployment

	// if (irw == irwp) return 0; // if trying to just keep current job, agent doesn't search for any other jobs

	age = (double) (it + p->t0);
	
	xb =  p->Bpn_0 + p->Bpn_a*age; 	
	
	if (is != 0) xb=xb+p->Bpn_s[is-1]; // educational dummies
	
	if (irwp ==(p->nrwork-1))  	xb=xb+p->Bpn_unemp; 	// Dummy for unemployment in last period
	
	xb = xb + p->Bpn_jobdensity * p->jobdensity[irw + (p->nrwork-1) * is]; // Effect of jobdensity. 
	
	if (irwp != (p->nrwork-1)) xb = xb + p->Bpn_ttime_rwp*p->ttime[f_IDCHOICE(p, irw,irwp)];
	
	xb = xb + p->Bpn_ttime_rhp*p->ttime[f_IDCHOICE(p, irw,irhp)]; 
	// if (irw!= (p->nrwork-1))  	xb=xb+p->Bpn_rw[irw]; // Regional dummies (chosen work region)

	pn=1/(1+exp(-xb));
	return pn;
}

// probability of keeping current job. 
double f_pk(const param* p, declSTATES){
	
	double pk, xb, age; 
	if (irwp == (p->nrwork-1)) return 1; // Always possible to stay unemployed

	age = (double) (it + p->t0);
	
	xb = p->Bpk_0 + p->Bpk_a*age; 	
	
	if (is != 0)  xb=xb+p->Bpk_s[is-1]; // educational dummies
	
	xb = xb + p->Bpk_jobdensity * p->jobdensity[irwp + (p->nrwork-1) * is]; // Effect of jobdensity. 
	
	xb = xb + p->Bpk_ttime_rhp*p->ttime[f_IDCHOICE(p, irwp, irhp)];
	
	// if (irwp!= (p->nrwork-1))  	xb=xb+p->Bpk_rw[irwp]; // Regional dummies (chosen work region)

	pk=1/(1+exp(-xb)); 
	return pk;
}

double f_inc(const param* p, declSTATES, declOUTCOMES){
	// Expected income given age, region of work, income type, schooling and unemployment status. 
	// Region specific effects are for regions only (without non-employment)!

	double lninc;
	double tmp;
	double inc, age_poly;
	int idx_rw, idx_age_fe, age, year, poly_dg, age_p, i, j;

	// Age of individual
	age = it + p->t0;

	// Degree of polynomial in age
	poly_dg = p->n_beta_age_pol / p->ns;
	age_poly = 0.0;
	lninc = 0.0;

	// Income when employed
	// Wage equation is estimated in 1st stage in log(inc). 	
	if (irw != (p->nrwork-1))
	{

		// Index for work region
		idx_rw = irw + is*p->nrhome;
		
		// Add intercept regional effects and cohort effects
		lninc = p->beta_cons[is]
			+ p->beta_rw[idx_rw];

		// Unemployment in last period penalizes wages in current period.
		if (irwp == (p->nrwork-1))
		{
			lninc = lninc + p->beta_l_unemp[is]
					+ p->beta_l_unemp_rw[idx_rw];
		}

		// Adding age fixed effects if specified income regression permits them
		if (p->n_beta_age_fe > 0)
		{
			idx_age_fe = it + p->t0 - p->beta_min_age_fe + is*p->n_beta_age_fe/p->ns;
			lninc = lninc + p->beta_age_fe[idx_age_fe];
		}

		// Create polynomial in age recursively
		if (p->n_beta_age_pol > 0)
		{
			age_p = 1;
			for (i = 0; i < poly_dg; ++i)
			{
				age_p = age * age_p;
				age_poly = age_poly + age_p * p->beta_age_pol[i + is*poly_dg];
			}

			lninc = lninc + age_poly;
		}
	}
	else // When receiving non-work income
	{
		lninc = p->b_0[is];

		// Adding age fixed effects if specified income regression permits them
		if (p->n_b_age_fe > 0)
		{
			idx_age_fe = (int)(it + p->t0 - p->b_min_age_fe + is*(p->n_b_age_fe/p->ns));

			lninc = lninc + p->b_age_fe[idx_age_fe];
		}

		// Create polynomial in age recursively
		if (p->n_b_age_pol > 0)
		{
			age_p = 1;
			for (j = 0; j < poly_dg; ++j)
			{
				age_p = age * age_p;
				age_poly = age_poly + age_p * p->b_age_pol[j + is*poly_dg];
			}

			lninc = lninc + age_poly;
		}
	}

	if (lninc > 17)
		mexErrMsgTxt("f_inc: lninc is out of realistic bounds");

	inc = exp(lninc) / p->dkk_scale;

	return inc;
}

double f_aftertax(const param* p, const double inc, const int irh){
	// returns the after tax income, conditional on the region of residence
	// uncomment all six print statements to check the operation

	// return inc; // return income before tax (for testing)

	int ibr;
	double th, rt = 0, tax = 0, thp = 0;

	// printf("\nAfter tax calculator for income = %1.5f\n",inc);

	for (ibr = 0; ibr < p->ntaxbracket; ibr++) 
	{	// Loop over tax brackets

		th = p->thresholds[ibr*p->nrhome + irh];

		if (inc<th || ibr==p->ntaxbracket-1)
		{	// income below current threshold, or above maximum threshold
			tax += rt*(inc-thp);
			// printf("ibr=%d [thp,th]=[%1.3f,%1.3f] rt=%1.4f >> ",ibr,thp,th,rt);
			// printf("tax + remainder = %1.2f\n",tax);
			break;
		}
		else
		{	// income above current threshold
			tax += rt*(th-thp);
			// printf("ibr=%d [thp,th]=[%1.3f,%1.3f] rt=%1.4f >> ",ibr,thp,th,rt);
			// printf("tax + full bracket = %1.2f\n",tax);
		}

		rt = p->rates[ibr*p->nrhome + irh];
		thp = th;
	}
	// printf("After tax income = %1.5f (average tax rate %1.3f%%)\n",inc-tax,100*tax/inc);

	return inc-tax;
}

double f_get_inc(const param* p, declSTATES, declOUTCOMES){
    // Retrieve precomputed income by state-choice combination
    // Compute tax and return after tax income

    // return f_aftertax (p,p->inc[f_inc_idx(p, STATES, OUTCOMES)],irh);
    return f_aftertax (p,f_inc(p, STATES, OUTCOMES),irh);
}

void f_precomp_inc(param* p){
	// Precomputes before tax incomes for the given year.
	// These are cleaned for each iteration over years in evaluate_fct.

	// DANGER NOTE: if loop structure is changed, then according changes must be made in f_inc_idx in definitions.c  
	// DANGER NOTE: if f_inc is changed with respect to influencing state variables, then the loops below must be updated accordingly.

	int work, age, school, d_emp, emp_p, i, n_age;

	n_age = (p->T - p->t0 + 1);
	i = 0;

	// Looping over variables influencing income: 
	// work region choice -> age -> schooling -> nonemp status (binary)
	for (work = 0; work < p->nrwork; ++work)
	{
		for (age = 0; age < n_age; ++age)
		{
			for (school = 0; school < p->ns; ++school)
			{
				for (d_emp = 0; d_emp < 2; ++d_emp)
				{
					if (d_emp == 0)
						emp_p = 0;
					else
						emp_p = p->nrwork-1; // Non-employment in previous period influences income.

					// Compute income by standard function
					p->inc[i] = f_inc(p, age, 0, emp_p, 0, 0, 0, school, 0, 0, 0, work);	

					++i;
				}				
			}
		}
	}
}

double f_kappa(const param* p, declSTATES, declOUTCOMES, const double inc){
	// Marginal utility of money given income and macro state.
	// inc:		income net of monetary expenditures of moving and housing costs.

	double kappa;
	kappa = p->kappa_0 + p->kappa_y*inc 	+ p->kappa_y2*inc*inc; 
	
	if (p->reduced_form_sqm) // kappa holds reduced form parameters (alpha_P) that must be scaled
		return kappa=kappa*(SCALE_PHIH2*p->phi_h2)/p->psi_uc; 		// kappa=alpha_P*phi_h2/uc
	return kappa;
}

double f_phi(const param* p, declSTATES, declOUTCOMES, const double inc){
	// heterogeneous parameter, phi, affecting the baseline marginal utility of housing 
	double age, school, phi;
	age = (double) (it + p->t0);

	phi =  p->phi_0 	+ 	p->phi_a*age 	+  p->phi_a2*age*age/1000  
					 	+  	p->phi_y*inc 	+  p->phi_y2*inc*inc;


	if (icouple ==1)  	phi=phi+p->phi_ms; 
	if (ikids!= 0)  	phi=phi+p->phi_c[ikids-1]; //0 is base category, coef are indexed with -1
	// is stands for index of schooling
 	if (is!= 0)  		phi=phi+p->phi_s[is-1]; //0 is base category, coef are indexed with -1
 	if (irh!= 0)  		phi=phi+p->phi_rh[irh-1]; //0 is base category, coef are indexed with -1

 	if (p->reduced_form_sqm)
		return phi=-phi*(SCALE_PHIH2*p->phi_h2); // phi_h=-alpha_0*phi_h2 
	return phi;
}

double f_psqm(const param* p, declOUTCOMES){
	// Price per squre meter costs measured in yearly user costs of ownership.
	// Prices when parsed to the model are scaled by mp.dkk_scale (eg 100,000). 
	// The scaling takes place in rwmode.init(mp)
	return p->psi_uc * p->prices[irh];
}

double f_umvcost_rl(const param* p, declSTATES, declOUTCOMES){	
	// Utility moving costs, depending on couple status and kids. 
	
	double age, school, gamma;
	
	if (irhp == irh) // moving cost if staying
		return 0.0; 

	age = (double) (it + p->t0);
	
	gamma = p->gamma_0 + p->gamma_a*age +  p->gamma_a2*age*age/1000;
	if (icouple ==1)  	gamma=gamma+p->gamma_ms; 
	if (ikids!= 0)  	gamma=gamma+p->gamma_c[ikids-1]; 
 	if (is!= 0)  		gamma=gamma+p->gamma_s[is-1];  
 	
	return gamma;
}

double f_umvcost_wl(const param* p, declSTATES, declOUTCOMES){
	
	double pk, rho, age; 
	age = (double) (it + p->t0);

	if (irwp == irw)  // no moving cost if staying in job 
		return 0.0; 

	rho = p->rho_0 
		+ p->rho_a*age; 	
	if (is != 0)  					rho=rho+p->rho_s[is-1];  // educational dummies
	if (irwp ==(p->nrwork-1))  	rho=rho+p->rho_unemp;    // Dummy for moving out of unemployment
	
	// cost of searching job far from work location
	if ((irw !=(p->nrwork-1)) & (irwp !=(p->nrwork-1)))
		rho = rho + p->rho_ttime_rwp*p->ttime[f_IDCHOICE(p, irw,irwp)];

	// cost of searching job far from home location
	if ((irw !=(p->nrwork-1)) & (irhp !=(p->nrwork-1))) //second check is obsolete
		rho = rho + p->rho_ttime_rhp*p->ttime[f_IDCHOICE(p, irw,irhp)];
	// if (irw!= (p->nrwork-1))  	rho=rho+p->rho_rw[irw]; // Regional dummies (chosen work region)

	return rho;
}

double f_mswcost(const param* p, declSTATES, declOUTCOMES){
	// Monetary costs of moving. A fraction of the price of the previous home to reflect selling fees,
	// and a fraction of the new home to reflect transaction costs - lawyers, registration etc. 
	
	double mswcost;
	mswcost  = (irhp != irh)*(p->gamma_sell*p->prices[irhp] + p->gamma_buy*p->prices[irh]);
	return mswcost;
}

double f_ttimecost(const param* p, declOUTCOMES){

	double ttimecost; /* declare type */
	ttimecost = p->eta_ttime*p->ttime[IDCHOICE]+p->eta_ttime2*p->ttime[IDCHOICE]*p->ttime[IDCHOICE];
	return ttimecost;
}

double f_regtaste(const param* p, declSTATES, declOUTCOMES){
	// Regional specific taste-part of utility.
	double regtaste = 0.0; 

	regtaste =   regtaste 
			   // + p->alpha_dining*p->dining[irh];
			    + p->alpha_dining*p->dining[irh]
				+ p->alpha_cafes*p->cafes[irh]
				+ p->alpha_cafes_n*p->cafes_n[irh]
				+ p->alpha_cafes_sqm*p->cafes_sqm[irh]
				+ p->alpha_dining_n*p->dining_n[irh]
				+ p->alpha_dining_sqm*p->dining_sqm[irh]
				+ p->alpha_hotels*p->hotels[irh]
				+ p->alpha_hotels_n*p->hotels_n[irh]
				+ p->alpha_hotels_sqm*p->hotels_sqm[irh]
				+ p->alpha_s0*p->s0[irh]
				+ p->alpha_s0_dens*p->s0_dens[irh]
				+ p->alpha_s1*p->s1[irh]
				+ p->alpha_s1_dens*p->s1_dens[irh]
				+ p->alpha_s2*p->s2[irh]
				+ p->alpha_s2_dens*p->s2_dens[irh];

	if (irh!= 0)  regtaste = regtaste+p->alpha_rh[irh-1]; //0 is base category, coef are indexed with -1

	return regtaste;  
}

double f_sqm(const param* p, declSTATES, declOUTCOMES){
	// not part of utility, but part of solution (FIXME: CHECK)
	// Demand for housing in square meters given states. 

	double inc, kappa, phi, sqm;
	inc 	=  	f_get_inc(p, STATES, OUTCOMES);
 	kappa 	=  	f_kappa(p, STATES, OUTCOMES, inc);
	phi 	=  	f_phi(p, STATES, OUTCOMES, inc);
	sqm  	= 	MAX(0, 	((kappa*f_psqm(p, OUTCOMES)-phi)/(SCALE_PHIH2*p->phi_h2)));

	return sqm;
}

double f_u_sqm(const param* p, declSTATES, declOUTCOMES){
	// only used for output in evaluate_fct... not par of solution
	// utility of housing. 

	double inc, kappa, phi, sqm ;
	inc 	=  	f_get_inc(p, STATES, OUTCOMES);
 	kappa 	=  	f_kappa(p, STATES, OUTCOMES, inc);
	phi 	=  	f_phi(p, STATES, OUTCOMES, inc);
	sqm  	= 	MAX(0, 	((kappa*f_psqm(p, OUTCOMES)-phi)/(SCALE_PHIH2*p->phi_h2)));
return phi*sqm + (SCALE_PHIH2*p->phi_h2)/2*sqm*sqm;
}

double f_u_money(const param* p, declSTATES, declOUTCOMES){
	// only used for output in evaluate_fct... not par of solution
	// utility of money 
	double inc, kappa, phi, sqm ;
	inc 	=  	f_get_inc(p, STATES, OUTCOMES);
 	kappa 	=  	f_kappa(p, STATES, OUTCOMES, inc);
 	phi 	=  	f_phi(p, STATES, OUTCOMES, inc);
	sqm  	= 	MAX(0, 	((kappa*f_psqm(p, OUTCOMES)-phi)/(SCALE_PHIH2*p->phi_h2)));
return kappa * (inc - sqm*f_psqm(p, OUTCOMES) - f_mswcost(p, STATES, OUTCOMES));
}

double f_u(const param* p, declSTATES, declOUTCOMES){	// Current Utility
	// utility: returns utility for a given time,state point and decision	
	double inc, kappa, phi, sqm, utility;
	inc 	=  	f_get_inc(p, STATES, OUTCOMES);
 	kappa 	=  	f_kappa(p, STATES, OUTCOMES, inc);
	phi 	=  	f_phi(p, STATES, OUTCOMES, inc);
	sqm  	= 	MAX(0, 	((kappa*f_psqm(p, OUTCOMES)-phi)/(SCALE_PHIH2*p->phi_h2)));


	
	// Summed utility components.
	utility = 	phi*sqm + ((SCALE_PHIH2*p->phi_h2)/2)*sqm*sqm // utility of housing moving  
			+ 	f_regtaste(p, STATES, OUTCOMES) 	// regional taste diferences  
			- 	f_umvcost_rl(p, STATES, OUTCOMES) - f_umvcost_wl(p, STATES, OUTCOMES)	// utility cost of moving  
			-   f_ttimecost(p, OUTCOMES) 		// utility cost of commuting 
			+ 	kappa * (inc - sqm*f_psqm(p, OUTCOMES) - f_mswcost(p, STATES, OUTCOMES)); // utility value of income net of monetay cost
	
	// substract disutility of working
	if (irw != (p->nrwork-1)) 
		utility = utility - p->c_work;

	return utility;
}  

// JOB TRANSITION PROBAILITY
double f_tr_job(const param* p, declSTATES, const int dwl, const int job_transition, int *irwp1){
// dwl is CHOICE = intended job location
// irwp1 is OUTCOME
// job_transition can be 
// 1	: move to chosen job 
// 0	: stay in current job
// -1	: move to unemployment 
// f_job_tr must sum to 1 over job_transition
// pn=probability of succeeding in getting a new job. 
// pl=probability of keeping current job	
double job_tr;
	switch(job_transition){
		case 1	: // move to chosen job
	  		*irwp1=dwl;
			job_tr=f_pn(p, STATES, dwl); 
			break;
		case 0	: // stay in current job
	  		*irwp1=irwp;
			job_tr=(1-f_pn(p, STATES, dwl))*f_pk(p, STATES); 
			break;
		case -1	: // move to unemployment 
	  		*irwp1=(p->nrwork-1);
			job_tr=(1-f_pn(p, STATES, dwl))*(1-f_pk(p, STATES));
			break;
		default: 
			 job_tr=0;
			 printf("job_transition=%d\n",job_transition);
			 mexErrMsgTxt("Infeasible job-transtion: must be -1, 0 or 1"); 
			 break;
	}
	return job_tr;
}



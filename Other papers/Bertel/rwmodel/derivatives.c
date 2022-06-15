// The model specific derivatives restidental and location model

/*
PARAMS Not to be included 
beta_cons, beta_age, beta_age2, beta_unemp, b (never used since we estimate beta in 1st step)

PARAMS DONE (TESTED)
alpha_dining ,alpha_cafes ,alpha_cafes_n ,alpha_cafes_sqm ,alpha_dining_n ,alpha_dining_sqm ,
alpha_hotels ,alpha_hotels_n ,alpha_hotels_sqm ,
alpha_s0 ,alpha_s0_dens ,alpha_s1 ,alpha_s1_dens ,alpha_s2 ,alpha_s2_dens; 
kappa_0,kappa_y,kappa_y2
gamma_0,gamma_a,gamma_a2,gamma_ms,gamma_s,gamma_c
phi_0,phi_a,phi_a2,phi_ms,phi_y,phi_y2,phi_s,phi_c,phi_rh
eta_ttime,eta_ttime2

PARAMS IN UTILITY THAT STILL NEEDS DERIVS
psi_uc
gamma_sell,gamma_buy

PARAMS IN Transition THAT STILL NEEDS DERIVS
	// probability of succeeding in getting a new job. 
	mp.Bpn_0
	mp.Bpn_a
	mp.Bpn_s
	mp.Bpn_unemp

	% probability of succeeding in keeping current job. 
	mp.Bpk_0
	mp.Bpk_a
	mp.Bpk_s
	mp.Bpk_unemp


PARMS IN,BELLMAN THAT STILL NEEDS DERIVS
lambda_wl, lambda_rl 
delta
*/


double g_phi(const param* p, declSTATES, declOUTCOMES, const double inc, int pname, int j){
	// derivative of heterogeneous parameter, phi, affecting the baseline marginal utility of housing 
	double age, scale=1;

	if (p[0].reduced_form_sqm)
 		scale= -(SCALE_PHIH2*p->phi_h2);

	switch(pname){
		case IDX_phi_h2: 
			if (p[0].reduced_form_sqm)
				return SCALE_PHIH2*f_phi(p, STATES, OUTCOMES, inc)/(SCALE_PHIH2*p->phi_h2);			
			break;
		default:
			break;
	}
	return 0;
}

double g_kappa(const param* p, declSTATES, declOUTCOMES, const double inc, int pname, int j){
	// Derivative of heterogeneous parameter, phi, affecting the baseline marginal utility of housing 
	double age, scale=1;
	if (p[0].reduced_form_sqm)
 		scale=  (SCALE_PHIH2*p->phi_h2)/p[0].psi_uc; 

	switch(pname)
	{
	case IDX_phi_h2: 
		if (p[0].reduced_form_sqm)
			return SCALE_PHIH2*f_kappa(p, STATES, OUTCOMES, inc)/(SCALE_PHIH2*p->phi_h2);			
		break;
	default:
		break;
	}
	return 0;
}

double g_sqm(const param* p, declSTATES, declOUTCOMES, double inc, double sqm, int pname, int j){
	// Derivative of demand for housing in square meters given states. 

	double age, dphi, dkappa, dsqm;

	if (sqm<0.000000000000000000001) return 0.0; // housing demand truncated, derivative is zero

	dphi 	=	g_phi(p, STATES, OUTCOMES, inc, pname, j);
	dkappa 	=	g_kappa(p, STATES, OUTCOMES, inc, pname, j);
	dsqm=(f_psqm(p, OUTCOMES)*dkappa - dphi)/(SCALE_PHIH2*p->phi_h2); 
	if (pname==IDX_phi_h2)
		dsqm=dsqm-SCALE_PHIH2*sqm/(SCALE_PHIH2*p->phi_h2);
	return dsqm; 
}

double g_umvcost_rl(const param* p, declSTATES, declOUTCOMES, int pname, int j){	
	// Derivative of utility of moving costs, depending on couple status and kids. 
	
	double age;
	
	if (irhp == irh) // moving cost if staying
		return 0.0; 

	switch(pname){

		case IDX_gamma_0: 
			return 1; 
			break; 		
		
		case IDX_gamma_c: 
			if (j==(ikids-1)) return 1;
			break; 		
		
		case IDX_gamma_s: 
			if (j==(is-1)) return 1;		
			break; 		
		
		case IDX_gamma_ms: 
			if (icouple==1) return 1;		
			break; 		
		
		case IDX_gamma_a: 
			age=(double) (it + p[0].t0); 
			return age;
			break; 		
		
		case IDX_gamma_a2: 
			age=(double) (it + p[0].t0); 
			return age*age/1000;
			break; 

		default:
			return 0;
	        break;		
	}

	
	return 0;
}

double g_umvcost_wl(const param* p, declSTATES,  declOUTCOMES, int pname, int j){	
	
	double age;
	
	if (irwp == irw) // moving cost if staying
		return 0.0; 

	switch(pname){

		case IDX_rho_0: 
			return 1; 
			break; 		
		
		case IDX_rho_a: 
			age=(double) (it + p[0].t0); 
			return age;
			break; 		
	
		case IDX_rho_s: 
			if (j==(is-1)) return 1;		
			break; 		

		case IDX_rho_ttime_rwp: 
			if ((irw !=(p[0].nrwork-1)) & (irwp !=(p[0].nrwork-1)))
				return p[0].ttime[f_IDCHOICE(p, irw,irwp)]; 
			break; 		

		case IDX_rho_ttime_rhp: 
			if ((irw !=(p[0].nrwork-1)) & (irhp !=(p[0].nrwork-1))) 
				return p[0].ttime[f_IDCHOICE(p, irw,irhp)]; 
			break; 		

		case IDX_rho_unemp: 
			if (irwp ==(p[0].nrwork-1)) return 1; 
			break; 		
				
		default:
			return 0;
	        break;		
	}

	
	return 0;	
}

double g_regtaste(const param* p, declSTATES, declOUTCOMES, int pname, int j){	
	// Derivative of regional amenities
	// FIXME: Add derivatives of amenities	

	switch(pname){
		case IDX_alpha_dining:
	        return p[0].dining[irh];
	        break;
		case IDX_alpha_cafes:
			return p[0].cafes[irh];
			break;
		case IDX_alpha_cafes_n:
			return p[0].cafes_n[irh];
			break;
		case IDX_alpha_cafes_sqm:
			return p[0].cafes_sqm[irh];
			break;
		case IDX_alpha_dining_n:
			return p[0].dining_n[irh];
			break;
		case IDX_alpha_dining_sqm:
			return p[0].dining_sqm[irh];
			break;
		case IDX_alpha_hotels:
			return p[0].hotels[irh];
			break;
		case IDX_alpha_hotels_n:
			return p[0].hotels_n[irh];
			break;
		case IDX_alpha_hotels_sqm:
			return p[0].hotels_sqm[irh];
			break;
		case IDX_alpha_s0:
			return p[0].s0[irh];
			break;
		case IDX_alpha_s0_dens:
			return p[0].s0_dens[irh];
			break;
		case IDX_alpha_s1:
			return p[0].s1[irh];
			break;
		case IDX_alpha_s1_dens:
			return p[0].s1_dens[irh];
			break;
		case IDX_alpha_s2:
			return p[0].s2[irh];
			break;
		case IDX_alpha_s2_dens:
			return p[0].s2_dens[irh];
			break;
		
		case IDX_alpha_rh: 
			if (j==(irh-1)) return 1;
			break; 
		default:
			return 0;
	        break;		
	}

	return 0;
}

double g_ttimecost(const param* p, declSTATES, declOUTCOMES, int pname, int j){	
	// Derivative of utility travel time costs costs, depending on couple status and kids. 
	
	switch(pname){		
		case IDX_eta_ttime: 
			return p[0].ttime[IDCHOICE];
			break; 		
		case IDX_eta_ttime2: 
			return p[0].ttime[IDCHOICE]*p[0].ttime[IDCHOICE];
			break; 		
		default:
			return 0;
	        break;		
	}

	return 0;
}

double g_util(const param* p, declSTATES, declOUTCOMES, double inc, double sqm, int pname, int j){
	double phi,kappa, du,dphi, dsqm, dkappa, dregtaste, dttimecost, dumvcost;

	phi  		= 	f_phi(p, STATES, OUTCOMES, inc);
	kappa 		= 	f_kappa(p, STATES, OUTCOMES, inc);	

	dphi 		=	g_phi(p, STATES, OUTCOMES, inc, pname, j);
	dsqm		=	g_sqm(p, STATES, OUTCOMES, inc, sqm, pname, j);
	dkappa 		=	g_kappa(p, STATES, OUTCOMES, inc, pname, j);
	dregtaste	=	g_regtaste(p, STATES, OUTCOMES, pname, j);
	dumvcost 	=	g_umvcost_rl(p, STATES, OUTCOMES, pname, j)+g_umvcost_wl(p, STATES, OUTCOMES, pname, j);
	dttimecost 	=	g_ttimecost(p, STATES, OUTCOMES, pname, j);
	du 			=	dphi*sqm + (phi + (SCALE_PHIH2*p->phi_h2)*sqm - kappa*f_psqm(p, OUTCOMES))*dsqm
				+ dkappa * (inc - sqm*f_psqm(p, OUTCOMES) - f_mswcost(p, STATES, OUTCOMES))
				+ dregtaste - dumvcost - dttimecost;
	
	if (IDX_phi_h2==pname) 
		du=du+SCALE_PHIH2*sqm*sqm/2; 
	if (IDX_c_work==pname) // substract disutility of working
		if (irw != (p[0].nrwork-1)) du=du - 1;

	return du;
}


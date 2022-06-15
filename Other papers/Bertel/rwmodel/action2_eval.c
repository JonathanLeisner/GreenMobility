double f_logl_rl(const param *p, declSTATES, const int irh, const int irwp1, const double *evf, const double *devf,
                  double *v_rl1, double *P_rl, double *dP_rl)
{   // This function computes the probability of observing a particular transition 
    // from previous locations (irhp, irwp) to new (potentially realized) resideial location (irh)
    // CONDITIONAL on work locaton irwp1!!! 
    // IF irw == COMPUTEALL P_rl(irhp1) is computed for all residential locations irhp1 in (0, p.nrhome-1)
    // Outputs are in the last 3 arguments 
 
    //int STATES1;    // indices for states in the next period
    int irhp1;
    double lsum_rl;
    
    // for derivatives
    int ipar;
    double inc, sqm;
    double dlsum_rl,*dv_rl1;

    if (dP_rl!=NULL){
        dv_rl1 = (double*) calloc((p->n_param*p->nrhome) ,sizeof(double)); 
        for (ipar = 0; ipar < p->n_param; ++ipar)
            dP_rl[ipar]=0;
    }

    // Likelihood for end of period residential location (irhp1), conditional on work-location (irwp1)
    for (irhp1=0;irhp1<p->nrhome;irhp1++)
    {   //loop over residential location choices
        v_rl1[irhp1]= f_u(p, STATES, irhp1, irwp1) + p->delta*p->survival[SURVIDX]*evf[f_IDEVF(p,it, irhp1, irwp1, X)];
        
        if (dP_rl!=NULL){
            // derivatives of residential-choice specific values 
            inc  = f_get_inc(p, STATES, irhp1, irwp1);
            sqm  = f_sqm(p, STATES, irhp1, irwp1);

            for (ipar = 0; ipar < p->n_param; ++ipar){
                 dv_rl1[irhp1 + ipar*p->nrhome] 
                     = g_util(p, STATES, irhp1, irwp1, inc, sqm, p->iparkey[ipar], p->isubpar[ipar])
                         + p->delta*p->survival[SURVIDX]*devf[f_IDEVF(p,it, irhp1, irwp1, X) + ipar*p->nstate];
            }    
        }
    }
 
    lsum_rl=logsum(p->nrhome,v_rl1,p->lambda_rl); 

    // Compute all ccp's
    if ((P_rl!=NULL) && (v_rl1!=NULL) && ((dP_rl!=NULL)  || (irh==COMPUTEALL))) 
        for (irhp1=0;irhp1<p->nrhome;irhp1++)
            P_rl[irhp1]=exp((v_rl1[irhp1]-lsum_rl)/p->lambda_rl);

    //Now do the derivatives 
    if (dP_rl!=NULL)
    {

        for (ipar = 0; ipar < p->n_param; ++ipar){
            // derivative of log(P_rl)
            
            //  1.a compute derivative of lsum_rl: CCP weighted dvf
            dlsum_rl=0;
            for (irhp1=0;irhp1<p->nrhome;irhp1++){
                dlsum_rl+=P_rl[irhp1]*dv_rl1[irhp1+ ipar*p->nrhome];
            }
            dP_rl[ipar]=P_rl[irh]*(dv_rl1[irh+ ipar*p->nrhome]-dlsum_rl)/p->lambda_rl;

        }
    free(dv_rl1);
    }

    // return log(P_rl(irh|irw))
    if (irh==COMPUTEALL) 
        return -99;
    return (v_rl1[irh]-lsum_rl)/p->lambda_rl;
}

double f_logl_wl(const param *p, declSTATES, const int irw,  const double *evf, const double *devf,
                  double *P_wl, double *dP_wl)
{   // This function computes the probability of observing a particular job transition 
    // from previous job locations to new (potentially observed) locations
    // Outputs are in the last 2 arguments (scalar, scalar)
 
    int irwp1, dwl;        // indices for states in the next period
    int i_job_tr;   // indicator for job transition (-1,0,1)
    double lsum_wl, tr_job, ev_dwl[MAXREGIONS], ccp_wl[MAXREGIONS];
    
    // for derivatives
    int ipar;
    double *dlsum_wl, *dev_dwl;

    for (ipar = 0; ipar < p->n_param; ++ipar)
        dP_wl[ipar]=0;
    dlsum_wl = (double*) calloc(p->n_param,sizeof(double)); 
    dev_dwl = (double*) calloc(p->n_param*p->nrwork,sizeof(double));     
    
    //compute values from all choices (integrated over job search outcomes)
    lsum_wl=f_vf(p, STATES, evf, ev_dwl, devf, dev_dwl, dlsum_wl);

    for (dwl=0;dwl<p->nrwork;dwl++)
        ccp_wl[dwl]= exp((ev_dwl[dwl]-lsum_wl)/p->lambda_wl); 

    // Likelihood observed work-location (irw). 
    for (dwl=0;dwl<p->nrwork;dwl++) P_wl[dwl]=0; 

    for (i_job_tr =-1;i_job_tr<=1;i_job_tr++)
    {   //loop over job transitions
        for (dwl=0;dwl<p->nrwork;dwl++)
        {   //loop over work location choices
            // Next period work location   
            tr_job=f_tr_job(p, STATES, dwl, i_job_tr, &irwp1);          
            if ((irwp1==irw) || (irw==COMPUTEALL))
            {
                P_wl[irwp1]+=tr_job*ccp_wl[dwl]; 
                for (ipar = 0; ipar < p->n_param; ++ipar)
                    dP_wl[ipar]+=tr_job*ccp_wl[dwl]*
                        (dev_dwl[dwl + ipar*p->nrwork] - dlsum_wl[ipar])/p->lambda_wl;  
            } 
        }  // end of drl          
    }  
    free(dlsum_wl);
    free(dev_dwl);
    if (irw==COMPUTEALL) 
        return -99;
    return log(P_wl[irw]);
     
}  


/*********************************************************************************************************
 f_likehood calculates the probabilities of location outcomes and its derivatives
 *********************************************************************************************************
 
  SYNTAX             :   
 
 INPUTS
 ---------------------------------------------------------------------------------------------------------
    p               :   Parameter structure that contains model parameters, grids, quadrature points, etc. 
                        (see help in rwmodel.m)

    evf             :   Matrix with the solution to the model. Pre-calculated in solver.c. nstates x nchoices.  

    UPDATE!!!

**********************************************************************************************/

void   f_likehood(const param *p, declSTATES, declOUTCOMES, const double *evf, const double *devf,
                  double *P_wl, double *P_rl, double *dP_wl, double *dP_rl)
{   // This function computes the probability of observing a particular transition 
    // from previous locations to new (potentially observed) locations
    // Outputs are in the last 4 arguments (scalar, scalar, )
 
    double logl_rl, logl_wl, v_rl1[MAXREGIONS];    
    logl_wl=f_logl_wl(p, STATES,irw, evf, devf, P_wl, dP_wl);
    logl_rl=f_logl_rl(p, STATES,OUTCOMES, evf, devf, v_rl1, P_rl, dP_rl);
        
}  



/*******************************************************************************************************************
 f_logl computes log likelihood and its derivatives given sp[ip] and cp[ip].
 *******************************************************************************************************************

 SYNTAX             :   f_logl(&p, SP, CP, evf, devf, &logl[ip], sp_obs, &dlogl[ip])
 
 INPUTS
 ---------------------------------------------------------------------------------------------------------
    p               :   Pointer parameter structure that contains model parameters, grids, quadrature points, etc. 
                        (see help in rwmodel.m)

    evf             :   Matrix with the solution to the model. Pre-calculated in solver.c. nstates x 1.  

    devf            :   Matrix with derivatives of the solution to the model. Pre-calculated in solver.c. nstates x nparameters.  

    N               :   Number of observations in state data

    &logl           :   Pointer to log-likelihood contribution of ip

    &dlogl          :   Pointer to derivatibe of log-likelihood contribution of ip

**********************************************************************************************/

void   f_logl(param *p, declSTATES, declOUTCOMES, const double *evf, const double *devf, 
              double *logl, const int N, double *dlogl)
{   //This function computes the log likelihood and it's derivatives

    double logl_rl, logl_wl, P_wl[MAXREGIONS],P_rl[MAXREGIONS],v_rl1[MAXREGIONS];
    
    // for derivatives
    int ipar;
    double *dP_wl,*dP_rl;

    if (dlogl==NULL)
    {   //without derivatives
        p->n_param=0;
        logl_wl=f_logl_wl(p, STATES, irw, evf, NULL, P_wl, NULL);
        logl_rl=f_logl_rl(p, STATES,OUTCOMES, evf, NULL, v_rl1, NULL, NULL);
    }
    else
    {   //with derivatives
        dP_wl = (double*) calloc((p->n_param) ,sizeof(double)); 
        dP_rl = (double*) calloc((p->n_param) ,sizeof(double)); 

        logl_wl=f_logl_wl(p, STATES,irw, evf, devf, P_wl, dP_wl);
        logl_rl=f_logl_rl(p, STATES,OUTCOMES, evf, devf, v_rl1, P_rl, dP_rl);
        
        //Derivatives
        for (ipar = 0; ipar < p->n_param; ++ipar) //loop over parameters
            dlogl[ipar*N] = dP_wl[ipar]/P_wl[irw] + dP_rl[ipar]/P_rl[irh];

        free(dP_rl);
        free(dP_wl);
    }
    logl[0]=logl_wl+logl_rl;

}

/*******************************************************************************************************************
 ccp_eval_all computes outcome transition probabilities for given sp[ip] and cp[ip].
 *******************************************************************************************************************
 
 SYNTAX             :   ccp_eval_all(p, evf, sp, ccp)
 
 INPUTS
 ---------------------------------------------------------------------------------------------------------
    p               :   Parameter structure that contains model parameters, grids, etc. 

    evf             :   Matrix with the solution to the model. Pre-calculated in solver.c. nstates x nchoices.  

    sp              :   Vector of state points. 1 x p->nstatevar.

 OUTPUT
 ---------------------------------------------------------------------------------------------------------
    ccp             :   Pointer to CCP of ip   
**********************************************************************************************/

static void ccp_eval_all(const param* p, const double* evf, declSTATES,  double *ccp)
{
    int ichoice, irw, irh;
    double logl_wl, logl_rl, P_wl[MAXREGIONS], P_rl[MAXREGIONS], v_rl1[MAXREGIONS]; 

    logl_wl=f_logl_wl(p, STATES, COMPUTEALL, evf, NULL, P_wl, NULL);

    for (irw=0; irw < p[0].nrwork; irw++){
        // compute P_rl conditional on states and irw
        logl_rl=f_logl_rl(p, STATES, COMPUTEALL, irw, evf, NULL, v_rl1, P_rl, NULL);
        for (irh=0; irh < p[0].nrhome; irh++){
            ichoice=f_IDCHOICE(p, irh, irw); 
            ccp[ichoice]=P_rl[irh]*P_wl[irw]; 
            
        }
    }
}











//  The solver for the rwmodel with sequential choice of 1) jobsearch and 2) residence regions
//  Included in solver.c

double f_vf(const param *p, declSTATES, const double *evf, double *ev_dwl, 
        const double *devf, double *dev_dwl, double *dvf)
{   //This function computes the choice specific values 

    int STATES1, dwl, drl;    // indices for states and choices, see definitions.c 
    int i_job_tr;  // indexes of discrete states and decisions tomorrow
    double vf, ls_rl[MAXREGIONS], lsum_wl,ccp_rl,ccp_wl, v_drl[MAXREGIONS];
    double tr_job, ichoice;

    int ipar; 
    double inc, sqm, *dv_wl_rl;

    if (dev_dwl!=NULL)
        dv_wl_rl = (double*) calloc((p->n_param*p->nrhome*p->nrwork) ,sizeof(double)); 

    // initialize
    for (dwl=0;dwl<p->nrwork;dwl++) ev_dwl[dwl]=0; 
    for (ipar = 0; ipar < p->n_param; ++ipar){
        for (dwl=0;dwl<p->nrwork;dwl++){
            if (dev_dwl!=NULL) dev_dwl[dwl + ipar*p->nrwork]=0;
        }
    }

    for (irwp1=0;irwp1<p->nrwork;irwp1++)
    { // for each work location outcome, irwp1
        // compute residential-choice specific values condtional on irwp1
        for (drl=0;drl<p->nrhome;drl++)
            v_drl[drl]= f_u(p, STATES, drl, irwp1) + p->delta*p->survival[SURVIDX]*evf[f_IDEVF(p,it,drl,irwp1,X)];

        // Compute logsum over residential choices, drl conditional on irwp1
        ls_rl[irwp1]=logsum(p->nrhome,v_drl,p->lambda_rl); 
        
        
        // ... and their derivatives
        if (dev_dwl!=NULL)
        {
            for (drl=0;drl<p->nrhome;drl++) 
            {
                // compute derivatives of residential-choice specific values 
                inc  = f_get_inc(p, STATES, drl, irwp1);
                sqm  = f_sqm(p, STATES, drl, irwp1);
             
                for (ipar = 0; ipar < p->n_param; ++ipar)
                {
                    dv_wl_rl[drl + irwp1*(p->nrhome) + ipar*(p->nrwork)*(p->nrhome)]
                    =exp((v_drl[drl]-ls_rl[irwp1])/p[0].lambda_rl)
                    *(g_util(p, STATES, drl, irwp1, inc, sqm, p->iparkey[ipar], p->isubpar[ipar])
                        +  p->delta*p->survival[SURVIDX]*devf[f_IDEVF(p,it, drl, irwp1, X) + ipar*p->nstate]);
                }
            }
        }
    }  // end of irwp1     

    for (dwl=0;dwl<p->nrwork;dwl++)
        for (ipar = 0; ipar < p->n_param; ++ipar)
            dev_dwl[dwl + ipar*p->nrwork]=0;
    

    for (i_job_tr =-1;i_job_tr<=1;i_job_tr++){   //loop over job transitions
        for (dwl=0;dwl<p->nrwork;dwl++){   //loop over work location choices

            // Next period work location: stochastic choice dependent transition  
            // f_tr_job sets pointers to next period work location (irwp1) as a byproduct.   
            tr_job=f_tr_job(p, STATES, dwl, i_job_tr, &irwp1);

            // take expecation of over job outcome probability,  tr_job
            // to obtain expected value of work location choice, dwl  
            ev_dwl[dwl] += tr_job*ls_rl[irwp1]; // see eq 7 Maria/Christian thesis. 

            // ... and their derivatives 
            if (dev_dwl!=NULL)
            {     
                for (drl=0;drl<p->nrhome;drl++)
                {   // compute derivatives of expected residential-choice specific values
                    // conditional work location outcome                     
                    for (ipar = 0; ipar < p->n_param; ++ipar){
                        dev_dwl[dwl + ipar*p->nrwork]+=
                        tr_job*dv_wl_rl[drl + irwp1*(p->nrhome) + ipar*(p->nrwork)*(p->nrhome)];

                    }
                }
            }
        } // end of dwl
    } // end of i_job_tr

    vf=logsum(p->nrwork,ev_dwl,p->lambda_wl); // according to paper eq. 8 this is ev_w or ls_rw

    // finally their derivatives of vf
    if ((dvf!=NULL) && (dev_dwl!=NULL))
    {
        for (ipar = 0; ipar < p->n_param; ++ipar)
            dvf[ipar]=0;

        for (dwl=0;dwl<p->nrwork;dwl++)
        {
            ccp_wl=exp((ev_dwl[dwl]-vf)/p[0].lambda_wl);
            for (ipar = 0; ipar < p->n_param; ++ipar)
                dvf[ipar]+=ccp_wl*dev_dwl[dwl + ipar*p->nrwork];
        }
    }

    if (dev_dwl!=NULL)
            free(dv_wl_rl);

    return vf; 
}

void f_evf(const param* p, double *evf,double *devf,const int ievf) 
{   //This function computes the value function $v_{t+1}(wl,rl,x)$ and its derivative for all possible values of the state $(wl,rl,x)$ 
    //Only called for non-terminal periods
    //States are indexed by ievf 
   
    int STATES, X1;    // indices for states and choices, see definitions.c 
    double tr_x, evf_i, vf1, dv1[MAXNPARNAMES], ev_dwl[MAXREGIONS], *dev_dwl;
    int ipar;

    // int i_kids_tr, i_job_tr;  // indexes of discrete states and decisions tomorrow

    if (devf!=NULL)
        dev_dwl = (double*) calloc((p->n_param*p->nrwork) ,sizeof(double)); 
    else 
        dev_dwl = NULL;

    // STATE index
    // Set pointers of state variables to appropriate values given ievf (state space index, which dpsolver is looping over). 
    f_IDEVF_INV(p, ievf, pointSTATES);  // STATES now have current state values

    //Loop over next period states 
    evf_i=0.0; //initialize accumulating vars
    for (ipar = 0; ipar < p->n_param && devf!=NULL; ++ipar)
            devf[ievf+ ipar*p->nstate] = 0.0;
    
#ifdef FOOLPROOF
    // Check that transition probs sum up to one
    double trprsum = 0.0;
    // HERE: INCOMPATIBLE WITH OPENMP, search for 'pragma critical' if needs to work
#endif   

    // most of the states are time invariant (need to assign to fill out X1)
    im1 = im;
    is1 = is;
    iinc1 = iinc;
    imove1 = imove;
    // only kids and couple state can change (apart from locations/employment)
    for (ikids1=0;ikids1<p->nkids;ikids1++)
    {   //loop over next period kids state
        for (icouple1=0;icouple1<p->ncouple;icouple1++)
        {   //loop over next period couples state

            // Transition probability by table lookup
            tr_x = p->trpr[TRPRfromIDX + p->m_trpr * TRPRtoIDX];

#ifdef FOOLPROOF
    trprsum += tr_x;
#endif   

            // Skip zero probability paths
            if (tr_x<TOLERANCE) continue;

            // Next period value function
            // Note that irhp, irwp are post decision states and thus next period state is irhp, irwp
            vf1=f_vf(p, it+1, irhp, irwp, X1, evf, ev_dwl, devf, dev_dwl, dv1); 

            // accumulate expected value function with transition probability
            evf_i+=vf1*tr_x;

            //Now do the derivatives 
            for (ipar = 0; ipar < p->n_param && devf!=NULL; ++ipar){
                devf[ievf+ ipar*p->nstate] += tr_x*dv1[ipar];  
            }

        } // end of couples loop
    } // end of kids loop

#ifdef FOOLPROOF
    if (fabs(trprsum-1.0)>TOLERANCE) printf("PROBLEM: trpr does not sum up to one!\n");
#endif   

    // // i_kids_tr runs over pre-specified possible changes in number of kids to reduce computation. I.e. one less, the same, 1 more.
    // for (i_kids_tr=0;i_kids_tr<1;i_kids_tr++) // NO TRANSITIONS -- TEMP
    // {   //loop over next period states
    //     tr_x=f_tr_x(p, it, X, pointX1, i_kids_tr); // f_tr sets pointers in X1 as a byproduct.

    //     // Next period value function
    //     // Note that irhp, irwp are post decision states and thus next period state is irhp, irwp
    //     vf1=f_vf(p, it+1, irhp, irwp, X1, evf, ev_dwl, devf, dev_dwl, dv1); 
        
    //     // accumulate expected value function with transition probability
    //     evf_i+=vf1*tr_x;

    //     //Now do the derivatives 
    //     for (ipar = 0; ipar < p->n_param && devf!=NULL; ++ipar){
    //         devf[ievf+ ipar*p->nstate] += tr_x*dv1[ipar];  
    //     }

    // } // end of loop over id kids.

    evf[ievf]=evf_i;

    if (devf!=NULL) free(dev_dwl);
}

static void dpsolver(const param* p, double* evf, double* devf) 
{   // Calculates the expected value for every post-decision state
    // (all states except locations + locations after choice is made)
    int it, ievf, ipar;
    double ev; 
 
    for(it=p->T - p->t0;it>=0;it--)
    {   //time iteration

    #ifdef OPENMP
    #pragma omp parallel for default(shared) private(ievf) schedule(dynamic)
    #endif
        for (ievf=p->bl_t*it;ievf<p->bl_t*(it+1); ievf++) 
        { // loop over post-decision states
            if (it==p->T - p->t0){
                evf[ievf]=0;
                if (devf!=NULL){
                    for (ipar=0; ipar<p->n_param; ipar++) {
                        devf[ievf+ipar*p->nstate]=0.0;
                    }
                }
            } 
            else {
                f_evf(p,evf,devf,ievf); // compute expected value function
            }
        }
    }
}





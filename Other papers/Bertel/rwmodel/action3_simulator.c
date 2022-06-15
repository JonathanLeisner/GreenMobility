static int choice_simulator(const param* p,const double rnd,const double *ccp); 
bool survival_simulator(const param* p,declSTATES, const double rnd);
void state_transition_simulator(const param* p,declSTATES, declOUTCOMES, declpointSTATES1, const double rnd); 

static void simulator(const param* p, 
                    const double* evf, 
                    const statedata *sp0, 
                    double *randstream, 
                    const int iN, 
                    const int N, 
                    double *sim_id_out,
                    statedata *sp, choicedata *cp)

{

    int ip, ip1, t, ichoice, OUTCOMES, sim_id;
    int *istate, *iccp, *index;
    double *vfs, *array, sum_p;
    double lsum, vf;
    double *ccp;
    double rnd;
    bool surv;

    ccp= (double*) calloc(p->nchoice, sizeof(double));
   
    sim_id = iN;

    for (t = 0; t < (p->sim_T+1); ++t)
    {
        // Note that ip changes from being the individual counter - which it is outside simulator()'s domain' in simulator.c -
        // to being a counter of individual + time in simulation.
        ip=t+ iN*(p->sim_T+1);
        ip1=ip+1;  

        if (t==0) sp[ip]=sp0[iN]; 
        sp[ip].ok=1; // FIXME (call check function for sp an cp instead)

        // Record id of simulated individual and point of time in simulation (2nd column in output).
        sim_id_out[ip] = sim_id;
        sim_id_out[ip + N*(p->sim_T+1)] = t; // Point in simulation time.

        // compute all choice probabilities
        ccp_eval_all(p, evf, SP,  ccp);
        // simulate choice
        rnd=randstream[ip]; 
        ichoice=choice_simulator(p, rnd, ccp); 
        f_IDCHOICE_INV(p, ichoice, pointOUTCOMES);
        cp[ip].irh=irh;  // FIXME Manual indexing here. Need some function to do this
        cp[ip].irw=irw;


        if (t < (p->sim_T)){ // if not in the last period of simulation
            // simulate next period data. 
            rnd=randstream[ip + 1*(p->sim_T+1)]; // next shock is used for state simulator
            state_transition_simulator(p, SP, CP, pointSP1, rnd); 
            sp[ip1].it=sp[ip].it + 1;
            sp[ip1].irh=cp[ip].irh; 
            sp[ip1].irw=cp[ip].irw; 

            // simulate survival
            rnd=randstream[ip + 2*(p->sim_T+1)]; // next shock is used for survival simulation
            surv = survival_simulator(p,SP,rnd);

            // if next period age is above maximum age regenerate age index to 0, i.e. set age=t0;
            if (sp[ip1].it > (p->T-p->t0) || !surv)   {
                sp[ip1].it = 0;
                sim_id = sim_id + N; // Individual gets a new id when age is reset.
            }               
        }
    }
    free(ccp);
}

static int choice_simulator(const param* p,const double rnd,const double *ccp) 
{   //this function simulates choice id for one person in this person's sp_i (state point for individual i)
    //NOTE: infeasible choices should be guaranteed to have 0.0 prob in ccp, 
    // so no additional checks for feasibility are done here.
    int ip, ichoice;
    double sum_ccp;

    sum_ccp = 0.0; //accumulating sum for choice probs
    for (ichoice=0;ichoice<p->nchoice;ichoice++)
    {   //over all choices
        sum_ccp+=ccp[ichoice]; //see compute_all_ccp() for structure of ccp

        if (rnd<sum_ccp) 
        {   // return the choice
            return ichoice;
        } 
    } //ichoice

#ifdef FOOLPROOF    
    mexErrMsgTxt("Error in choice_simulator.");
#endif    
    return 0;
}

bool survival_simulator(const param* p,declSTATES, const double rnd) 
{   // Returns 1 if simulated survival
    // function call helps conver SP into STATES
    return rnd < p->survival[SURVIDX];
}

void state_transition_simulator(const param* p,declSTATES, declOUTCOMES, declpointSTATES1, const double rnd) 
{
    // This function sets pointers to next period states by transition probabilities. 
    
    double sum_tr=0.0;
    int ik1,ic1;

    // most of the states are time invariant (need to assign to fill out STATES1)
    *im1 = im;
    *is1 = is;
    *iinc1 = iinc;
    *imove1 = imove;
    // only kids and couple state can change (apart from locations/employment)
    for (ik1=0;ik1<p->nkids;ik1++)
    {   //loop over next period kids state
        for (ic1=0;ic1<p->ncouple;ic1++)
        {   //loop over next period couples state

            // Transition probability by table lookup
            // TRPRtoIDX = (p->ncouple*ik1+ ic1 )
            sum_tr += p->trpr[TRPRfromIDX + p->m_trpr * (p->ncouple*ik1+ ic1)];
            
            if (rnd<sum_tr)
            {   
                // Break out of function when pointers to next period states are correctly assigned  
                // by transition probabilities and rnd. 
                *ikids1 = ik1;  // set pointers in STATES1
                *icouple1 = ic1;
                return;
            } 

        } // end of couples loop
    } // end of kids loop

    // if didn't yet return
    mexErrMsgTxt("Error in state_transition_simulator.");
}


/**********************************************************************************************
 Function that evaluates variuous funtions and its derivatives at given state points
 **********************************************************************************************
 
 SYNTAX             :   [f, g]=evaluate_fct(p, sp, cp, fct_name)
 
 INPUTS
 ---------------------------------------------------------------------------------------------------------
    p               :   Parameter structure that contains model parameters, grids, 
                        (see help in rwmodel.m)

    fct_name:       :   String with function ot evaluate. 
                        fct_name can take values 'inc', 'sqm', 'pk', 'pn'...
                        If for example fct_name = 'f_util', then f=f_util(&p, SP, CP) and G returns its
                        derivatives g=g_util(&p, SP, CP, inc, p.iparkey[ipar], p.isubpar[ipar]);

    sp              :   Matrix with individuals' state points. N x p.nstatevar.

    cp:             :   Matrix with individuals' choice points. N x p.nchoicevar.

    yr:             :   Matrix-column with individuals' years points. N x 1.
 
 OUTPUTS
 ---------------------------------------------------------------------------------------------------------
    f           :   (Nx1) function evaluated at data 
    g           :   (N x p.n_param) derivative of f wrt. p.n_param parameters speified in p.param  
**********************************************************************************************/

//include other code files in the correct order
#include "definitions.c"
#include "foolproof.c"
#include "modelpart.c"
#include "derivatives.c"
#include "logit.c"
#include "action1_solver.c"
#include "action2_eval.c"

// declaration of functions
static mxArray* create_out_structure(const param* p, const mxArray *inputstate, const mxArray *inputchoice);

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////      Gate called from Matlab      ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   /* The gate to be called from MatLab */
    param p;
    bool compute_derivatives;
    int i, sp_obs, cp_obs, yr_obs, ip, ipar, compute_derivs, rl;
    int ip_y0, ip_yN, iy;
    double *out=NULL, *dout=NULL, inc, sqm;
    char *fct_name;
    statedata *sp;
    choicedata *cp;
    int *year_row_idx;
    double *evf, *devf, *logl, *dlogl, vf;//for likelihoods
    double tmp,buff_regions[MAXREGIONS],buff_regions2[MAXREGIONS]; //buffers for region specific values
    double buff_param[MAXNPARNAMES]; //buffer for single gradient
    double buff_regions_param[MAXNPARNAMES*MAXREGIONS]; //buffer for region specific derivatives

    /* Check the inputs and number of outputs */
    if (nrhs!=5) mexErrMsgTxt("Error: wrong number of inputs! Need 5 inputs");
    if (nlhs>2) mexErrMsgTxt("Error: wrong number of outputs! Need max 2 outputs");

#ifdef FOOLPROOF
    mexPrintf("FOOLPROOF checks are performed, this should be off in production runs\n");
#endif

#ifdef OPENMP
    #if defined (FOOLPROOF)
        printf("Evaluate_fct running OPENMP with %d cores out of %d available\n",omp_get_max_threads(),omp_get_num_procs());
        mexErrMsgTxt("Run OPENMP with SKIP_FOOLPROOF");
    #endif
    #ifdef NUM_CORES
        omp_set_num_threads(NUM_CORES);
    #else
        omp_set_num_threads(omp_get_num_procs());
    #endif
#endif

    //-------------------------------------------------------------------------------------------
    // INPUT: p: read in parameters common to all years by setting year to NULL
    //-------------------------------------------------------------------------------------------
    compute_derivatives = nlhs > 1;
    f_parseparameters(&p,prhs[0],compute_derivatives,INITIAL_PARSE);

    //-------------------------------------------------------------------------------------------
    // INPUT: fct_name
    //-------------------------------------------------------------------------------------------
    fct_name = (char*) mxArrayToString(prhs[1]);

   //-------------------------------------------------------------------------------------------
    // INPUT: sp and cp and year
    //-------------------------------------------------------------------------------------------

    // Check that there is equally many observations in input data on choices and states.
    sp_obs =  (int) mxGetM(prhs[2]); // Dimension of state point matrix
    cp_obs =  (int) mxGetM(prhs[3]); // Dimension of choice point vector
    yr_obs =  (int) mxGetM(prhs[4]); // Dimension of choice point vector

    if (sp_obs != cp_obs) mexErrMsgTxt("Need equally many observations in input data on choices and states.");
    if (sp_obs != yr_obs) mexErrMsgTxt("Need equally many observations in input data on states and years.");
    
    // Read in data on states to sp and cp
    sp = statedata_input(prhs[2]);       // statedata_input allocates memory to sp
    cp = choicedata_input(prhs[3]);     // choicedata_input allocates memory to cp
    year_row_idx = yeardata_input(&p,prhs[4]);  // choicedata_input allocates memory to cp

    //-------------------------------------------------------------------------------------------
    // PREPARE OUTPUT
    //-------------------------------------------------------------------------------------------
    if (nlhs>=1)
    {   // Column output for values
        if (strcmp(fct_name,"demand_h")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(p.nrhome, 1, mxREAL);
            out=(double*) mxGetPr(plhs[0]);
            for (i=0;i<p.nrhome;out[i++]=mxGetNaN()); //init with NaNs
        }
        else if (strcmp(fct_name,"ccp")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(sp_obs, p.nrhome, mxREAL);
            out=(double*) mxGetPr(plhs[0]);            
            for (i=0;i<p.nrhome;out[i++]=mxGetNaN()); //init with NaNs
        }
        else
        {
            plhs[0] = mxCreateDoubleMatrix(sp_obs, 1, mxREAL);
            out=(double*) mxGetPr(plhs[0]); 
            for (i=0;i<sp_obs;out[i++]=mxGetNaN()); //init with NaNs
        }
    }
    if (nlhs>=2)   
    {   // Compute derivatives of out put
        if (strcmp(fct_name,"demand_h")==0)
        {
            plhs[1] = mxCreateDoubleMatrix(p.nrhome, p.n_param, mxREAL);
            dout=(double*) mxGetPr(plhs[1]);
            // for (i=0;i<p.nrhome*p.n_param;dout[i++]=mxGetNaN()); //init with NaNs
        }
        else
        {
            plhs[1] = mxCreateDoubleMatrix(sp_obs, p.n_param, mxREAL);
            dout=(double*) mxGetPr(plhs[1]);
            // for (i=0;i<sp_obs*p.n_param;dout[i++]=mxGetNaN()); //init with NaNs
        }
    }
    else{ // skip computation of derivatives
        p.n_param=0;    
    }


    //-------------------------------------------------------------------------------------------
    // Loop over years in data 
    //-------------------------------------------------------------------------------------------

    for (iy = 0; iy < p.n_years_data; ++iy)
    {
        // First and last row of current year in data
        ip_y0 = year_row_idx[2*iy];
        ip_yN = year_row_idx[2*iy+1];

        if (ip_yN-ip_y0 == 0) continue; //skip the years with no data

        //-------------------------------------------------------------------------------------------
        // INPUT: p - parsed once for each year to account for time varying data and parameters
        //-------------------------------------------------------------------------------------------
    
        f_parseparameters(&p, prhs[0],compute_derivatives,iy);
        f_adjust_statedata(&p, sp, ip_y0, ip_yN);
        f_adjust_choicedata(&p, cp, ip_y0, ip_yN);

        // Precomputation of all possible incomes for year iy
        f_precomp_inc(&p);

        //-------------------------------------------------------------------------------------------
        // BRANCH OUT BY FUNCTION NAME
        //-------------------------------------------------------------------------------------------

        if (strcmp(fct_name,"f_inc")==0)
        {
            if (nlhs>=1)
            {   //first output: value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++) out[ip]=f_get_inc(&p, SP, CP);
            }
            if (nlhs>=2)
            {   //second output: derivative
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar)  schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++)
                {
                    inc     =   f_get_inc(&p, SP, CP);
                    for (ipar = 0; ipar < p.n_param; ++ipar)
                        dout[ip+ipar*sp_obs]=mxGetNaN(); // Not implemented, income params are estimated in first step.
                }        
            }
        }
        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"f_kappa")==0)
        {
            if (nlhs>=1)
            {   //first output: value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++)
                {
                    inc = f_get_inc(&p, SP, CP);
                    out[ip]=f_kappa(&p, SP, CP,inc);
                }        
            }
            if (nlhs>=2)
            {   //second output: derivative
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++)
                {
                    inc = f_get_inc(&p, SP, CP);
                    for (ipar = 0; ipar < p.n_param; ++ipar)
                        dout[ip+ipar*sp_obs]=g_kappa(&p, SP, CP, inc, p.iparkey[ipar], p.isubpar[ipar]);
                }        
            }
        }
        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"f_phi")==0)
        {
            if (nlhs>=1)
            {   //first output: value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++)
                {
                    inc = f_get_inc(&p, SP, CP);
                    out[ip]=f_phi(&p, SP, CP,inc);
                }        
            }
            if (nlhs>=2)
            {   //second output: derivative
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++)
                {
                    inc = f_get_inc(&p, SP, CP);
                    for (ipar = 0; ipar < p.n_param; ++ipar)
                        dout[ip+ipar*sp_obs]=g_phi(&p, SP, CP, inc, p.iparkey[ipar], p.isubpar[ipar]);
                }        
            }
        }
        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"f_sqm")==0)
        {
            if (nlhs>=1)
            {   //first output: value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++) out[ip]=f_sqm(&p, SP, CP);
            }
            if (nlhs>=2)
            {   //second output: derivative
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++)
                {
                    inc = f_get_inc(&p, SP, CP);
                    sqm = f_sqm(&p, SP, CP);
                    for (ipar = 0; ipar < p.n_param; ++ipar)
                        dout[ip+ipar*sp_obs]=g_sqm(&p, SP, CP, inc, sqm, p.iparkey[ipar], p.isubpar[ipar]);
                }        
            }
        }
        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"f_psqm")==0)
        {
            if (nlhs>=1)
            {   //first output: value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++) out[ip]=f_psqm(&p, CP);
            }
            if (nlhs>=2)
            {   //second output: derivative
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++)
                    for (ipar = 0; ipar < p.n_param; ++ipar)
                        dout[ip+ipar*sp_obs]=mxGetNaN(); //not yet implemented
            }
        }
        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"f_umvcost")==0)
        {
            if (nlhs>=1)
            {   //first output: value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++) out[ip]=f_umvcost_rl(&p, SP, CP)+f_umvcost_wl(&p, SP, CP);;
            }
            if (nlhs>=2)
            {   //second output: derivative
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++)
                    for (ipar = 0; ipar < p.n_param; ++ipar)
                        dout[ip+ipar*sp_obs]=   g_umvcost_rl(&p, SP, CP, p.iparkey[ipar], p.isubpar[ipar])
                                            +   g_umvcost_wl(&p, SP, CP, p.iparkey[ipar], p.isubpar[ipar]) ;
            }
        }
        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"f_mswcost")==0)
        {
            if (nlhs>=1)
            {   //first output: value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++) out[ip]=f_mswcost(&p, SP, CP);
            }
            if (nlhs>=2)
            {   //second output: derivative
                for(ip = ip_y0; ip < ip_yN; ip++)
                    for (ipar = 0; ipar < p.n_param; ++ipar)
                        dout[ip+ipar*sp_obs]=mxGetInf(); //not yet implemented
            }
        }
        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"f_ttimecost")==0)
        {
            if (nlhs>=1)
            {   //first output: value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++) out[ip]=f_ttimecost(&p, CP);
            }
            if (nlhs>=2)
            {   //second output: derivative
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++)
                    for (ipar = 0; ipar < p.n_param; ++ipar)
                        dout[ip+ipar*sp_obs]=g_ttimecost(&p, SP, CP, p.iparkey[ipar], p.isubpar[ipar]);
            }
        }
        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"f_regtaste")==0)
        {
            if (nlhs>=1)
            {   //first output: value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++) out[ip]=f_regtaste(&p, SP, CP);
            }
            if (nlhs>=2)
            {   //second output: derivative
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++)
                    for (ipar = 0; ipar < p.n_param; ++ipar)
                        dout[ip+ipar*sp_obs]=g_regtaste(&p, SP, CP, p.iparkey[ipar], p.isubpar[ipar]);
            }
        }
        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"f_u")==0 || strcmp(fct_name,"f_util")==0)
        {
            if (nlhs>=1)
            {   //first output: value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++) 
                    out[ip]=f_u(&p, SP, CP);
            }
            if (nlhs>=2)
            {   //second output: derivative
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++)
                {
                    inc = f_get_inc(&p, SP, CP);
                    sqm = f_sqm(&p, SP, CP);
                    for (ipar = 0; ipar < p.n_param; ++ipar)
                        dout[ip+ipar*sp_obs]=g_util(&p, SP, CP, inc, sqm, p.iparkey[ipar], p.isubpar[ipar]);
                }        
            }
        }
        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"f_u_sqm")==0)
        {
            if (nlhs>=1)
            {   //first output: value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++) 
                    out[ip]=f_u_sqm(&p, SP, CP);
            }
            if (nlhs>=2)
            {   //second output: derivative
                mexErrMsgTxt("Derivatives not implemented for f_u_sqm()!");
            }
        }
        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"f_u_money")==0)
        {
            if (nlhs>=1)
            {   //first output: value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++) 
                    out[ip]=f_u_money(&p, SP, CP);
            }
            if (nlhs>=2)
            {   //second output: derivative
                mexErrMsgTxt("Derivatives not implemented for f_u_money()!");
            }
        }
        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"f_pn")==0)
        {
            if (nlhs>=1)
            {   //first output: value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++) out[ip]=f_pn(&p, SP, cp[ip].irw);
            }
            if (nlhs>=2)
            {   //second output: derivative
                for(ip = ip_y0; ip < ip_yN; ip++)
                    for (ipar = 0; ipar < p.n_param; ++ipar)
                        dout[ip+ipar*sp_obs]=mxGetInf(); //not yet implemented
            }
        }
        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"f_pk")==0)
        {
            if (nlhs>=1)
            {   //first output: value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip, inc, sqm, ipar) schedule(dynamic)
                #endif
                for(ip = ip_y0; ip < ip_yN; ip++) out[ip]=f_pk(&p, SP);
            }
            if (nlhs>=2)
            {   //second output: derivative
                for(ip = ip_y0; ip < ip_yN; ip++)
                    for (ipar = 0; ipar < p.n_param; ++ipar)
                        dout[ip+ipar*sp_obs]=mxGetInf(); //not yet implemented
            }
        }
        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"expected_logsums_home_given_work_intention")==0)
        {   // logsum of the residence location choices, integrated over labor transitions,
            // thus conditional on work location intention (rather than outcomes)
            // computed at work outcomes in the data, as if they were observed intentions
            evf = (double*) calloc(p.nstate, sizeof(double));
            devf = (double*) calloc(p.n_param*p.nstate, sizeof(double));
            dpsolver(&p,evf,devf);
            if (nlhs>=1)
            {   //first output: value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip,inc,sqm,ipar,buff_regions) schedule(dynamic)
                #endif
                for (ip = ip_y0; ip < ip_yN; ++ip)
                {   //loop over all observations
                    f_vf(&p,SP,evf,buff_regions,devf,NULL,NULL);
                    out[ip] = buff_regions[sp[ip].irw];
                }
            }
            if (nlhs>=2)
            {   //second output: derivative
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip,inc,sqm,ipar,buff_regions,buff_regions_param) schedule(dynamic)
                #endif
                for (ip = ip_y0; ip < ip_yN; ++ip)
                {   //loop over all observations
                    f_vf(&p,SP,evf,buff_regions,devf,buff_regions_param,NULL);
                    out[ip] = buff_regions[sp[ip].irw];
                    for (ipar = 0; ipar < p.n_param; ++ipar) 
                        dout[ip+sp_obs*ipar]=buff_regions_param[sp[ip].irw + ipar*p.nrwork];
                }
            }
            free(evf);
            free(devf);
        }
        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"prob_observed_work_location")==0)
        {   // P_wl, the probability of observed work location outcome
            // computed at work outcomes in the data, part of the likelihood
            evf = (double*) calloc(p.nstate, sizeof(double));
            devf = (double*) calloc(p.n_param*p.nstate, sizeof(double));
            dpsolver(&p,evf,devf);
            if (nlhs==1)
            {   //first output: value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip,tmp,buff_regions) schedule(dynamic)
                #endif
                for (ip = ip_y0; ip < ip_yN; ++ip)
                {   //loop over all observations
                    tmp=f_logl_wl(&p, SP, cp[ip].irw, evf, NULL, buff_regions, NULL);
                    out[ip]=exp(tmp); 
                }
            }
            if (nlhs==2)
            {   //second output: derivative
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip,ipar, buff_regions, buff_param)  schedule(dynamic)
                #endif
                for (ip = ip_y0; ip < ip_yN; ++ip)
                {   //loop over all observations
                    tmp=f_logl_wl(&p, SP, cp[ip].irw, evf, devf, buff_regions, buff_param);
                    out[ip]=exp(tmp); 
                    for (ipar = 0; ipar < p.n_param; ++ipar) 
                         dout[ip+sp_obs*ipar] = buff_param[ipar];
                }
            }
            free(evf);
            free(devf);
        }
        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"prob_observed_home_location")==0)
        {   // P_wl, the probability of observed work location outcome
            // computed at work outcomes in the data, part of the likelihood
            evf = (double*) calloc(p.nstate, sizeof(double));
            devf = (double*) calloc(p.n_param*p.nstate, sizeof(double));
            dpsolver(&p,evf,devf);
            if (nlhs==1)
            {   //first output: value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip,inc,sqm,ipar,buff_regions, buff_regions2) schedule(dynamic)
                #endif
                for (ip = ip_y0; ip < ip_yN; ++ip)
                {   //loop over all observations
                    //logl_rl=f_logl_rl(p, STATES,OUTCOMES, evf, NULL, v_rl1, NULL, NULL);

                    tmp=f_logl_rl(&p, SP, CP, evf, NULL, buff_regions, NULL, NULL);
                    out[ip]=exp(tmp); 
                  }
            }
            if (nlhs==2)
            {   //second output: derivative
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip,ipar, buff_regions, buff_regions2,buff_param) schedule(dynamic)
                #endif
                for (ip = ip_y0; ip < ip_yN; ++ip)
                {   //loop over all observations
                    // f_likehood(&p,SP,CP,evf,devf,&tmp,&out[ip],NULL,buff_param);

                    tmp=f_logl_rl(&p, SP, CP, evf, devf, buff_regions, buff_regions2, buff_param);
                    out[ip]=exp(tmp); 
      
                    for (ipar = 0; ipar < p.n_param; ++ipar) 
                        dout[ip+sp_obs*ipar] = buff_param[ipar];
                }
            }
            free(evf);
            free(devf);
        }
        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"loglikelihood")==0 || strcmp(fct_name,"loglike")==0 || strcmp(fct_name,"logl")==0)
        {   // calculator of log-likelihood
            evf = (double*) calloc(p.nstate, sizeof(double));
            devf = (double*) calloc(p.n_param*p.nstate, sizeof(double));
            dpsolver(&p,evf,devf);
            if (nlhs==1)
            {   //only output value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip) schedule(dynamic)
                #endif
                for (ip = ip_y0; ip < ip_yN; ++ip)
                {   //loop over all observations
                    f_logl(&p, SP, CP, evf, devf, &out[ip], sp_obs, NULL);
                }
            }
            else if (nlhs==2)
            {   //value and derivatives
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip) schedule(dynamic)
                #endif
                for (ip = ip_y0; ip < ip_yN; ++ip)
                {   //loop over all observations
                    f_logl(&p, SP, CP, evf, devf, &out[ip], sp_obs, &dout[ip]);
                }
            }
            free(evf);
            free(devf);
        }
        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"demand_h")==0)
        {
            evf = (double*) calloc(p.nstate, sizeof(double));
            devf = (double*) calloc(p.n_param*p.nstate, sizeof(double));
            // dpsolver(&p,evf,devf);
            if (nlhs>=1)
            {   //first output: value
                // #ifdef OPENMP
                // #pragma omp parallel for default(shared) private(ip,inc,sqm,ipar,buff_regions, buff_regions2) schedule(dynamic)
                // #endif
                for (rl = 0; rl < p.nrwork; ++rl) out[rl]=0;
                // for (ip = ip_y0; ip < ip_yN; ++ip)
                // {   //loop over all observations
                //     tmp=f_logl_rl(&p, SP, COMPUTEALL, COMPUTEALL, evf, NULL, buff_regions, buff_regions2, NULL);
                //     for (rl = 0; rl < p.nrwork; ++rl) 
                //     {
                //         tmp=f_sqm(&p, SP, rl, cp[ip].irw);
                //         out[rl]+=tmp*buff_regions2[rl]; 
                //     }
                // }
        
            }
            if (nlhs>=2)
            {   //second output: derivative

            }
            free(evf);
            free(devf);
        }

        else if (strcmp(fct_name,"ccp")==0)
        {
            evf = (double*) calloc(p.nstate, sizeof(double));
            devf = (double*) calloc(p.n_param*p.nstate, sizeof(double));
            dpsolver(&p,evf,devf);
            if (nlhs>=1)
            {   //first output: value
                #ifdef OPENMP
                #pragma omp parallel for default(shared) private(ip,rl,buff_regions, buff_regions2) schedule(dynamic)
                #endif
                for (ip = ip_y0; ip < ip_yN; ++ip)
                {   //loop over all observations
                    tmp=f_logl_rl(&p, SP, COMPUTEALL, cp[ip].irw, evf, NULL, buff_regions, buff_regions2, NULL);
                    for (rl = 0; rl < p.nrhome; ++rl)
                    {
                        out[ip + rl*sp_obs]=buff_regions2[rl];
                    }
                }
            }
            if (nlhs>=2)
            {   //second output: derivative

            }
            free(evf);
            free(devf);
        }

        //-------------------------------------------------------------------------------------------
        else if (strcmp(fct_name,"***")==0)
        {
            if (nlhs>=1)
            {   //first output: value

            }
            if (nlhs>=2)
            {   //second output: derivative

            }
        }
        //-------------------------------------------------------------------------------------------
        else
            mexErrMsgTxt("Unknown function to run in evaluate_fct()!");

        //Free space allocated in parseparameters
        f_cleanparameters(&p,compute_derivatives,iy);

    } // End of computing on year iy 

    // Free data resources.
    free(sp);
    free(cp);
    free(year_row_idx);
    f_cleanparameters(&p,compute_derivatives,INITIAL_PARSE);
}


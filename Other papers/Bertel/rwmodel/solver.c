/**********************************************************************************************
 Solver for residential and work-location model
 By: Maria Juul-Hansen, Fedor Iskhakov, Christian Langholz, John Rust and Bertel Schjerning
 **********************************************************************************************
 
 SYNTAX             :   [ev, dev]=solver(p);
 
 INPUTS
 ---------------------------------------------------------------------------------------------------------
    p               :   Parameter structure that contains model parameters, grids, quadrature points, etc. 
                        (see help in rwmodel.m)
 
 
 OUTPUTS
 ---------------------------------------------------------------------------------------------------------
    Cell arrays by years in mp.year_data for:

    ev              :   r= by 1 matrix that holds the expected value function, where r=p.nstate 
                        is the number of discrete states. 
                        
                        Rows are indexed by ievf = 0,...p.nstate-1 and ordered by state combinations 
                        defined in f_IDEVF_INV and f_IDEVF (see definitions.c)  
                        

    dev (optional):     rxc array of derivatives of ev with respect to parameters in p.pnames 

                        r is the number of discrete states

                        c=p[0].n_param = number of parameters to differentiate EV wrt. 
                        (= total number of parameters implied by p.pnames and the size of each parameter vector). 
                        The order of the derivatives follows the ordering specified in p.pnames.
 
    Notes: ev fully characterizes the solution of the model
                         
**********************************************************************************************/

//include other code files in the correct order
#include "definitions.c"
#include "foolproof.c"
#include "modelpart.c"
#include "derivatives.c"
#include "logit.c"
#include "action1_solver.c"

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////      Gate called from Matlab      ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   /* The gate to be called from MatLab */
    param p;
    double *ev, *dev; 
    int i, iy, year;
    mxArray *tmp;
    bool compute_derivatives;

    /* Check the inputs and number of outputs */
    if (nrhs!=1) mexErrMsgTxt("Error: wrong number of inputs! Need 1 input: mp structure");
    if (nlhs>2) mexErrMsgTxt("Error: wrong number of outputs! Maximum 2 output is allowed");

#ifdef DEBUGMSG
    mexPrintf("\nSolver is run with DEBUGMSG and will print A LOT MORE diagnostics messages.\n");
#endif

#ifdef FOOLPROOF
    if (!mxIsStruct(prhs[0])) mexErrMsgTxt("Error: expecting structure of parameters for first argument");
    mexPrintf("FOOLPROOF checks are performed, this should be off in production runs\n");
#endif

#ifdef OPENMP
    #if defined (FOOLPROOF)
        mexErrMsgTxt("Run OPENMP with SKIP_FOOLPROOF and without DEBUGMSG");
    #endif
    #ifdef NUM_CORES
        omp_set_num_threads(NUM_CORES);
    #else
        omp_set_num_threads(omp_get_num_procs());
    #endif
    // printf("Solver running OPENMP with %d cores out of %d available\n",omp_get_max_threads(),omp_get_num_procs());
#endif

    // Parse parameters (function in defitions). Note: f_parse.. allocates memory for kids transitions, so free
    // Compute derivatives if more than 1 output argument 
    compute_derivatives = nlhs > 1;
    f_parseparameters(&p,prhs[0],compute_derivatives,INITIAL_PARSE);

#ifdef FOOLPROOF
    // Check if indices that retrieve points in the state and (state,choice) space are consistently calculated 
    // and inverted back to the states and choices they were based on. (function in foolproof)
    check_index(&p);
#endif
    
    //-------------------------------------------------------------------------------------------
    // OUTPUT 1: ev, expected value function (cell array by year)
    //-------------------------------------------------------------------------------------------

    // cell array corresponding to years
    plhs[0]=mxCreateCellMatrix(1,p.n_years_data);

    //-------------------------------------------------------------------------------------------
    // OUTPUT 2: dev, derivative of expected value function (optional)
    //-------------------------------------------------------------------------------------------
    if (compute_derivatives) plhs[1]=mxCreateCellMatrix(1,p.n_years_data);
    else dev=NULL; // do not compute gradients

    // Loop over years in mp.years_data 
    for (iy = 0; iy < p.n_years_data; iy++)
    {
        // Output space in the cell
        tmp = mxCreateDoubleMatrix(p.nstate,1,mxREAL);
        mxSetCell(plhs[0],iy,tmp);
        ev=(double*) mxGetPr(tmp);

        if (compute_derivatives)
        {
            if (p.n_pnames==0)
            {
                tmp=mxCreateDoubleMatrix(1,1,mxREAL);
                mxSetCell(plhs[1],iy,tmp);
                dev=(double*) mxGetPr(tmp);
                dev=NULL;
                if (iy==0) mexWarnMsgTxt("p.pnames is empty, cannot compute derivatives");
            }
            else
            {
                tmp = mxCreateDoubleMatrix(p.nstate,p.n_param,mxREAL);
                mxSetCell(plhs[1],iy,tmp);
                dev=(double*) mxGetPr(tmp);
            }
        }

        // Fetch year specific parameters
        f_parseparameters(&p,prhs[0],compute_derivatives,iy);

        // Precomputation of all possible incomes for year iy
        f_precomp_inc(&p);

#ifdef FOOLPROOF
        // Initialize every entry in ev to inf to visualize if entries are skipped during computation.
        for (i=0;i<p.nstate;i++) ev[i]=-mxGetInf();
#endif

        //call the solver
        dpsolver(&p,ev,dev);

        //Free space allocated in parseparameters per year
        f_cleanparameters(&p,compute_derivatives,iy);

    } // loop over years

    //Free space allocated in parseparameters independent of years
    f_cleanparameters(&p,compute_derivatives,INITIAL_PARSE);

}


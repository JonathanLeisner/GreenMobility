 /*******************************************************
 * simulate.c: Simulates a panel data set from the model *
 *******************************************************
 
 SYNTAX:            [sp, cp, yp]=simulator(p, evf, sp0, y0, randvec);  

 INPUTS:
  ---------------------------------------------------------------------------------------------------------
    p               :   Parameters. Same format as above (see definitions)
    sp0             :   N x p.nstates matrix of initial state points. Same format as sp under action_2_eval
                        (see STATES in definitions.c for order of states)
                        First argument is always age in {t0, T-1} not it in {0, T-t0-1}
                        [it,  irhp,  irwp,  ikids,  icouple,  im,  is,  iinc,  imove]
    yp0              :  N x 1 matrix of year in initial state points 
    randvec         :   uniform random numbers one for each simulated random event (currently T*N*3). FIXME     
 
 OUTPUTS:
 ---------------------------------------------------------------------------------------------------------
    sp            :   ((p.sim_T+1)*nsobs) x nstatevar: data set with columns described in p.s_names 
                        See actio3_sim.c for more help
    cp            :   ((p.sim_T+1)*nsobs) x nchoicevar: data set with columns described in p.c_names 
                       
    yp:           :   ((p.sim_T+1)*nsobs) x 1:  year associated wih each simulation

    Notes: When household gets older than the maximum age of households, T, the household dies and
    is replaced by a new household of age t0. The new household inherits the location of the old household. 
 
**********************************************************************************************/

#define RANDPERSONYEAR 3
// random numbers per 1 person-year simulation:
//    1 = simulate job/residence location decisions
//    2 = transition of children and couple state
//    3 = survival

#define RANDYEAR 0
// We also need per year for macro variables (such as common macro shock)
//    1 = No shocks yet 

//include other code files in the correct order
#include "definitions.c"
// #include "foolproof.c"
#include "modelpart.c"
#include "derivatives.c"
#include "logit.c"

#include "action1_solver.c"
#include "action2_eval.c"
#include "action3_simulator.c"


////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////      Gate called from Matlab      ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   /* The gate to be called from MatLab */
    param p;
    // int nscol, nccol, ntsp, nsobs, ncobs, nstate, nchoice;
    int nscol, nccol, nsobs, nrandstream;
    mxArray *evcell;
    // int *states, *statepoints, *choicepoints;
    // mxArray* outstruct;
    // double *evf, *vf, *lsum, *ccp;
    double *evf, *randstream;
    double *sim_states_out, *sim_choice_out, *sim_id_out; 
    statedata *sp0, *sp;
    choicedata *cp;
    int ip, i;
    int *year_row_idx;
    int ip_y0, ip_yN, iy;



    /* Check the inputs and number of outputs */
    if (nrhs!=5) mexErrMsgTxt("Error: wrong number of inputs! Need 5 inputs");
    if (nlhs!=3) mexErrMsgTxt("Error: wrong number of outputs! Need 3 outputs");

#ifdef FOOLPROOF
    if (!mxIsStruct(prhs[0])) mexErrMsgTxt("Error: expecting structure of parameters for first argument");
    mexPrintf("FOOLPROOF checks are performed, this should be off in production runs\n");
#endif

#ifdef OPENMP
    #if defined (FOOLPROOF)
        printf("Simulator running OPENMP with %d cores out of %d available\n",omp_get_max_threads(),omp_get_num_procs());
        mexErrMsgTxt("Run OPENMP with SKIP_FOOLPROOF");
    #endif
    #ifdef NUM_CORES
        omp_set_num_threads(NUM_CORES);
    #else
        omp_set_num_threads(omp_get_num_procs());
    #endif
#endif



    //**********************************
    // SECTION 1: Read inputs         *
    //**********************************

    // Input 1: Parse parameters ( function in definitions).
    // f_parseparameters(&p,prhs[0], false); 
    f_parseparameters(&p,prhs[0],false,INITIAL_PARSE);

    // do not compute derivatives when simulating model
    p.n_param=0;
    if (p.sim_T>1) mexWarnMsgTxt("BEWARE, the simulator is inconsistent with sim_T>1!");

    // Input 2: Read in data on model solution, evf
    if (!prhs[1]) mexErrMsgTxt("Input error in solution, evf");

    // Input 3: Read in data on states to sp0 for use in simulator. 
    if (!prhs[2]) mexErrMsgTxt("Input error in initial input states, sp0");
    nsobs =  (int) mxGetM(prhs[2]);
    nscol =  (int) mxGetN(prhs[2]);
    if (nscol != p.nstatevar) mexErrMsgTxt("Number of columns in sp0 do not match number of state variables in p");
    sp0 = statedata_input(prhs[2]);  

    // Input 4: 
    year_row_idx = yeardata_input(&p,prhs[3]);

    // Input 5: Read in data on rand stream 
    if (!prhs[4]) mexErrMsgTxt("Input error in random stream in the simulator");
    nrandstream  =  (int) mxGetNumberOfElements(prhs[4]);
    randstream= (double*)  mxGetPr(prhs[4]);


    if (nrandstream < (RANDPERSONYEAR*nsobs+RANDYEAR)*(p.sim_T+1)) 
        mexErrMsgTxt("Not enough random numbers in randstream, need (RANDPERSONYEAR*N+RANDYEAR)*(T+1)");

    //**********************************
    //* SECTION 2: Prepare outputs     *
    //**********************************

    // Allocate memory
    plhs[0] = mxCreateDoubleMatrix((p.sim_T+1)*nsobs, p.nstatevar, mxREAL); 
    sim_states_out = (double*)  mxGetPr(plhs[0]); 
    plhs[1] = mxCreateDoubleMatrix((p.sim_T+1)*nsobs, p.nchoicevar, mxREAL); 
    sim_choice_out = (double*)  mxGetPr(plhs[1]);     
    plhs[2] = mxCreateDoubleMatrix((p.sim_T+1)*nsobs, 2, mxREAL); // cols: pers_id, time_id
    sim_id_out = (double*)  mxGetPr(plhs[2]); 

    for (i = 0; i < (nsobs*(p.sim_T+1))*p.nstatevar; ++i) // loop over all persons and timer periods
        sim_states_out[i]=mxGetNaN();

    // Simulated state data
    sp= (statedata*)  calloc(nsobs*(p.sim_T+1)*p.nstatevar,sizeof(statedata));    // Space for the array of statedata structs 
    cp= (choicedata*) calloc(nsobs*(p.sim_T+1)*p.nchoicevar,sizeof(choicedata));    // Space for the array of statedata structs 
   
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
    
        evcell = mxGetCell(prhs[1], iy); //pointer to the solution for given year 
        evf  = mxGetPr(evcell);          //data from the solution

        f_parseparameters(&p, prhs[0],false,iy);
        f_adjust_statedata(&p, sp0, ip_y0, ip_yN);

        // // Precomputation of all possible incomes for year iy
        f_precomp_inc(&p);

#ifdef OPENMP
#pragma omp parallel for default(shared) private(ip) schedule(dynamic)
#endif
        for(ip = ip_y0; ip < ip_yN; ip++) // loop over all persons in sp0
        {
            simulator(&p, evf, sp0, randstream, ip, nsobs, sim_id_out, sp, cp);
        }
        
        //Free space allocated in parseparameters
        f_cleanparameters(&p,false,iy);

    } // End of computing on year iy 

    //copy over to the output space
    for (i = 0; i < nsobs*(p.sim_T+1); i++) // loop over all persons and timer periods
    {
        save_statedata(&p, sp, i, (p.sim_T+1)*nsobs, sim_states_out);
        save_choicedata(&p, cp, i,(p.sim_T+1)*nsobs, sim_choice_out);
    }

    // Free resources.
    free(sp);
    free(cp);
    free(year_row_idx);

    f_cleanparameters(&p,false,INITIAL_PARSE);


}



//  The definitions and some common functions for carmodel
//  Is included in carmodel.c

#include <string.h>
#include <stdlib.h>
#include "stdio.h"
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include <float.h>
#include <errno.h>


//OpenMP will be added when compiled with OpenMP support (-lopenmp)
#ifdef OPENMP
    #ifdef ONPC
        #include <omp.h>
    #else
        #include "omp.h"
    #endif
#endif

//checks that should not be performed in the fastest version
#ifndef SKIP_FOOLPROOF
 #define FOOLPROOF   //run checks everywhere
#endif

#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define SIGN(X) ((X)>0?1:-1)
#define DBL2INT(X) ((int)((X) + 0.05*SIGN((X))))
#define TOLERANCE 1e-10
#define SCALE_PHIH2 1e-4
//constant to signify infeasible choice (in value functions) MUST BE LARGE NEGATIVE
#define INFEASIBLECHOICEVF -998
//should be above the max length of vector parameter
#define MAXREGIONS 20

// set irw=COMPUTEALL or irh=COMPUTEALL when computing ccps for all alternatives
// COMPUTEALL cannot be between 0 and MAXREGIONS
#define COMPUTEALL -99

// order of states in all in and output matrices
#define ARGit       0  
#define ARGirhp     1
#define ARGirwp     2
#define ARGikids    3
#define ARGicouple  4
#define ARGim       5
#define ARGis       6
#define ARGiinc     7
#define ARGimove    8

// order of arguments in all in and output matrices
#define ARGirh      0
#define ARGirw      1

#define STATES      it,  irhp,  irwp,  ikids,  icouple,  im,  is,  iinc,  imove 
#define OUTCOMES    irh,   irw
#define STATES1     irhp1, irwp1, ikids1, icouple1, im1, is1, iinc1, imove1
// all states except post-decision states
#define X           ikids, icouple, im, is, iinc, imove 
#define X1          ikids1, icouple1, im1, is1, iinc1, imove1 

#define pointSTATES  &it, &irhp,  &irwp,  &ikids,  &icouple,  &im,  &is,  &iinc,  &imove 
#define pointSTATES1      &irhp1, &irwp1, &ikids1, &icouple1, &im1, &is1, &iinc1, &imove1 
#define pointX       &ikids,  &icouple,  &im,  &is,  &iinc,  &imove 
#define pointX1      &ikids1, &icouple1, &im1, &is1, &iinc1, &imove1 


#define pointOUTCOMES &irh, &irw

#define declX           const int ikids, const int icouple, const int im, const int is, const int iinc, const int imove 
#define declpointX      int *ikids,  int *icouple,  int *im,  int *is,  int *iinc,  int *imove
#define declpointX1     int *ikids1, int *icouple1, int *im1, int *is1, int *iinc1, int *imove1

#define declSTATES const int it, int irhp, int irwp, const int ikids, const int icouple, const int im, const int is, const int iinc, const int imove 
#define declpointSTATES int *it, int *irhp,  int *irwp,  int *ikids,  int *icouple,  int *im,  int *is,  int *iinc,  int *imove
#define declpointSTATES1         int *irhp1, int *irwp1, int *ikids1, int *icouple1, int *im1, int *is1, int *iinc1, int *imove1

#define declOUTCOMES  int irh, int irw
#define declpointOUTCOMES int *irh, int *irw

#define SP          sp[ip].it,  sp[ip].irh, sp[ip].irw, sp[ip].ikids, sp[ip].icouple, sp[ip].im, sp[ip].is, sp[ip].iinc, sp[ip].imove
#define pointSP1               &sp[ip1].irh, &sp[ip1].irw, &sp[ip1].ikids, &sp[ip1].icouple, &sp[ip1].im, &sp[ip1].is, &sp[ip1].iinc, &sp[ip1].imove
#define SP1                     sp[ip1].irh, sp[ip1].irw, sp[ip1].ikids, sp[ip1].icouple, sp[ip1].im, sp[ip1].is, sp[ip1].iinc, sp[ip1].imove

#define CP cp[ip].irh, cp[ip].irw

// Keys for parameters (used for computing derivatives)
// Ordering of keys is irrelevant, but you cannot have duplicates. 

#define IDX_gamma_0                 9
#define IDX_gamma_a                 10
#define IDX_gamma_a2                11
#define IDX_gamma_ms                12
#define IDX_gamma_s                 13
#define IDX_gamma_c                 14

#define IDX_eta_ttime               24
#define IDX_eta_ttime2              25

#define IDX_phi_h2                  33
#define IDX_rho_0                   34
#define IDX_rho_a                   35
#define IDX_rho_s                   36
#define IDX_rho_ttime_rhp           37
#define IDX_rho_ttime_rwp           38
#define IDX_rho_unemp               39

#define IDX_c_work                  40
#define IDX_alpha_rh                41

#define IDX_alpha_dining            42
#define IDX_alpha_cafes             43
#define IDX_alpha_cafes_n           44
#define IDX_alpha_cafes_sqm         45
#define IDX_alpha_dining_n          46
#define IDX_alpha_dining_sqm        47
#define IDX_alpha_hotels            48
#define IDX_alpha_hotels_n          49
#define IDX_alpha_hotels_sqm        50
#define IDX_alpha_s0                51
#define IDX_alpha_s0_dens           52
#define IDX_alpha_s1                53
#define IDX_alpha_s1_dens           54
#define IDX_alpha_s2                55
#define IDX_alpha_s2_dens           56
// Parameters above this point have certified working derivatives

#define MAXNPARNAMES                100
#define BADKEY                      -1

// key for 
#define INITIAL_PARSE               -99


int param_key(char *pname){
    // this function should contain a list of all the parameters we compute derivatives wrt. 

    if (strcmp(pname, "gamma_0") ==0)               return  IDX_gamma_0;       
    if (strcmp(pname, "gamma_a") ==0)               return  IDX_gamma_a;        
    if (strcmp(pname, "gamma_a2") ==0)              return  IDX_gamma_a2;       
    if (strcmp(pname, "gamma_ms") ==0)              return  IDX_gamma_ms;       
    if (strcmp(pname, "gamma_s") ==0)               return  IDX_gamma_s;        
    if (strcmp(pname, "gamma_c") ==0)               return  IDX_gamma_c;  

    if (strcmp(pname, "rho_0") ==0)                 return  IDX_rho_0;       
    if (strcmp(pname, "rho_a") ==0)                 return  IDX_rho_a;       
    if (strcmp(pname, "rho_s") ==0)                 return  IDX_rho_s;       
    if (strcmp(pname, "rho_unemp") ==0)             return  IDX_rho_unemp;       
    if (strcmp(pname, "rho_ttime_rhp") ==0)         return  IDX_rho_ttime_rhp;       
    if (strcmp(pname, "rho_ttime_rwp") ==0)         return  IDX_rho_ttime_rwp;       

    if (strcmp(pname, "c_work") ==0)                return  IDX_c_work;       
   
    if (strcmp(pname, "eta_ttime") ==0)             return  IDX_eta_ttime;      
    if (strcmp(pname, "eta_ttime2") ==0)            return  IDX_eta_ttime2; 

    if (strcmp(pname, "alpha_rh") ==0)              return  IDX_alpha_rh;

    if (strcmp(pname, "alpha_dining") == 0)         return  IDX_alpha_dining;
    if (strcmp(pname, "alpha_cafes") == 0)          return  IDX_alpha_cafes ;
    if (strcmp(pname, "alpha_cafes_n") == 0)        return  IDX_alpha_cafes_n;
    if (strcmp(pname, "alpha_cafes_sqm") == 0)      return  IDX_alpha_cafes_sqm;
    if (strcmp(pname, "alpha_dining_n") == 0)       return  IDX_alpha_dining_n ;
    if (strcmp(pname, "alpha_dining_sqm") == 0)     return  IDX_alpha_dining_sqm;
    if (strcmp(pname, "alpha_hotels") == 0)         return  IDX_alpha_hotels;
    if (strcmp(pname, "alpha_hotels_n") == 0)       return  IDX_alpha_hotels_n ;
    if (strcmp(pname, "alpha_hotels_sqm") == 0)     return  IDX_alpha_hotels_sqm;
    if (strcmp(pname, "alpha_s0") == 0)             return  IDX_alpha_s0 ;
    if (strcmp(pname, "alpha_s0_dens") == 0)        return  IDX_alpha_s0_dens;
    if (strcmp(pname, "alpha_s1") == 0)             return  IDX_alpha_s1 ;
    if (strcmp(pname, "alpha_s1_dens") == 0)        return  IDX_alpha_s1_dens;
    if (strcmp(pname, "alpha_s2") == 0)             return  IDX_alpha_s2 ;
    if (strcmp(pname, "alpha_s2_dens") == 0)        return  IDX_alpha_s2_dens;
 
    if (strcmp(pname, "phi_h2") ==0)                return  IDX_phi_h2;
    
    printf("Derivative for %s is not implemented", pname);
    mexErrMsgTxt("Need to add parameter to list of derivatives in definitions.c");
    return BADKEY;
}


typedef struct paramstruct{
    // Reduced form parameters indexing driving

// ############### parameters indexing dimensions ###############

    int T; // Max possible age. At age T, the decision problem ends and the current state is the abosorbing state. Consider using T_i instead: 65-data.age
    int t0; // lowest possible age
    int lifespan; //number of time periods
    int sim_T; // number of calendar years to simulate
    int y0; // initial calendar year
    int nmovec; // number of individual moving cost types
    int ninc; // number of unobserved income types
    int ns; // number of schooling levels
    int nm; // number of macro states
    int nkids; // max number of kids
    int ncouple; // number of couple states
            
    int nrwork;  // Number of work regions + unemployment choice
    int nrhome; // number of home regions (no unemployment choice)
    int nchoice; // number of coice combos
    int nstate; // number of states
    int nstatevar; // number of state variables: age, rhome_prev, rwork_prev ...
    int nchoicevar; // number of choice variables: rh, rw

    double dkk_scale; // scale of monetary units. E.g. 100,000s. 

    char *pnames[MAXNPARNAMES];
    int n_pnames, ni_param[MAXNPARNAMES], n_param, *iparname, *isubpar, *iparkey;

    int *years_data; // Years to solve the model in, should be in the data to make it possible
    int n_years_data, year, iy;
                                  
 // ############### model parameters ##################

    // travel time parameters
    double eta_ttime;          // time +ones(mp.Rhome*mp.Rwork)*0.05; +[0;0.4,0;0.8;0,0]; %ones(mp.Rhome,mp.Rwork);%transport cost parameter. For all (rh,rw)
    double eta_ttime2;       // time^2
    double *ttime; // travel time data. FIXME: is it correct that it is in the param vector (mp in rwmodel.m)?
    double delta;     // discount factor
    double *mymatrix; 
 
    // parameters of the income equations
    double *beta_cons;
    double *beta_age_pol;
    double *beta_age_fe;
    double *beta_rw;
    double *beta_l_unemp;
    double *beta_l_unemp_rw;
    double *b_0;
    double *b_age_pol;
    double *b_age_fe;
    int beta_min_age_fe, b_min_age_fe;
    // corresponding dimensions (when used outside of this file)
    int n_beta_age_pol, n_beta_age_fe, n_b_age_pol, n_b_age_fe;

    // parameters, job transitions
    double Bpn_0, Bpn_a, Bpn_unemp, *Bpn_s, Bpn_jobdensity, Bpn_ttime_rwp, Bpn_ttime_rhp;
    double Bpk_0, Bpk_a, Bpk_unemp, *Bpk_s, Bpk_jobdensity, Bpk_ttime_rhp;    
    int n_Bpn_s, n_Bpk_s;

    // parameters and data for amenities
    double alpha_dining ,alpha_cafes ,alpha_cafes_n ,alpha_cafes_sqm ,alpha_dining_n ,alpha_dining_sqm ,alpha_hotels ,alpha_hotels_n ,alpha_hotels_sqm ,alpha_s0 ,alpha_s0_dens ,alpha_s1 ,alpha_s1_dens ,alpha_s2 ,alpha_s2_dens; 
    double *dining, *cafes, *cafes_n, *cafes_sqm, *dining_n, *dining_sqm, *hotels, *hotels_n, *hotels_sqm, *s0, *s0_dens, *s1, *s1_dens, *s2, *s2_dens;
    double *alpha_rh;
    double *jobdensity; // Density of jobs in each region (by education) 

    //parameters, housing costs
    double *prices;
    double *rp; // per-period rate of price increases;
    double *sqm;  // size of homes in regions.   
    double psi_uc; // Translates prices into yearly usercosts of housing. 
    double psi_uc_base; // The predefined phi_h2 parameter in 1st stage demand model. Untouched during estimation of choice.

    // parameters indexing utility of residential moving cost 
    double gamma_0,gamma_a,gamma_a2,gamma_ms;
    double *gamma_c, *gamma_s; 
    int n_gamma_c, n_gamma_s;    

    // parameters indexing utility of job moving cost 
    double rho_0, rho_a, rho_unemp, *rho_s, rho_ttime_rwp, rho_ttime_rhp;

    double c_work; // disutility of working

    //parameters, housing demand (phi)
    double phi_0,phi_a,phi_a2,phi_ms,phi_y,phi_y2;
    double *phi_c, *phi_s, *phi_rh; 
    int n_phi_c, n_phi_s, n_phi_rh;

    double phi_h2   ; // parameter on quadratic term in utility of house size;  
    double phi_h2_base; // The predefined phi_h2 parameter in 1st stage demand model. Untouched during estimation of choice.

    // parameters, marginal utility of money
    double kappa_0; // intercept
    double kappa_y; // coef. on inc
    double kappa_y2; // coef. on inc squared

    double gamma_buy; // transaction cost parameter
    double gamma_sell; // monetary moving cost parameter

    double lambda_wl; // scale parameter, work location choice
    double lambda_rl; // scale parameter, residential location choice

    //  Flag indicating whether parsed phi and kappa params are scaled by phi_h2_base (and also psi_uc_base for Kappa).
    int reduced_form_sqm;

    // transition matrices
    double *trpr;
    int m_trpr,n_trpr;

    double *survival;
    int m_survival;

    // parameters tax equations
    double *thresholds, *rates;
    int ntaxbracket;

    // Precomputed income
    double *inc;

    // parameters only to be used in c-code
    int nds; // number of discrete states
    int ndc; // number of discrete choices
    int nevf; // number of elements in ev

    // blocks to summarize offset of state index (see f_IDEVF )
    int bl_inc, bl_s, bl_m, bl_couple, bl_kids, bl_rwp, bl_rhp, bl_t;

} param; //structure of model parameters (for argument passing)


// Structure containing points in the state space; e.g. for an empirical data matrix. 
typedef struct statedatastruct
{
    int it;         // Age,t.
    int irh;        // Home region.
    int irw;        // Work region.
    int ikids;      // Kids state.
    int icouple;    // Couple state.
    int im;         // Macro state.
    int is;         // School state.
    int iinc;       // Income type.
    int imove;      // Move type.
    int ok;         // Data within range
} statedata; 

// Structe containing points in choice space, e.g. for an empirical data matrix.
typedef struct choicedatastruct
{
    int irh; // Region chosen for residence.
    int irw; // Region chosen for work. 
} choicedata;


// Declaration of functions in this file

static void f_cleanparameters(param* p,const bool compute_derivatives, const int iy);
static int f_IDEVF(const param* p, declSTATES);
static int f_IDEVF_INV(const param* p, const int ievf, declpointSTATES);
static int f_IDCHOICE(const param* p, declOUTCOMES);
static void f_IDCHOICE_INV(const param* p, const int ichoice, declpointOUTCOMES);
static statedata* statedata_input(const mxArray *input);
static choicedata* choicedata_input(const mxArray *input);
static void save_statedata(const param* p, const statedata *sp, const int i, const int N, double *sim_states_out);
static void save_choicedata(const param* p, const choicedata *cp, const int i, const int N, double *sim_choice_out);
int f_get_ninc(const param* p);
int f_inc_idx(const param* p, declSTATES, declOUTCOMES);


static void f_parseparameters(param* p, const mxArray *mxparam, const bool compute_derivatives, const int iy) {

    // FIXME: NEEDS comments
    // updates parameters p and return int to reflect of this happened. 
    // if f_parseparameters are called with iy == INITIAL_PARSE it returns 0 and only initial params are parsed
    // if f_parseparameters are called with iy >=0  it returns 1 and all parameters are updated

    char *name1, *name2;
    mxArray *tmp, *param_tmp;
    double tmp_dbl, *full_Pk, *full_Pk_label, tmp_scalar;
    double *tmp_dbl_pt;
    int ntr, max_age_Pk, min_age_Pk, max_kids_Pk, i, j, idx, c, k, m, index0, ninc;

    if (iy==INITIAL_PARSE)
    {   // INITIAL PARSE FOR ALL THE YEARS
        
        // ******************   
        // model dimensions
        // ****************** 

        //  Time horizon
        tmp=mxGetField(mxparam,0,"T");
        if (!tmp) mexErrMsgTxt("T: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->T = DBL2INT(tmp_dbl);

        // initial age
        tmp=mxGetField(mxparam,0,"t0");
        if (!tmp) mexErrMsgTxt("t0: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->t0 = DBL2INT(tmp_dbl);
        
        // number of time periods
        p->lifespan=p->T-p->t0+1;

        // number of simulated years
        tmp=mxGetField(mxparam,0,"sim_T");
        if (!tmp) mexErrMsgTxt("sim_T: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->sim_T = DBL2INT(tmp_dbl);

        // initial calendar year (simulated period)
        tmp=mxGetField(mxparam,0,"y0");
        if (!tmp) mexErrMsgTxt("y0: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->y0 = DBL2INT(tmp_dbl);

        // number of individual moving cost types
        tmp=mxGetField(mxparam,0,"nmovec");
        if (!tmp) mexErrMsgTxt("nmovec: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->nmovec = DBL2INT(tmp_dbl);

        // number of income types
        tmp=mxGetField(mxparam,0,"ninc");
        if (!tmp) mexErrMsgTxt("ninc: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->ninc = DBL2INT(tmp_dbl);

        // number of schooling levels
        tmp=mxGetField(mxparam,0,"ns");
        if (!tmp) mexErrMsgTxt("ns: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->ns = DBL2INT(tmp_dbl);

        // number of macro states
        tmp=mxGetField(mxparam,0,"nm");
        if (!tmp) mexErrMsgTxt("nm: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->nm = DBL2INT(tmp_dbl);
     
        // max number of kids
        tmp=mxGetField(mxparam,0,"nkids");
        if (!tmp) mexErrMsgTxt("nkids: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->nkids = DBL2INT(tmp_dbl);

        // number of couple states
        tmp=mxGetField(mxparam,0,"ncouple");
        if (!tmp) mexErrMsgTxt("ncouple: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->ncouple = DBL2INT(tmp_dbl);

        // Number of work regions + unemployment choice
        tmp=mxGetField(mxparam,0,"nrwork");
        if (!tmp) mexErrMsgTxt("nrwork: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->nrwork = DBL2INT(tmp_dbl);
        if (p->nrwork>MAXREGIONS) mexErrMsgTxt("Update MAXREGIONS in definitions.c constant before proceeded to avoid crash!");

        // number of home regions (no unemployment choice)
        tmp=mxGetField(mxparam,0,"nrhome");
        if (!tmp) mexErrMsgTxt("nrhome: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->nrhome = DBL2INT(tmp_dbl);
        if (p->nrhome>MAXREGIONS) mexErrMsgTxt("Update MAXREGIONS in definitions.c constant before proceeded to avoid crash!");
        
        // number of state variables
        tmp=mxGetField(mxparam,0,"nstatevar");
        if (!tmp) mexErrMsgTxt("nstatevar: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->nstatevar = DBL2INT(tmp_dbl);
       
        // number of choice variables
        tmp=mxGetField(mxparam,0,"nchoicevar");
        if (!tmp) mexErrMsgTxt("nchoicevar: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->nchoicevar = DBL2INT(tmp_dbl);

        // ====================
        // State space
        // ====================
        /*dimentions for the blocks in the matrix of value function*/
        p->bl_inc     = p->nmovec;               // p->nmovec 
        p->bl_s       = p->ninc*p->bl_inc;       // p->ninc*p->nmovec
        p->bl_m       = p->ns*p->bl_s;           // p->ns*p->ninc*p->nmovec
        p->bl_couple  = p->nm*p->bl_m;           // p->nm*p->ns*p->ninc*p->nmovec
        p->bl_kids    = p->ncouple*p->bl_couple; // p->ncouple*p->nm*p->ns*p->ninc*p->nmovec
        p->bl_rwp     = p->nkids*p->bl_kids;     // p->nkids*p->ncouple*p->nm*p->ns*p->ninc*p->nmovec
        p->bl_rhp     = p->nrwork*p->bl_rwp;     // p->nrwork*p->nkids*p->ncouple*p->nm*p->ns*p->ninc*p->nmovec
        p->bl_t       = p->nrhome*p->bl_rhp;     // p->nrhome*p->nrwork*p->nkids*p->ncouple*p->nm*p->ns*p->ninc*p->nmovec
        p->nstate     = p->lifespan*p->bl_t;     // p->lifespan*p->nrhome*p->nrwork*p->nkids*p->ncouple*p->nm*p->ns*p->ninc*p->nmovec*

        // Decision space
        p->nchoice = p->nrhome*p->nrwork; // number of choice combos

        // ******************   
        // data dimensions
        // ****************** 

        // Calendar years found in data set passed through run_evaluate_fct.
        tmp=mxGetField(mxparam,0,"years_data");
        if (!tmp) mexErrMsgTxt("years_data: Input structure is incomplete!");
        p->n_years_data = (int) mxGetNumberOfElements(tmp);
        tmp_dbl_pt = (double*) mxGetPr(tmp);
        p->years_data = (int*) calloc(p->n_years_data,sizeof(int));
        for (i=0; i<p->n_years_data; i++)
            p->years_data[i] = DBL2INT(tmp_dbl_pt[i]);

        // *****************************   
        // derivatives, pnames, etc
        // ***************************** 

        //  pnames
        tmp=mxGetField(mxparam,0,"pnames");
        if (!tmp) mexErrMsgTxt("pnames: Input structure is incomplete!");
        p->n_pnames  = (int) mxGetNumberOfElements(tmp);

        p->n_param=0;
        if (compute_derivatives)
        {
            for (i = 0; i < p->n_pnames; ++i){
                p->pnames[i]   = (char*) mxArrayToString(mxGetCell(tmp, i));    
                param_tmp=mxGetField(mxparam,0,p->pnames[i]);
                if (!param_tmp) {
                    p->ni_param[i]=0;
                    p->n_param+=p->ni_param[i];
                }
                else{
                    p->ni_param[i]=mxGetNumberOfElements(param_tmp);
                    p->n_param+=p->ni_param[i];
                }
            }        
        }
        else{
            p->n_param  = 0;     // the switch for whether to compute gradients
            p->n_pnames = 0; 
        }

    }
    else {

        // ***********************************************************
        // When f_parse_parameters is called for particular year
        // ***********************************************************

        // Get current calendar year from mp given iy
        p->year = p->years_data[iy];
        p->iy = iy;

        p->iparname =(int*) calloc(p->n_param, sizeof(int));
        p->iparkey =(int*) calloc(p->n_param, sizeof(int));
        p->isubpar =(int*) calloc(p->n_param, sizeof(int));

        idx=0;
        for (i = 0; i < p->n_pnames; ++i){
            for (j = 0; j < p->ni_param[i]; ++j){
                p->iparname[idx]=i;
                p->isubpar[idx]=j; 
                p->iparkey[idx]=param_key(p->pnames[i]); 
                idx+=1;
            }    
        }

        // Flag indicating whether parsed phi and kappa params are scaled by phi_h2_base and psi_uc_base
        tmp=mxGetField(mxparam,0,"reduced_form_sqm");
        if (!tmp) mexErrMsgTxt("reduced_form_sqm: Input structure is incomplete!");
        tmp_dbl = (double) mxGetScalar(tmp);
        p->reduced_form_sqm = DBL2INT(tmp_dbl);

        // scale of monetary units. 100,000s etc.
        tmp=mxGetField(mxparam,0,"dkk_scale");
        if (!tmp) mexErrMsgTxt("dkk_scale: Input structure is incomplete!");
        p->dkk_scale = (double) mxGetScalar(tmp); 

        // scaling param phi_h2 in 1st stage estimation. Serves to rescale sqm demand during estimation. 
        tmp=mxGetField(mxparam,0,"phi_h2_base");
        if (!tmp) mexErrMsgTxt("phi_h2_base: Input structure is incomplete!");
        p->phi_h2_base = (double) mxGetScalar(tmp); 

         // ******************   
         // model parameters
         // ******************   

        // parameter with travel time 
        tmp=mxGetField(mxparam,0,"eta_ttime");
        if (!tmp) mexErrMsgTxt("eta_ttime: Input structure is incomplete!");
        p->eta_ttime = (double) mxGetScalar(tmp);

        // parameter with time^2 in minutes
        tmp=mxGetField(mxparam,0,"eta_ttime2");
        if (!tmp) mexErrMsgTxt("eta_ttime2: Input structure is incomplete!");
        p->eta_ttime2 = (double) mxGetScalar(tmp);

        // discount factor
        tmp=mxGetField(mxparam,0,"delta");
        if (!tmp) mexErrMsgTxt("delta: Input structure is incomplete!");
        p->delta = (double) mxGetScalar(tmp);


        // **** Parameters of income equation **********
        // nested structure beta_empl
        
        tmp=mxGetField(mxparam,0,"beta_empl");
        if (!tmp) mexErrMsgTxt("beta_cons: Input structure is incomplete!"); //once!
        tmp=mxGetField(tmp,0,"beta_cons");
        if (!tmp) mexErrMsgTxt("beta_cons: Input structure is incomplete!");
        p->beta_cons = (double*) mxGetPr(tmp);
        if ((int) mxGetNumberOfElements(tmp) != p->ns) //by schooling levels
             mexErrMsgTxt("beta_cons: Input structure is corrupt, wrong input size!");

        tmp=mxGetField(mxparam,0,"beta_empl");
        tmp=mxGetField(tmp,0,"beta_rw");
        if (!tmp) mexErrMsgTxt("beta_rw: Input structure is incomplete!");
        p->beta_rw = (double*) mxGetPr(tmp);
        if ((int) mxGetNumberOfElements(tmp) != p->nrhome*p->ns) //by school-locations
             mexErrMsgTxt("beta_rw: Input structure is corrupt, wrong input size!");

        tmp=mxGetField(mxparam,0,"beta_empl");
        tmp=mxGetField(tmp,0,"beta_l_unemp");
        if (!tmp) mexErrMsgTxt("beta_l_unemp: Input structure is incomplete!");
        p->beta_l_unemp = (double*) mxGetPr(tmp);
        if ((int) mxGetNumberOfElements(tmp) != p->ns) //by schooling levels
             mexErrMsgTxt("beta_l_unemp: Input structure is corrupt, wrong input size!");

        tmp=mxGetField(mxparam,0,"beta_empl");
        tmp=mxGetField(tmp,0,"beta_l_unemp_rw");
        if (!tmp) mexErrMsgTxt("beta_l_unemp_rw: Input structure is incomplete!");
        p->beta_l_unemp_rw = (double*) mxGetPr(tmp);
        if ((int) mxGetNumberOfElements(tmp) != p->nrhome*p->ns) //by school-locations
             mexErrMsgTxt("beta_l_unemp_rw: Input structure is corrupt, wrong input size!");

        tmp=mxGetField(mxparam,0,"beta_empl");
        tmp=mxGetField(tmp,0,"n_beta_age_pol");
        if (!tmp) mexErrMsgTxt("n_beta_age_pol: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->n_beta_age_pol = DBL2INT(tmp_dbl);

        if (p->n_beta_age_pol > 0)
        {
            tmp=mxGetField(mxparam,0,"beta_empl");
            tmp=mxGetField(tmp,0,"beta_age_pol");
            if (!tmp) mexErrMsgTxt("beta_age_pol: Input structure is incomplete!");
            p->beta_age_pol = (double*) mxGetPr(tmp);
            if ((int) mxGetNumberOfElements(tmp) != p->n_beta_age_pol)
                mexErrMsgTxt("beta_age_pol: Input structure is corrupt, wrong input size!");
        }

        tmp=mxGetField(mxparam,0,"beta_empl");
        tmp=mxGetField(tmp,0,"n_beta_age_fe");
        if (!tmp) mexErrMsgTxt("n_beta_age_fe: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->n_beta_age_fe = DBL2INT(tmp_dbl);
        // check dimension?

        if (p->n_beta_age_fe > 0)
        {
            tmp=mxGetField(mxparam,0,"beta_empl");
            tmp=mxGetField(tmp,0,"beta_age_fe");
            if (!tmp) mexErrMsgTxt("beta_age_fe: Input structure is incomplete!");
            p->beta_age_fe = (double*) mxGetPr(tmp);     
            if ((int) mxGetNumberOfElements(tmp) != p->n_beta_age_fe)
                mexErrMsgTxt("beta_age_fe: Input structure is corrupt, wrong input size!");
        }

        tmp=mxGetField(mxparam,0,"beta_empl");
        tmp=mxGetField(tmp,0,"beta_min_age_fe");
        if (!tmp) mexErrMsgTxt("beta_min_age_fe: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->beta_min_age_fe = DBL2INT(tmp_dbl);


        // **** Parameters of non-employment income **********
        // nested structure beta_unempl

        tmp=mxGetField(mxparam,0,"beta_unempl");
        if (!tmp) mexErrMsgTxt("beta_unempl: Input structure is incomplete!"); //once
        tmp=mxGetField(tmp,0,"b_cons");
        if (!tmp) mexErrMsgTxt("b_cons: Input structure is incomplete!");
        p->b_0 = (double*) mxGetPr(tmp);
        if ((int) mxGetNumberOfElements(tmp) != p->ns) //by schooling levels
             mexErrMsgTxt("beta_cons: Input structure is corrupt, wrong input size!");

        tmp=mxGetField(mxparam,0,"beta_unempl");
        tmp=mxGetField(tmp,0,"n_b_age_pol");
        if (!tmp) mexErrMsgTxt("n_b_age_pol: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->n_b_age_pol = DBL2INT(tmp_dbl);

        if (p->n_b_age_pol > 0)
        {
            tmp=mxGetField(mxparam,0,"beta_unempl");
            tmp=mxGetField(tmp,0,"b_age_pol");
            if (!tmp) mexErrMsgTxt("b_age_pol: Input structure is incomplete!");
            p->b_age_pol = (double*) mxGetPr(tmp);
            if ((int) mxGetNumberOfElements(tmp) != p->n_b_age_pol)
                mexErrMsgTxt("b_age_pol: Input structure is corrupt, wrong input size!");
        }

        tmp=mxGetField(mxparam,0,"beta_unempl");
        tmp=mxGetField(tmp,0,"n_b_age_fe");
        if (!tmp) mexErrMsgTxt("n_b_age_fe: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->n_b_age_fe = DBL2INT(tmp_dbl);

        if (p->n_b_age_fe > 0)
        {
            tmp=mxGetField(mxparam,0,"beta_unempl");
            tmp=mxGetField(tmp,0,"b_age_fe");
            if (!tmp) mexErrMsgTxt("b_age_fe: Input structure is incomplete!");
            p->b_age_fe = (double*) mxGetPr(tmp);
            if ((int) mxGetNumberOfElements(tmp) != p->n_b_age_fe)
                mexErrMsgTxt("b_age_fe: Input structure is corrupt, wrong input size!");
        }

        tmp=mxGetField(mxparam,0,"beta_unempl");
        tmp=mxGetField(tmp,0,"b_min_age_fe");
        if (!tmp) mexErrMsgTxt("b_min_age_fe: Input structure is incomplete!");
        tmp_dbl = mxGetScalar(tmp);
        p->b_min_age_fe = DBL2INT(tmp_dbl);


        // **** Parameters for job transitions **********

        // coefficients indexing probability of succeeding in getting a new job.
        // coefficients prob. of new job, constant
        tmp=mxGetField(mxparam,0,"Bpn_0");
        if (!tmp) mexErrMsgTxt("Bpn_0: Input structure is incomplete!");
        p->Bpn_0 = (double) mxGetScalar(tmp);

         // coefficients prob. of new job, age
        tmp=mxGetField(mxparam,0,"Bpn_a");
        if (!tmp) mexErrMsgTxt("Bpn_a: Input structure is incomplete!");
        p->Bpn_a = (double) mxGetScalar(tmp);

        // coefficients prob. of new job, age
        tmp=mxGetField(mxparam,0,"Bpn_unemp");
        if (!tmp) mexErrMsgTxt("Bpn_a: Input structure is incomplete!");
        p->Bpn_unemp = (double) mxGetScalar(tmp);

        // coefficients prob. of new job, education dummies
        tmp=mxGetField(mxparam,0,"Bpn_s");
        if (!tmp) mexErrMsgTxt("Bpn_s: Input structure is incomplete!");
        p->Bpn_s = (double*)  mxGetPr(tmp);
        p->n_Bpn_s = (int) mxGetNumberOfElements(tmp); 

         // coefficients prob. of new job, travel time between irwp and irw
        tmp=mxGetField(mxparam,0,"Bpn_ttime_rwp");
        if (!tmp) mexErrMsgTxt("Bpn_ttime_rwp: Input structure is incomplete!");
        p->Bpn_ttime_rwp = (double) mxGetScalar(tmp);

        tmp=mxGetField(mxparam,0,"Bpn_ttime_rhp");
        if (!tmp) mexErrMsgTxt("Bpn_ttime_rhp: Input structure is incomplete!");
        p->Bpn_ttime_rhp = (double) mxGetScalar(tmp);

        // coefficients indexing probability of keeping current job   
        // coefficients prob. of keeping job, constant
        tmp=mxGetField(mxparam,0,"Bpk_0");
        if (!tmp) mexErrMsgTxt("Bpk_0: Input structure is incomplete!");
        p->Bpk_0 = (double) mxGetScalar(tmp);

         // coefficients prob. of keeping job, age
        tmp=mxGetField(mxparam,0,"Bpk_a");
        if (!tmp) mexErrMsgTxt("Bpk_a: Input structure is incomplete!");
        p->Bpk_a = (double) mxGetScalar(tmp);

        // coefficients prob. of keeping job, age
        tmp=mxGetField(mxparam,0,"Bpk_unemp");
        if (!tmp) mexErrMsgTxt("Bpk_a: Input structure is incomplete!");
        p->Bpk_unemp = (double) mxGetScalar(tmp);

        // coefficient on job density in prb of getting a job
        tmp=mxGetField(mxparam,0,"Bpn_jobdensity");
        if (!tmp) mexErrMsgTxt("Bpn_jobdensity: Input structure is incomplete!");
        p->Bpn_jobdensity = (double)  mxGetScalar(tmp);    

        // coefficients prob. of keeping job, education dummies
        tmp=mxGetField(mxparam,0,"Bpk_s");
        if (!tmp) mexErrMsgTxt("Bpk_s: Input structure is incomplete!");
        p->Bpk_s = (double*)  mxGetPr(tmp);
        p->n_Bpk_s = (int) mxGetNumberOfElements(tmp); 

        // coefficient on job density in prb of keeping a job
        tmp=mxGetField(mxparam,0,"Bpk_jobdensity");
        if (!tmp) mexErrMsgTxt("Bpk_jobdensity: Input structure is incomplete!");
        p->Bpk_jobdensity = (double)  mxGetScalar(tmp);

        tmp=mxGetField(mxparam,0,"Bpk_ttime_rhp");
        if (!tmp) mexErrMsgTxt("Bpk_ttime_rhp: Input structure is incomplete!");
        p->Bpk_ttime_rhp = (double) mxGetScalar(tmp);

        // coefficients for housing demand, residential region dummies
        tmp=mxGetField(mxparam,0,"alpha_rh");
        if (!tmp) mexErrMsgTxt("alpha_rh: Input structure is incomplete!");
        p->alpha_rh = (double*)  mxGetPr(tmp);

        // -------------------------------------------
        //  Amenities
        //  nested structure "amenities" -> cell array by years
        // -------------------------------------------

        tmp=mxGetField(mxparam,0,"alpha_dining");
        if (!tmp) mexErrMsgTxt("alpha_dining: Input structure is incomplete!");
        p->alpha_dining = (double) mxGetScalar(tmp);

        // amenitity: dining
        tmp=mxGetField(mxparam,0,"amenities");
        if (!tmp) mexErrMsgTxt("amenities: Input structure is incomplete!"); //once
        tmp=mxGetField(tmp,0,"dining");
        if (!tmp) mexErrMsgTxt("dining: Input structure is incomplete!");
        tmp = mxGetCell(tmp, iy); 
        p->dining = (double*) mxGetPr(tmp); 
        if ((int) mxGetNumberOfElements(tmp) != p->nrhome) //by regions
             mexErrMsgTxt("amenities.dining: Input structure is corrupt, wrong input size!");

        // coefficients to amenities: cafes
        tmp=mxGetField(mxparam,0,"alpha_cafes");
        if (!tmp) mexErrMsgTxt("alpha_cafes: Input structure is incomplete!");
        p->alpha_cafes = (double) mxGetScalar(tmp);

        // amenitity: cafes
        tmp=mxGetField(mxparam,0,"amenities");
        tmp=mxGetField(tmp,0,"cafes");
        if (!tmp) mexErrMsgTxt("cafes: Input structure is incomplete!");
        tmp = mxGetCell(tmp, iy); 
        p->cafes = (double*) mxGetPr(tmp); 
        if ((int) mxGetNumberOfElements(tmp) != p->nrhome) //by regions
             mexErrMsgTxt("amenities.cafes: Input structure is corrupt, wrong input size!");


        // coefficients to amenities: cafes_n
        tmp=mxGetField(mxparam,0,"alpha_cafes_n");
        if (!tmp) mexErrMsgTxt("alpha_cafes_n: Input structure is incomplete!");
        p->alpha_cafes_n = (double) mxGetScalar(tmp);

        // amenitity: cafes_n
        tmp=mxGetField(mxparam,0,"amenities");
        tmp=mxGetField(tmp,0,"cafes_n");
        if (!tmp) mexErrMsgTxt("cafes_n: Input structure is incomplete!");
        tmp = mxGetCell(tmp, iy); 
        p->cafes_n = (double*) mxGetPr(tmp); 
        if ((int) mxGetNumberOfElements(tmp) != p->nrhome) //by regions
             mexErrMsgTxt("amenities.cafes_n: Input structure is corrupt, wrong input size!");

        // coefficients to amenities: cafes_sqm
        tmp=mxGetField(mxparam,0,"alpha_cafes_sqm");
        if (!tmp) mexErrMsgTxt("alpha_cafes_sqm: Input structure is incomplete!");
        p->alpha_cafes_sqm = (double) mxGetScalar(tmp);

        // amenitity: cafes_sqm
        tmp=mxGetField(mxparam,0,"amenities");
        tmp=mxGetField(tmp,0,"cafes_sqm");
        if (!tmp) mexErrMsgTxt("cafes_sqm: Input structure is incomplete!");
        tmp = mxGetCell(tmp, iy); 
        p->cafes_sqm = (double*) mxGetPr(tmp); 
        if ((int) mxGetNumberOfElements(tmp) != p->nrhome) //by regions
             mexErrMsgTxt("amenities.cafes_sqm: Input structure is corrupt, wrong input size!");

        // coefficients to amenities: dining_n
        tmp=mxGetField(mxparam,0,"alpha_dining_n");
        if (!tmp) mexErrMsgTxt("alpha_dining_n: Input structure is incomplete!");
        p->alpha_dining_n = (double) mxGetScalar(tmp);

        // amenitity: dining_n
        tmp=mxGetField(mxparam,0,"amenities");
        tmp=mxGetField(tmp,0,"dining_n");
        if (!tmp) mexErrMsgTxt("dining_n: Input structure is incomplete!");
        tmp = mxGetCell(tmp, iy); 
        p->dining_n = (double*) mxGetPr(tmp); 
        if ((int) mxGetNumberOfElements(tmp) != p->nrhome) //by regions
             mexErrMsgTxt("amenities.dining_n: Input structure is corrupt, wrong input size!");

        // coefficients to amenities: dining_sqm
        tmp=mxGetField(mxparam,0,"alpha_dining_sqm");
        if (!tmp) mexErrMsgTxt("alpha_dining_sqm: Input structure is incomplete!");
        p->alpha_dining_sqm = (double) mxGetScalar(tmp);

        // amenitity: dining_sqm
        tmp=mxGetField(mxparam,0,"amenities");
        tmp=mxGetField(tmp,0,"dining_sqm");
        if (!tmp) mexErrMsgTxt("dining_sqm: Input structure is incomplete!");
        tmp = mxGetCell(tmp, iy); 
        p->dining_sqm = (double*) mxGetPr(tmp); 
        if ((int) mxGetNumberOfElements(tmp) != p->nrhome) //by regions
             mexErrMsgTxt("amenities.dining_sqm: Input structure is corrupt, wrong input size!");

        // coefficients to amenities: hotels
        tmp=mxGetField(mxparam,0,"alpha_hotels");
        if (!tmp) mexErrMsgTxt("alpha_hotels: Input structure is incomplete!");
        p->alpha_hotels = (double) mxGetScalar(tmp);

        // amenitity: hotels
        tmp=mxGetField(mxparam,0,"amenities");
        tmp=mxGetField(tmp,0,"hotels");
        if (!tmp) mexErrMsgTxt("hotels: Input structure is incomplete!");
        tmp = mxGetCell(tmp, iy); 
        p->hotels = (double*) mxGetPr(tmp); 
        if ((int) mxGetNumberOfElements(tmp) != p->nrhome) //by regions
             mexErrMsgTxt("amenities.hotels: Input structure is corrupt, wrong input size!");

        // coefficients to amenities: hotels_n
        tmp=mxGetField(mxparam,0,"alpha_hotels_n");
        if (!tmp) mexErrMsgTxt("alpha_hotels_n: Input structure is incomplete!");
        p->alpha_hotels_n = (double) mxGetScalar(tmp);

        // amenitity: hotels_n
        tmp=mxGetField(mxparam,0,"amenities");
        tmp=mxGetField(tmp,0,"hotels_n");
        if (!tmp) mexErrMsgTxt("hotels_n: Input structure is incomplete!");
        tmp = mxGetCell(tmp, iy); 
        p->hotels_n = (double*) mxGetPr(tmp); 
        if ((int) mxGetNumberOfElements(tmp) != p->nrhome) //by regions
             mexErrMsgTxt("amenities.hotels_n: Input structure is corrupt, wrong input size!");

        // coefficients to amenities: hotels_sqm
        tmp=mxGetField(mxparam,0,"alpha_hotels_sqm");
        if (!tmp) mexErrMsgTxt("alpha_hotels_sqm: Input structure is incomplete!");
        p->alpha_hotels_sqm = (double) mxGetScalar(tmp);

        // amenitity: hotels_sqm
        tmp=mxGetField(mxparam,0,"amenities");
        tmp=mxGetField(tmp,0,"hotels_sqm");
        if (!tmp) mexErrMsgTxt("hotels_sqm: Input structure is incomplete!");
        tmp = mxGetCell(tmp, iy); 
        p->hotels_sqm = (double*) mxGetPr(tmp); 
        if ((int) mxGetNumberOfElements(tmp) != p->nrhome) //by regions
             mexErrMsgTxt("amenities.hotels_sqm: Input structure is corrupt, wrong input size!");

        // coefficients to amenities: s0
        tmp=mxGetField(mxparam,0,"alpha_s0");
        if (!tmp) mexErrMsgTxt("alpha_s0: Input structure is incomplete!");
        p->alpha_s0 = (double) mxGetScalar(tmp);

        // amenitity: s0
        tmp=mxGetField(mxparam,0,"amenities");
        tmp=mxGetField(tmp,0,"s0");
        if (!tmp) mexErrMsgTxt("s0: Input structure is incomplete!");
        tmp = mxGetCell(tmp, iy); 
        p->s0 = (double*) mxGetPr(tmp); 
        if ((int) mxGetNumberOfElements(tmp) != p->nrhome) //by regions
             mexErrMsgTxt("amenities.s0: Input structure is corrupt, wrong input size!");

        // coefficients to amenities: s0_dens
        tmp=mxGetField(mxparam,0,"alpha_s0_dens");
        if (!tmp) mexErrMsgTxt("alpha_s0_dens: Input structure is incomplete!");
        p->alpha_s0_dens = (double) mxGetScalar(tmp);

        // amenitity: s0_dens
        tmp=mxGetField(mxparam,0,"amenities");
        tmp=mxGetField(tmp,0,"s0_dens");
        if (!tmp) mexErrMsgTxt("s0_dens: Input structure is incomplete!");
        tmp = mxGetCell(tmp, iy); 
        p->s0_dens = (double*) mxGetPr(tmp); 
        if ((int) mxGetNumberOfElements(tmp) != p->nrhome) //by regions
             mexErrMsgTxt("amenities.s0_dens: Input structure is corrupt, wrong input size!");

        // coefficients to amenities: s1
        tmp=mxGetField(mxparam,0,"alpha_s1");
        if (!tmp) mexErrMsgTxt("alpha_s1: Input structure is incomplete!");
        p->alpha_s1 = (double) mxGetScalar(tmp);

        // amenitity: s1
        tmp=mxGetField(mxparam,0,"amenities");
        tmp=mxGetField(tmp,0,"s1");
        if (!tmp) mexErrMsgTxt("s1: Input structure is incomplete!");
        tmp = mxGetCell(tmp, iy); 
        p->s1 = (double*) mxGetPr(tmp); 
        if ((int) mxGetNumberOfElements(tmp) != p->nrhome) //by regions
             mexErrMsgTxt("amenities.s1: Input structure is corrupt, wrong input size!");

        // coefficients to amenities: s1_dens
        tmp=mxGetField(mxparam,0,"alpha_s1_dens");
        if (!tmp) mexErrMsgTxt("alpha_s1_dens: Input structure is incomplete!");
        p->alpha_s1_dens = (double) mxGetScalar(tmp);

        // amenitity: s1_dens
        tmp=mxGetField(mxparam,0,"amenities");
        tmp=mxGetField(tmp,0,"s1_dens");
        if (!tmp) mexErrMsgTxt("s1_dens: Input structure is incomplete!");
        tmp = mxGetCell(tmp, iy); 
        p->s1_dens = (double*) mxGetPr(tmp); 
        if ((int) mxGetNumberOfElements(tmp) != p->nrhome) //by regions
             mexErrMsgTxt("amenities.s1_dens: Input structure is corrupt, wrong input size!");

        // coefficients to amenities: s2
        tmp=mxGetField(mxparam,0,"alpha_s2");
        if (!tmp) mexErrMsgTxt("alpha_s2: Input structure is incomplete!");
        p->alpha_s2 = (double) mxGetScalar(tmp);

        // amenitity: s2
        tmp=mxGetField(mxparam,0,"amenities");
        tmp=mxGetField(tmp,0,"s2");
        if (!tmp) mexErrMsgTxt("s2: Input structure is incomplete!");
        tmp = mxGetCell(tmp, iy); 
        p->s2 = (double*) mxGetPr(tmp); 
        if ((int) mxGetNumberOfElements(tmp) != p->nrhome) //by regions
             mexErrMsgTxt("amenities.s2: Input structure is corrupt, wrong input size!");

        // coefficients to amenities: s2_dens
        tmp=mxGetField(mxparam,0,"alpha_s2_dens");
        if (!tmp) mexErrMsgTxt("alpha_s2_dens: Input structure is incomplete!");
        p->alpha_s2_dens = (double) mxGetScalar(tmp);

        // amenitity: s2_dens
        tmp=mxGetField(mxparam,0,"amenities");
        tmp=mxGetField(tmp,0,"s2_dens");
        if (!tmp) mexErrMsgTxt("s2_dens: Input structure is incomplete!");
        tmp = mxGetCell(tmp, iy); 
        p->s2_dens = (double*) mxGetPr(tmp); 
        if ((int) mxGetNumberOfElements(tmp) != p->nrhome) //by regions
             mexErrMsgTxt("amenities.s2_dens: Input structure is corrupt, wrong input size!");

        // job-density
        tmp=mxGetField(mxparam,0,"amenities");
        tmp=mxGetField(tmp,0,"jobdensity");
        if (!tmp) mexErrMsgTxt("jobdensity: Input structure is incomplete!");
        tmp = mxGetCell(tmp, iy); 
        p->jobdensity = (double*)  mxGetPr(tmp);

        // house prices per sqm
        tmp=mxGetField(mxparam,0,"amenities");
        tmp=mxGetField(tmp,0,"prices");
        if (!tmp) mexErrMsgTxt("prices: Input structure is incomplete!");
        tmp = mxGetCell(tmp, iy); 
        p->prices = (double*) mxGetPr(tmp); 

        // -------------------------------------------
        //      Housing costs
        // -------------------------------------------
           
        // parameter that translates housing price into yearly user costs. 
        tmp=mxGetField(mxparam,0,"psi_uc");
        if (!tmp) mexErrMsgTxt("psi_uc: Input structure is incomplete!");
        p->psi_uc = (double) mxGetScalar(tmp);

        // base value of psi_uc in 1st stage estimation. Serves to rescale sqm demand during estimation. 
        tmp=mxGetField(mxparam,0,"psi_uc_base");
        if (!tmp) mexErrMsgTxt("psi_uc_base: Input structure is incomplete!");
        p->psi_uc_base = (double) mxGetScalar(tmp); 

        // coefficients switching costs, moving costs
        tmp=mxGetField(mxparam,0,"gamma_sell");
        if (!tmp) mexErrMsgTxt("gamma_sell: Input structure is incomplete!");
        p->gamma_sell = (double)  mxGetScalar(tmp);

        // coefficients switching costs, transaction costs
        tmp=mxGetField(mxparam,0,"gamma_buy");
        if (!tmp) mexErrMsgTxt("gamma_buy: Input structure is incomplete!");
        p->gamma_buy = (double) mxGetScalar(tmp);


        //--- utility residential moving cost  ---------------------------------------------------------------- //

        // coefficients switching costs, constant
        tmp=mxGetField(mxparam,0,"gamma_0");
        if (!tmp) mexErrMsgTxt("gamma_0: Input structure is incomplete!");
        p->gamma_0 = (double) mxGetScalar(tmp);

        // coefficients switching costs, age
        tmp=mxGetField(mxparam,0,"gamma_a");
        if (!tmp) mexErrMsgTxt("gamma_a: Input structure is incomplete!");
        p->gamma_a = (double) mxGetScalar(tmp);

        // coefficients switching costs, age2
        tmp=mxGetField(mxparam,0,"gamma_a2");
        if (!tmp) mexErrMsgTxt("gamma_a: Input structure is incomplete!");
        p->gamma_a2 = (double) mxGetScalar(tmp);

        // coefficients switching costs, couple
        tmp=mxGetField(mxparam,0,"gamma_ms");
        if (!tmp) mexErrMsgTxt("gamma_ms: Input structure is incomplete!");
        p->gamma_ms = (double) mxGetScalar(tmp);

        // coefficients switching costs, education dummies
        tmp=mxGetField(mxparam,0,"gamma_s");
        if (!tmp) mexErrMsgTxt("gamma_s: Input structure is incomplete!");
        p->gamma_s = (double*)  mxGetPr(tmp);
        p->n_gamma_s = (int) mxGetNumberOfElements(tmp); 

        // coefficients switching costs, child dummies
        tmp=mxGetField(mxparam,0,"gamma_c");
        if (!tmp) mexErrMsgTxt("gamma_c: Input structure is incomplete!");
        p->gamma_c = (double*)  mxGetPr(tmp);
        p->n_gamma_c = (int) mxGetNumberOfElements(tmp); 


        //--- job moving cost (measured in utility)  ---------------------------------------------------------------- //

        // coefficients of new job, constant
        tmp=mxGetField(mxparam,0,"rho_0");
        if (!tmp) mexErrMsgTxt("rho_0: Input structure is incomplete!");
        p->rho_0 = (double) mxGetScalar(tmp);

         // coefficients prob. of new job, age
        tmp=mxGetField(mxparam,0,"rho_a");
        if (!tmp) mexErrMsgTxt("rho_a: Input structure is incomplete!");
        p->rho_a = (double) mxGetScalar(tmp);

        // coefficients prob. of new job, age
        tmp=mxGetField(mxparam,0,"rho_unemp");
        if (!tmp) mexErrMsgTxt("rho_unemp: Input structure is incomplete!");
        p->rho_unemp = (double) mxGetScalar(tmp);

        // coefficients prob. of new job, education dummies
        tmp=mxGetField(mxparam,0,"rho_s");
        if (!tmp) mexErrMsgTxt("rho_s: Input structure is incomplete!");
        p->rho_s = (double*)  mxGetPr(tmp);

         // coefficients prob. of new job, travel time between irwp and irw
        tmp=mxGetField(mxparam,0,"rho_ttime_rwp");
        if (!tmp) mexErrMsgTxt("rho_ttime_rwp: Input structure is incomplete!");
        p->rho_ttime_rwp = (double) mxGetScalar(tmp);

        tmp=mxGetField(mxparam,0,"rho_ttime_rhp");
        if (!tmp) mexErrMsgTxt("rho_ttime_rhp: Input structure is incomplete!");
        p->rho_ttime_rhp = (double) mxGetScalar(tmp);

        //--- disutility of working , cw ---------------------------------------------------------------- //

        // coefficients of new job, constant 
        tmp=mxGetField(mxparam,0,"c_work");
        if (!tmp) mexErrMsgTxt("c_work: Input structure is incomplete!");
        p->c_work = (double) mxGetScalar(tmp);

        //--- utility of house size, h ---------------------------------------------------------------- //

        // coefficients for housing demand, intercept
        tmp=mxGetField(mxparam,0,"phi_0");
        if (!tmp) mexErrMsgTxt("phi_0: Input structure is incomplete!");
        p->phi_0 = (double) mxGetScalar(tmp);

        // coefficients for housing demand, residential region dummies
        tmp=mxGetField(mxparam,0,"phi_rh");
        if (!tmp) mexErrMsgTxt("phi_rh: Input structure is incomplete!");
        p->phi_rh = (double*)  mxGetPr(tmp);
        p->n_phi_rh = (int) mxGetNumberOfElements(tmp); 

        // coefficients for housing demand, age
        tmp=mxGetField(mxparam,0,"phi_a");
        if (!tmp) mexErrMsgTxt("phi_a: Input structure is incomplete!");
        p->phi_a = (double) mxGetScalar(tmp);

        // coefficients for housing demand, age^2
        tmp=mxGetField(mxparam,0,"phi_a2");
        if (!tmp) mexErrMsgTxt("phi_a2: Input structure is incomplete!");
        p->phi_a2 = (double) mxGetScalar(tmp);

        // coefficients for housing demand, marital status
        tmp=mxGetField(mxparam,0,"phi_ms");
        if (!tmp) mexErrMsgTxt("phi_ms: Input structure is incomplete!");
        p->phi_ms = (double) mxGetScalar(tmp);

        // coefficients for housing demand, income
        tmp=mxGetField(mxparam,0,"phi_y");
        if (!tmp) mexErrMsgTxt("phi_y: Input structure is incomplete!");
        p->phi_y = (double) mxGetScalar(tmp);

        // coefficients for housing demand, income^2
        tmp=mxGetField(mxparam,0,"phi_y2");
        if (!tmp) mexErrMsgTxt("phi_y2: Input structure is incomplete!");
        p->phi_y2 = (double) mxGetScalar(tmp);

        // coefficients for housing demand, education dummies
        tmp=mxGetField(mxparam,0,"phi_s");
        if (!tmp) mexErrMsgTxt("phi_s: Input structure is incomplete!");
        p->phi_s = (double*)  mxGetPr(tmp);
        p->n_phi_s = (int) mxGetNumberOfElements(tmp); 

        // coefficients for housing demand, child dummies
        tmp=mxGetField(mxparam,0,"phi_c");
        if (!tmp) mexErrMsgTxt("phi_c: Input structure is incomplete!");
        p->phi_c = (double*)  mxGetPr(tmp);
        p->n_phi_c = (int) mxGetNumberOfElements(tmp); 

        // coefficients for housing demand,  parameter on qudratic term in utiliy of house size  
        tmp=mxGetField(mxparam,0,"phi_h2");
        if (!tmp) mexErrMsgTxt("phi_h2: Input structure is incomplete!");
        p->phi_h2 = (double) mxGetScalar(tmp); 

        // If optimizer tries phi_h2 in positive domain, then transform back to negative.
        // if (p->is_estimating & p->phi_h2 > 0) p->phi_h2 = (-1)*abs((-1)*p->phi_h2); 
         // if (p->phi_h2 > 0) p->phi_h2 = (-1)*abs((-1)*p->phi_h2); 
        p->phi_h2 = (-1)*fabs((-1)*p->phi_h2); 

        //--- marginal utility of money  ---------------------------------------------------------------- //

         // coefficients marginal utility of money, intercept
        tmp=mxGetField(mxparam,0,"kappa_0");
        if (!tmp) mexErrMsgTxt("kappa_0: Input structure is incomplete!");
        p->kappa_0 = (double) mxGetScalar(tmp);

         // coefficients marginal utility of money, income
        tmp=mxGetField(mxparam,0,"kappa_y");
        if (!tmp) mexErrMsgTxt("kappa_y:  structure is incomplete!");
        p->kappa_y = (double) mxGetScalar(tmp);

         // coefficients marginal utility of money, income^2
        tmp=mxGetField(mxparam,0,"kappa_y2");
        if (!tmp) mexErrMsgTxt("kappa_y2: Input structure is incomplete!");
        p->kappa_y2 = (double) mxGetScalar(tmp);

        //--- taste shocks  ---------------------------------------------------------------- //
        
        // scale parameter, work location choice. 
        tmp=mxGetField(mxparam,0,"lambda_wl");
        if (!tmp) mexErrMsgTxt("lambda_wl: Input structure is incomplete!");
        p->lambda_wl = (double) mxGetScalar(tmp);

        // scale parameter for iid extreme value shocks, residental location hoice. 
        tmp=mxGetField(mxparam,0,"lambda_rl");
        if (!tmp) mexErrMsgTxt("lambda_rl: Input structure is incomplete!");
        p->lambda_rl = (double) mxGetScalar(tmp);

        //--- Transition transportation cost  ---------------------------------------------------------------- //

        tmp=mxGetField(mxparam,0,"traveltime");
        if (!tmp) mexErrMsgTxt("traveltime: Input structure is incomplete!");
        tmp = mxGetCell(tmp, iy); 
        p->ttime = (double*) mxGetPr(tmp); 
        if ((int) mxGetNumberOfElements(tmp)!=p->nrwork*p->nrhome)
            mexErrMsgTxt("traveltime: Input structure is corrupt, wrong dimension of travel time matrix!");
        // The way to index ttime data is:
        // IDCHOICE or f_IDCHOICE(const param* p, declOUTCOMES)

        //--- Transition matrices ----------------------------------------------------------------------------- //

        tmp=mxGetField(mxparam,0,"trprobs");
        if (!tmp) mexErrMsgTxt("trprobs: Input structure is incomplete!");
        p->trpr = (double*) mxGetPr(tmp); 
        p->m_trpr = (int) mxGetM(tmp); 
        p->n_trpr = (int) mxGetN(tmp); 

        if (p->m_trpr!=p->lifespan*p->ns*p->nkids*p->ncouple)
            mexErrMsgTxt("trprob: Input structure is corrupt, wrong number of rows of transition probilities!");
        if (p->n_trpr!=p->nkids*p->ncouple)
            mexErrMsgTxt("trprob: Input structure is corrupt, wrong number of columns of transition probilities!");

        //Indexing of transition probabilities:
#define TRPRfromIDX  (p->ncouple*p->nkids*p->ns*it+ p->ncouple*p->nkids*is + p->ncouple*ikids+ icouple )
#define TRPRtoIDX  (p->ncouple*ikids1+ icouple1 )

#ifdef FOOLPROOF
        for (i=0;i<p->m_trpr;i++){
            tmp_dbl = 0.0;
            for (j=0;j<p->n_trpr;j++){
                tmp_dbl += p->trpr[i + p->m_trpr * j];
            }
            if (fabs(tmp_dbl-1.0)>TOLERANCE) {
                printf("trpr does not sum up to one in row %d of transitions table!\n",i);
                mexErrMsgTxt("Error: transition probabilities matrix corrupted");
            }
        }
#endif   

        //--- Survival probabilities ----------------------------------------------------------------------------- //

        tmp=mxGetField(mxparam,0,"survival");
        if (!tmp) mexErrMsgTxt("survival: Input structure is incomplete!");
        p->survival = (double*) mxGetPr(tmp); 
        p->m_survival = (int) mxGetM(tmp); 

        if (p->m_survival!=p->lifespan*p->ns*p->ncouple)
            mexErrMsgTxt("survival: Input structure is corrupt, wrong number of survival probabilities!");

        //Indexing of survival probabilities:
#define SURVIDX  (p->ncouple*p->ns*it+ p->ncouple*is+ icouple )

        //--- Tax functions data ----------------------------------------------------------------------------- //

        tmp=mxGetField(mxparam,0,"tax_data");
        if (!tmp) mexErrMsgTxt("tax_data: Input structure is incomplete!");
        if ((int) mxGetM(tmp) != p->nrhome)
             mexErrMsgTxt("tax_data: Input structure is corrupt, wrong number of regions");
        p->ntaxbracket = (int) mxGetN(tmp)/2;
        p->thresholds = (double*) mxGetPr(tmp); 
        p->rates = p->thresholds + p->nrhome*p->ntaxbracket;
        // How to index into tax_data: 
        // itaxbracket*p->nrhome + irh 


        // Check that there are enough parameters
        // when they must correspond to the dimensions of the problem
        // if (p->ninc > ***) mexErrMsgTxt("Error in parameters: dimension mismatch in ***");
        if (p->ns > (p->n_phi_s+1)) mexErrMsgTxt("Error in parameters: dimension mismatch in phi_s");
        // if (p->nm > ***) mexErrMsgTxt("Error in parameters: dimension mismatch in ***");
        if (p->nkids > (p->n_phi_c+1)) mexErrMsgTxt("Error in parameters: dimension mismatch in phi_c");
        // if (p->ncouple > ***) mexErrMsgTxt("Error in parameters: dimension mismatch in ***");
        // if (p->nrwork > ***) mexErrMsgTxt("Error in parameters: dimension mismatch in ***");
        // if (p->nrhome > ***) mexErrMsgTxt("Error in parameters: dimension mismatch in ***");
        // if (p->nstatevar > ***) mexErrMsgTxt("Error in parameters: dimension mismatch in ***");
        // if (p->nchoicevar > ***) mexErrMsgTxt("Error in parameters: dimension mismatch in ***");

        if (p->ns > (p->n_Bpn_s+1)) mexErrMsgTxt("Error in parameters: dimension mismatch in Bpn_s");
        if (p->ns > (p->n_Bpk_s+1)) mexErrMsgTxt("Error in parameters: dimension mismatch in Bpk_s");
        if (p->ns > (p->n_gamma_s+1)) mexErrMsgTxt("Error in parameters: dimension mismatch in gamma_s");
        // p->n_gamma_c
        // p->n_phi_rh
        
        // FIXME: Write all checks for the dimensions of states and passed parameters

        // -- Calculated dependecies ------------------------------------------------------------------------------------------------------
        if (p->T - p->t0 < 0) mexErrMsgTxt("Error in f_parseparameters: T is smaller than minimum age t0");

        //--- Precompute state and choice specific incomes  ---------------------------------------------------------------- //
        ninc = f_get_ninc(p);
        p->inc = calloc(ninc, sizeof(double));

    } // end of year specific

}  // end of f_parseparameters

static void f_cleanparameters(param* p,const bool compute_derivatives, const int iy){
    // Free memory allocated in f_parseparameters.

    if (iy==INITIAL_PARSE)
    {   //free space allocated when parsing for all years
        free(p->years_data);
    }
    else
    {
        free(p->iparname);
        free(p->iparkey);
        free(p->isubpar);
        free(p->inc);
    }
}

int f_get_ninc(const param* p){
    // Retrieve the number of unique incomes in the state x decision dimensions.  

    int ninc;

    // Number of states relevant for income: N_age x N_schooling x 2 (Binary state, rw(t-1) == non_employment)
    // Number of choices relevant for income: N_work
    ninc = ((p->T - p->t0 + 1)*p->ns*2)*p->nrwork; 

    return ninc;
}

int f_inc_idx(const param* p, declSTATES, declOUTCOMES){
    // Calculate index of state-choice combination in precomputed income array, p->inc

    int idx, n_age, n_nonemp, nonemp;

    nonemp = (irwp == (p->nrwork - 1));
    n_nonemp = 2;
    n_age = p->T - p->t0 + 1;

    // Loop: irw -> it -> is -> nonemp
    idx = irw*(n_age*p->ns*n_nonemp) + it*(p->ns*n_nonemp) + is*n_nonemp + nonemp;

    return idx; 
}

#define IDEVF  (p->bl_t*it+ p->bl_rhp*irhp+ p->bl_rwp*irwp + p->bl_kids*ikids+ p->bl_couple*icouple+ p->bl_m*im+ p->bl_s*is+ p->bl_inc*iinc+ imove )

static int f_IDEVF(const param* p, declSTATES){
    //This returns the index in the EVF table, can be called as function or
    //as macro with current values of indices
return IDEVF; 
}

static int f_IDEVF_INV(const param* p, const int ievf, declpointSTATES){
    int idx;  
    idx=ievf; 

    *imove=idx%p->nmovec;
    idx=(idx-*imove)/p->nmovec;

    *iinc=idx%p->ninc;
    idx=(idx-*iinc)/p->ninc;

    *is=idx%p->ns;
    idx=(idx-*is)/p->ns;

    *im=idx%p->nm;
    idx=(idx-*im)/p->nm;

    *icouple=idx%p->ncouple;
    idx=(idx-*icouple)/p->ncouple;

    *ikids=idx%p->nkids;
    idx=(idx-*ikids)/p->nkids;

    *irwp=idx%p->nrwork;
    idx=(idx-*irwp)/p->nrwork;

    *irhp=idx%p->nrhome;
    idx=(idx-*irhp)/p->nrhome;

    *it=idx%p->lifespan;

    return idx;

}

// IDCHOICE: columns sorted by irh (slow) and irw (fast)
#define IDCHOICE (p->nrwork*irh + irw)

static int f_IDCHOICE(const param* p, declOUTCOMES){
    // This returns ichoice given OUTCOMES irw, irh
    return IDCHOICE; 
}

static void f_IDCHOICE_INV(const param* p, const int ichoice, declpointOUTCOMES){
    int idx;  
    idx=ichoice; 

    *irw=idx%p->nrwork;
    idx=(idx-*irw)/p->nrwork;

    *irh=idx%p->nrhome;

}

static int* yeardata_input(const param *p, const mxArray *input) {
    int i,j0,j1,ncols,nrows;
    int year,year0;
    double *years;
    int *year_row_idx;

    // Number of columns in year data input.
    ncols = (mwSize)mxGetN(input);

#ifdef FOOLPROOF
    // Check if the inputted choice data follows convention of a column for work and home respectively.
    ncols = (int) mxGetN(input); // Dimension of choice point vector
    if (ncols != 1) {
        printf("Number of columns are: %d \n",ncols );
        mexErrMsgTxt("Need one column for year data. ");
    }
#endif   

    // Number of observations in data input.
    nrows = (mwSize)mxGetM(input);

    // Pointer to input data.
    years = (double*) mxGetPr(input);

    // Allocate space for output structure.
    year_row_idx=(int*) calloc(2*p->n_years_data,sizeof(int));

#ifdef FOOLPROOF
    if (!year_row_idx) mexErrMsgTxt("yeardata_input: Failed to allocate memory for year_row_idx.");
#endif

    // ASSUME both p.years_data and years in the input data are sorted!
#ifdef FOOLPROOF
    // check 
    for (i = 0; i<p->n_years_data; i++)
    {
        year = p->years_data[i];
        if (i>1 && year<year0)
            mexErrMsgTxt("Years are not sorted in mp.years_data!");
        year0 = year;
    }
    for (i = 0; i<nrows; i++)
    {
        year = DBL2INT(years[i]);
        if (i>1 && year<year0)
            mexErrMsgTxt("Years are not sorted in the years data!");
        year0 = year;
    }
#endif

    j1 = 0; //j is running index in input data: j0=start, j1=end

    // Loop through the model years
    for (i = 0; i < p->n_years_data; i++)
    {   
        for (j0=j1; p->years_data[i] >  DBL2INT(years[j0]) && j0<nrows; j0++ );
        for (j1=j0; p->years_data[i] >= DBL2INT(years[j1]) && j1<nrows; j1++ );

        year_row_idx[2*i] = j0;   //start of year p->years_data[i]
        year_row_idx[2*i+1] = j1; //end of year p->years_data[i]
    }
    return year_row_idx;
}

static void f_adjust_statedata(const param *p, statedata *sp, const int ip_0, const int ip_N) {
    int ip;

    // Adjust data points to proper indices in C. 
    for (ip = ip_0; ip < ip_N; ++ip)
    {
        sp[ip].it -= p->t0;  // Age to index it.
        sp[ip].ikids = MIN(sp[ip].ikids, (p->nkids-1)); // If data points on kids exceed max defined in params, use the max (eg 4).
    }

#ifdef FOOLPROOF
    // If there is an index error, use ok flag to print out what line in data has an issue. 
    for (ip = ip_0; ip < ip_N; ++ip)
    {
        sp[ip].ok = 1;

        if (sp[ip].it       < 0 || sp[ip].it    > (p->T-p->t0)   || 
            sp[ip].irh      < 0 || sp[ip].irh   >= p->nrhome      ||
            sp[ip].irw      < 0 || sp[ip].irw   >= p->nrwork      ||  
            sp[ip].ikids    < 0 || sp[ip].ikids >= p->nkids       ||  
            sp[ip].icouple  < 0 || sp[ip].icouple >= p->ncouple   ||  
            sp[ip].im       < 0 || sp[ip].im    >= p->nm          || 
            sp[ip].is       < 0 || sp[ip].is    >= p->ns          || 
            sp[ip].iinc     < 0 || sp[ip].iinc  >= p->ninc        || 
            sp[ip].imove    < 0 || sp[ip].imove >= p->nmovec)
        {
            sp[ip].ok = 0;
        }

        if (!sp[ip].ok)
        {
            printf("There is an issue with data indices at line %d. \n", ip + 1);
        }


    // Check that data read in is in accordance with  specifications of model parsed to p.
    // Note that the parameters on state variables in p give the number of possible state. By indexing from 0, 
    // indices should therefore max be nstates-1.
        if (sp[ip].it < 0       || sp[ip].it    > (p->T - p->t0))   mexErrMsgTxt("Error in statedata_input: Age index is out of range.");
        if (sp[ip].irh < 0      || sp[ip].irh   >= p->nrhome)        mexErrMsgTxt("Error in statedata_input: Home region index is out of range");
        if (sp[ip].irw < 0      || sp[ip].irw   >= p->nrwork)        mexErrMsgTxt("Error in statedata_input: Work region index is out of range");
        if (sp[ip].ikids < 0    || sp[ip].ikids >= p->nkids)         mexErrMsgTxt("Error in statedata_input: Kids state index is out of range");
        if (sp[ip].icouple < 0  || sp[ip].icouple >= p->ncouple)     mexErrMsgTxt("Error in statedata_input: Couple index is out of range");
        if (sp[ip].im < 0       || sp[ip].im    >= p->nm)            mexErrMsgTxt("Error in statedata_input: Macro state index is out of range");
        if (sp[ip].is < 0       || sp[ip].is    >= p->ns)            mexErrMsgTxt("Error in statedata_input: Schooling state index is out of range");
        if (sp[ip].iinc < 0     || sp[ip].iinc  >= p->ninc)          mexErrMsgTxt("Error in statedata_input: Income type index is out of range");
        if (sp[ip].imove < 0    || sp[ip].imove >= p->nmovec)        mexErrMsgTxt("Error in statedata_input: Move type index is out of range");

    }
#endif

}


static statedata* statedata_input(const mxArray *input){
// This function takes as input a data matrix of observations passed from Matlab and 
// returns a struct of type statedata which holds a field for each state variable. The fields
// are arrays of observations of length nobs.

    int ip;
    mwSize nrows,ncols;
    double *inpoints; 
    statedata *sp;

    // Number of columns in data input.
    ncols = (mwSize)mxGetN(input);
    ncols = ncols;

    // Number of observations in data input.
    nrows = (mwSize)mxGetM(input);
   
    // Pointer to input data.
    inpoints = (double*) mxGetPr(input);

    // Allocate space for output structure.
    sp=(statedata*) calloc(nrows,sizeof(statedata));


#ifdef FOOLPROOF
    if (!sp) mexErrMsgTxt("statedata_input: Failed to allocate memory for input data.");
#endif

    // Loop over rows of data matrix and cast from double to int one by one (necessary). 
    // See definitions statedatastruct for order (same as STATES).
    for (ip = 0; ip < nrows; ++ip)
    {
        sp[ip].it       = DBL2INT(inpoints[ip + 0*nrows]); // Age
        sp[ip].irh      = DBL2INT(inpoints[ip + 1*nrows]); // Home region
        sp[ip].irw      = DBL2INT(inpoints[ip + 2*nrows]); // Work region
        sp[ip].ikids    = DBL2INT(inpoints[ip + 3*nrows]); // Kids state
        sp[ip].icouple  = DBL2INT(inpoints[ip + 4*nrows]); // Couple state
        sp[ip].im       = DBL2INT(inpoints[ip + 5*nrows]); // Macro state
        sp[ip].is       = DBL2INT(inpoints[ip + 6*nrows]); // School type
        sp[ip].iinc     = DBL2INT(inpoints[ip + 7*nrows]); // Income type
        sp[ip].imove    = DBL2INT(inpoints[ip + 8*nrows]); // Move type
    }

    return sp;

} 

static void f_adjust_choicedata(const param *p, choicedata *cp, const int ip_0, const int ip_N) {

    int ip, ok;

#ifdef FOOLPROOF
    for (ip = ip_0; ip < ip_N; ++ip)
    {
        ok=1;
 
        if (cp[ip].irh < 0 || cp[ip].irh >= (p->nrhome) || cp[ip].irw < 0 || cp[ip].irw >= p->nrwork){
            ok=0;
        } 

        if (!ok)
        {
            printf("There is an issue with data indices at line %d. \n", ip + 1);
        }

        // Check that data read in is in accordance with  specifications of model parsed to p.
        // Note that the parameters on state variables in p give the number of possible choices. By indexing from 0, 
        // indices should therefore max be nstates-1.
        if (cp[ip].irh < 0 || cp[ip].irh >= p->nrhome) mexErrMsgTxt("Error in choicedata input: Home region index is out of range.");
        if (cp[ip].irw < 0 || cp[ip].irw >= p->nrwork) mexErrMsgTxt("Error in choicedata input: Work region index is out of range");
    }
#endif
}

static choicedata* choicedata_input(const mxArray *input){
    // This function creates a structure containing choices of individuals. Used with statedata structure. 

#ifdef FOOLPROOF
    int ok, cp_col;
#endif   
    int ip;
    mwSize nrows,ncols;
    double *inpoints; 
    choicedata *cp;

#ifdef FOOLPROOF
    // Check if the inputted choice data follows convention of a column for work and home respectively.
    cp_col = (int) mxGetN(input); // Dimension of choice point vector
    if (cp_col != 2) {
        printf("Number of columns are: %d \n",cp_col );
        mexErrMsgTxt("Need two columns for choices, rhome-rwork. ");
    }
#endif   

    // Number of columns in data input.
    ncols = (mwSize)mxGetN(input);

    // Number of observations in data input.
    nrows = (mwSize)mxGetM(input);

    // Pointer to input data.
    inpoints = (double*) mxGetPr(input);

    // Allocate space for output structure.
    cp=(choicedata*) calloc(nrows,sizeof(choicedata));

#ifdef FOOLPROOF
    if (!cp) mexErrMsgTxt("choicedata_input: Failed to allocate memory for input data.");
#endif

    // Note: important that order follows that of choicedatastruct! work region should come in last column.
    for (ip = 0; ip < nrows; ++ip)
    {
        cp[ip].irh = DBL2INT(inpoints[ip]        );
        cp[ip].irw = DBL2INT(inpoints[ip + nrows]);
    }

    return cp;

}

static void save_statedata(const param* p, const statedata *sp, const int i, const int N, double *sim_states_out){
    sim_states_out[ARGit       *N + i]=sp[i].it+p[0].t0;
    sim_states_out[ARGirhp     *N + i]=sp[i].irh   ;
    sim_states_out[ARGirwp     *N + i]=sp[i].irw   ;
    sim_states_out[ARGikids    *N + i]=sp[i].ikids  ;
    sim_states_out[ARGicouple  *N + i]=sp[i].icouple;
    sim_states_out[ARGim       *N + i]=sp[i].im     ;
    sim_states_out[ARGis       *N + i]=sp[i].is     ;
    sim_states_out[ARGiinc     *N + i]=sp[i].iinc   ;
    sim_states_out[ARGimove    *N + i]=sp[i].imove  ;
}

static void save_choicedata(const param* p, const choicedata *cp, const int i, const int N, double *sim_choice_out){
    sim_choice_out[ARGirh  *N + i]=cp[i].irh;
    sim_choice_out[ARGirw  *N + i]=cp[i].irw;
}


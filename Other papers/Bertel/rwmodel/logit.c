/*********************************************************************************************************
 vfs_allchoices calculates the value function vfs for all choices at STATES
 *********************************************************************************************************
 
 SYNTAX             :   vfs_allchoices(p, declSTATES, *evf, *vfs)
 
 INPUTS
 ---------------------------------------------------------------------------------------------------------
    p               :   Parameter structure that contains model parameters, grids, quadrature points, etc. 
                        (see help in rwmodel.m)

    declSTATES      :   State that the value function should be calculated at.

    evf             :   Pointer to next period's values.

    vfs             :   Pointer to array of value functions for all choices of ip at state sp that is to be filled out. 
**********************************************************************************************
 Know issues (FIXME):
**********************************************************************************************/

static double logsum(const int n,const double* x,const double lambda){
    //calcualtes the logsum of the n first values of x
    int i;
    double res,mm;
    mm=x[0];
    for (i=1;i<n;i++)
    {
        if (fabs(x[i]-INFEASIBLECHOICEVF)<1e-10) continue;//skip infeasible choices
        if (x[i]>mm || fabs(mm-INFEASIBLECHOICEVF)<1e-10) mm=x[i];
    }
#ifdef FOOLPROOF
    if (fabs(mm-INFEASIBLECHOICEVF)<1e-10) mexErrMsgTxt("Logsum for empty choiceset can not be computed!");
#endif
    res=0.0;
    for (i=0;i<n;i++)
    {
        if (fabs(x[i]-INFEASIBLECHOICEVF)<1e-10) continue;//skip infeasible choices
        res+=exp((x[i]-mm)/lambda);
    }
    return mm+lambda*log(res);
}

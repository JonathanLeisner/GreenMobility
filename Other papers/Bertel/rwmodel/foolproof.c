#define NOPRINT // PRINTIDEFV, NOPRINT, PRINTCHOICE

static void check_index(const param* p){
    // ckeck of f_ID...
    int ievf, idout, STATES, OUTCOMES, ichoice;
    int ievf_chk, idout_chk;
    int it_chk, irhp_chk, irwp_chk, ikids_chk, icouple_chk, im_chk, is_chk, iinc_chk, imove_chk, ichoice_chk;

    #ifdef PRINTCHOICE    
    for (ichoice=0;(ichoice < (p[0].nchoice)); ichoice++) {  // loop over all ichoice

        f_IDCHOICE_INV(p, ichoice, pointOUTCOMES);
        ichoice_chk=f_IDCHOICE(p, OUTCOMES);         
        printf("f_IDCHOICE_INV ichoice=%d, ichoice_chk=%d, irh=%d, irw=%d \n", ichoice, ichoice_chk, irh, irw);
        if (ichoice_chk!=ichoice)  {
            printf("f_IDCHOICE=%d, IDCHOICE=%d, ichoice=%d, IDCHOICE-ichoice=%d\n", ichoice_chk, IDCHOICE, ichoice, IDCHOICE-ichoice);
            mexErrMsgTxt("f_IDCHOICE_INV and f_IDCHOICE do not match");     
        }
    }
    #endif
    
    #ifdef PRINTIDEFV // FIXME: should be PRINIDEVF? 
    for (ievf=0;(ievf<p[0].nstate); ievf++) {// loop over all idevf (equal to numer of states)
        
        f_IDEVF_INV(p, ievf, &it_chk, &irhp, &irwp, &ikids, &icouple, &im, &is, &iinc, &imove);

        printf("f_IDEVF_INV: ievf=%d, it=%d, irhp=%d, irwp=%d, ikids=%d, icouple=%d, im=%d, is=%d, iinc=%d, imove=%d\n",
            ievf,    it,    irhp,    irwp,    ikids,    icouple,    im,    is,    iinc,    imove) ; 
    }
    #endif

    for (it = 0; (it < (p[0].T-p[0].t0+1)); it++) { // loop over time periods

        for (ievf = p[0].bl_t*it; (ievf < p[0].bl_t*(it+1)); ievf++) { // loop over states

            // ckeck consistency of IDEVF index functions
            f_IDEVF_INV(p, ievf, &it_chk, &irhp, &irwp, &ikids, &icouple, &im, &is, &iinc, &imove);

            ievf_chk=f_IDEVF(p, it_chk, irhp, irwp, ikids, icouple, im, is, iinc, imove);

            if (ievf_chk!=ievf)  {
                printf("f_IDEVF=%d, IDEVF=%d, IDEVF-ievf=%d\n", ievf_chk, IDEVF, IDEVF-ievf);
                mexErrMsgTxt("f_IDEVF_INV and f_IDEVF do not match");     
            }
        }     
    }
    printf("Checked IDEVF index for all %d states and %d choices, p[0].nstate*p[0].nchoice = %d \n", p[0].nstate, p[0].nchoice, p[0].nstate*p[0].nchoice);

}

#include <stdio.h>
#include <stdlib.h>


#include "cdvode-internal.h"


static CDVODE_ODE_RETVAL istate(CDVODE_ISTATE_OUT istat)
{
    switch (istat) {
    case NOTHING_DONE:
        return 2;
    case SUCCESS_DONE:
        return 0;
    case MXSTEPS_EXCEEDED:
    case TO_MUCH_ACCURACY:
    case ILLEGAL_INPUT:
    case ERROR_TEST_FAILURES:
    case CONVERGENCE_FAILURES:
    case ZERO_ERR_TOLERANCE:
        return -1;
    default :
        return -20;
    }

}

static CDVODE_METHOD_FLAG mf(CDVODE_METH meth, CDVODE_JSV jsv, CDVODE_MITER miter)
{
    return jsv * (10 * meth + miter);
}



static int lrw(CDVODE_METHOD_FLAG mf, int neq, int ml, int mu)
{
    
    switch (mf) {
    case 10:
        return 20 * 16 * neq;
    case 11:
    case 12:
        return 22 + 16 * neq + 2 * neq * neq;
    case -11:
    case -12:
        return 22 + 16 * neq + neq * neq;
    case 13:
        return 22 + 17 * neq;
    case 14:
    case 15 :
        return 22 + 18 * neq + (3*ml + 2*mu)*neq;
    case -14:
    case -15 :
        return 22 + 17 * neq + (2*ml + mu)*neq;
    case 20:
        return 20 + 9 * neq;
    case 21:
    case 22:
        return 22 + 9 * neq+ 2 * neq * neq;
    case -21:
    case -22:
        return 22 + 9 * neq + neq * neq;
    case 23 :
        return 22 + 10 * neq;
    case 24:
    case 25:
        return 22 + 11 * neq + (3*ml + 2*mu)*neq;
    case -24:
    case -25:
        return 22 + 10 * neq + (2*ml + mu)*neq;

    };
    return 0;
}

static int liw(CDVODE_METHOD_FLAG mf, int neq)
{
    switch (mf) {
    case 10:
    case 13:
    case 20:
    case 23:
        return 30;
    default:
        return 30 + neq;
    }
}

static CDVODE_ODE_RETVAL valid_mf(CDVODE_METH meth, CDVODE_JSV jsv, CDVODE_MITER miter)
{
    CDVODE_METHOD_FLAG mf_ = mf(meth, jsv, miter);
    switch (mf_) {
    case 10:
    case 11:
    case 12:
    case 13:
    case 14:
    case 15 :

    case 20:
    case 21:
    case 22:
    case 23 :
    case 24:
    case 25:        

    case -11:
    case -12:

    case -14:
    case -15 :

    case -21:
    case -22:

    case -24:
    case -25:
        return SUCCESS;
    default :
        return ERROR;
    };
    
}

static cdvode_params* create_basic_cdvode_params(int neq, cdvode_ode_func f_func, CDVODE_METH meth, CDVODE_JSV jsv, CDVODE_MITER miter)
{
    const double rtol = 0.0;
    const double atol = 1e-12;
    cdvode_params* dsp = calloc(1, sizeof(cdvode_params));
    if(dsp == NULL) {
        return NULL;
    }
    dsp->ode_func = DVODE;
    dsp->f_func = f_func;
    dsp->neq    = neq;
    dsp->itol   = ALL_SCALAR;
    dsp->rtol.rtol_val = rtol;
    dsp->atol.atol_val = atol;
    dsp->itask = NORMAL;
    dsp->istate.istate_in = FIRST_CALL;
    dsp->iopt   = NO_OPTIONAL_INPUTS;
    dsp->lrw    = lrw(mf(meth, jsv, miter), neq, 0, 0);
    dsp->rwork  = (double*) calloc(dsp->lrw, sizeof(double));
    dsp->liw    = liw(mf(meth, jsv, miter), neq);
    dsp->iwork  = (int *)  calloc(dsp->liw, sizeof(int));
    dsp->jac    = NULL;
    dsp->method_params.jsv   = jsv;
    dsp->method_params.miter = miter;
    dsp->method_params.meth  = meth;
    dsp->mf     = mf(meth, jsv, miter);
    dsp->ipar   = NULL;
    dsp->rpar   = NULL;
/*
    fprintf(stderr, "%d %d %d\n", dsp->mf , dsp->lrw, dsp->liw);
*/  
    return dsp;
}

static void free_params(cdvode_params* dsp)
{
    free(dsp->iwork);
    free(dsp->rwork);
    free(dsp);
}


static CDVODE_ISTATE_OUT dvode_w(cdvode_params *dlsodap, double tnext, double *t, double *q) 
{

    //fprintf(stderr, "calling dlsoda_\n");
    dvode(dlsodap->f_func,
          &dlsodap->neq,
          q,
          t,
          &tnext,
          &dlsodap->itol,
          &dlsodap->rtol.rtol_val,
          &dlsodap->atol.atol_val,
          (const int*) &dlsodap->itask,
          &dlsodap->istate.istate,
          (const int*) &dlsodap->iopt,
          dlsodap->rwork,
          &dlsodap->lrw,
          dlsodap->iwork,
          &dlsodap->liw,
          dlsodap->jac,
          (const int*)&dlsodap->mf,
          dlsodap->rpar,
          dlsodap->ipar);

    return dlsodap->istate.istate_out;
}

CDVODE_ODE_RETVAL dvode_basic(double* stack, double* q, cdvode_ode_func f_func,int neq, double t0, double tf, double dt, CDVODE_METH meth)
{
    double t;
    CDVODE_ODE_RETVAL ode_ret = SUCCESS;
    int index = 0;
    if (valid_mf(meth, SAVE_COPY, NO_JACOBIAN) != SUCCESS) {
        fprintf(stderr, "invalid method flag...\n");
        return ERROR;
    }
    cdvode_params* dlsop = create_basic_cdvode_params(neq, f_func, meth, SAVE_COPY, NO_JACOBIAN);
    t = t0;
    stack = write_to_stack(stack, neq, &index, t, q);
    while(t < tf){    
        CDVODE_ISTATE_OUT ret = dvode_w(dlsop, t + dt, &t, q);
        ode_ret = istate(ret);
        if (ode_ret < 0) {
            break;
        }
        stack = write_to_stack(stack, neq, &index, t, q);
    }
    free_params(dlsop);
    return ode_ret;
}

void dvode(void (*f)(const int *neq, const double *t, const double *y, double *ydot, double* rpar, int* ipar),
           const int *neq,
           double *y,
           double *t,
           const double *tout,
           const int *itol,
           const double *rtol,
           const double *atol,
           const int *itask,
           int *istate,
           const int *iopt,
           double *rwork,
           const int *lrw,
           int *iwork,
           const int *liw,
           void (*jac)(const int *neq, const double *t, const double *y, const int *ml, const int *mu, double *pd, const int *nrowpd, double* rpar, int* ipar),
           const int *mf,
           double* rpar,
           int* ipar)
{
    dvode_(f,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf,rpar,ipar);
}


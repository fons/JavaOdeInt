
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "czvode-internal.h"


static CZVODE_ODE_RETVAL istate(CZVODE_ISTATE_OUT istat)
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

static CZVODE_METHOD_FLAG mf(CZVODE_METH meth, CZVODE_JSV jsv, CZVODE_MITER miter)
{
    return jsv * (10 * meth + miter);
}



static int lzw(CZVODE_METHOD_FLAG mf, int neq, int ml, int mu)
{
    
    switch (mf) {
    case 10:
        return 15 * neq;
    case 11:
    case 12:
        return 15 * neq + 2 * neq * neq;
    case -11:
    case -12:
        return 15 * neq + neq * neq;
    case 13:
        return 16 * neq;
    case 14:
    case 15 :
        return 17 * neq + (3*ml + 2*mu)*neq;
    case -14:
    case -15 :
        return 16 * neq + (2*ml + mu)*neq;
    case 20:
        return 8 * neq;
    case 21:
    case 22:
        return 8 * neq+ 2 * neq * neq;
    case -21:
    case -22:
        return 8 * neq + neq * neq;
    case 23 :
        return 9 * neq;
    case 24:
    case 25:
        return 10 * neq + (3*ml + 2*mu)*neq;
    case -24:
    case -25:
        return 9 * neq + (2*ml + mu)*neq;

    };
    return 0;
}

static int lrw(CZVODE_METHOD_FLAG mf, int neq, int ml, int mu)
{
    
    return 20 + neq;
}

static int liw(CZVODE_METHOD_FLAG mf, int neq)
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

static CZVODE_ODE_RETVAL valid_mf(CZVODE_METH meth, CZVODE_JSV jsv, CZVODE_MITER miter)
{
    CZVODE_METHOD_FLAG mf_ = mf(meth, jsv, miter);
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

static czvode_params* create_basic_czvode_params(int neq, czvode_ode_func f_func, CZVODE_METH meth, CZVODE_JSV jsv, CZVODE_MITER miter)
{
    const double rtol = 0.0;
    const double atol = 1e-12;
    czvode_params* dsp = calloc(1, sizeof(czvode_params));
    if(dsp == NULL) {
        return NULL;
    }
    dsp->ode_func = ZVODE;
    dsp->f_func = f_func;
    dsp->neq    = neq;
    dsp->itol   = ALL_SCALAR;
    dsp->rtol.rtol_val = rtol;
    dsp->atol.atol_val = atol;
    dsp->itask = NORMAL;
    dsp->istate.istate_in = FIRST_CALL;
    dsp->iopt   = NO_OPTIONAL_INPUTS;

    dsp->lzw    = lzw(mf(meth, jsv, miter), neq, 0, 0);
    dsp->zwork  = (ZVODE_COMPLEX*) calloc(dsp->lzw, sizeof(ZVODE_COMPLEX));

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
    fprintf(stderr, "\n===>%d %d %d %d\n", dsp->mf , dsp->lzw, dsp->lrw, dsp->liw);
*/  
    return dsp;
}

static void free_params(czvode_params* dsp)
{
    free(dsp->zwork);
    free(dsp->iwork);
    free(dsp->rwork);
    free(dsp);
}

static CZVODE_ISTATE_OUT zvode_w(czvode_params *dlsodap, double tnext, double *t, ZVODE_COMPLEX *q)
{
    /* fprintf(stderr, "\n\ncalling zvode_..\n");/
    fprintf(stderr, "%lf %lf \n", q[0].dr, q[1].dr);
    */
    assert(dlsodap->f_func != NULL);
    
    zvode(dlsodap->f_func,
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
          dlsodap->zwork,
          &dlsodap->lzw,
          dlsodap->rwork,
          &dlsodap->lrw,
          dlsodap->iwork,
          &dlsodap->liw,
          dlsodap->jac,
          (const int*)&dlsodap->mf,
          dlsodap->rpar,
          dlsodap->ipar);
    /* fprintf(stderr, "DONE calling zvode_..\n");*/
    return dlsodap->istate.istate_out;
}

CZVODE_ODE_RETVAL zvode_basic(double* stack, ZVODE_COMPLEX* q, czvode_ode_func f_func,int neq, double t0, double tf, double dt, CZVODE_METH meth)
{
    double t = t0;
    CZVODE_ODE_RETVAL ode_ret = SUCCESS;

    int index = 0;
    if (valid_mf(meth, SAVE_COPY, NO_JACOBIAN) != SUCCESS) {
        fprintf(stderr, "invalid method flag...\n");
        return ERROR;
    }
 
    czvode_params* dlsop = create_basic_czvode_params(neq, f_func, meth, SAVE_COPY, NO_JACOBIAN);
    stack = write_to_stack(stack, neq, &index, t, q);
    while(t < tf){    
        CZVODE_ISTATE_OUT ret = zvode_w(dlsop, t + dt, &t, q);
        ode_ret = istate(ret);
        if (ode_ret < 0) {
            return ode_ret;
        }
        stack = write_to_stack(stack, neq, &index, t, q);
    }
    free_params(dlsop);
    return ode_ret;

}

void zvode(void (*f)(const int *neq, const double *t, const ZVODE_COMPLEX *y, ZVODE_COMPLEX *ydot, ZVODE_COMPLEX* rpar, int* ipar),
                   const int *neq,
                   ZVODE_COMPLEX *y,
                   double *t,
                   const double *tout,
                   const int *itol,
                   const double *rtol,
                   const double *atol,
                   const int *itask,
                   int *istate,
                   const int *iopt,
                   ZVODE_COMPLEX* zwork,
                   const int *lzw,
                   double *rwork,
                   const int *lrw,
                   int *iwork,
                   const int *liw,
                   void (*jac)(const int *neq, const double *t, const ZVODE_COMPLEX *y, const int *ml, const int *mu, ZVODE_COMPLEX *pd, const int *nrowpd, ZVODE_COMPLEX* rpar, int* ipar),
                   const int *mf,
                   ZVODE_COMPLEX* rpar,
                   int* ipar)
{
    
    zvode_(f, neq, y, t, tout,itol, rtol, atol, itask, istate, iopt, zwork, lzw, rwork, lrw, iwork, liw, jac, mf, rpar, ipar);


}

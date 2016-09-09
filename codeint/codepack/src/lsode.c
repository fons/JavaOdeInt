#include <stdio.h>
#include <stdlib.h>

#include "codepack-internal.h"


static int lrw(CODEPACK_METHOD_FLAG mf, int neq, int ml, int mu)
{
    switch (mf) {
    case ADAMS_BASIC :
        return 20 + 16 * neq;
    case ADAMS_USER_FULL_JAC:
    case ADAMS_INTERNAL_FULL_JAC:
        return 22 + 16 * neq + neq*neq;

    case ADAMS_INTERNAL_DIAG_JAC:
        return 22 + 17 * neq;

    case ADAMS_USER_BAND_JAC:

    case ADAMS_INTERNAL_BAND_JAC:
        if (mu > 0 && mu > 0) {
            return 22 + 17 * neq + (2*ml + mu)*neq;
        }
        else {
            return 0;
        }
    case BDF_BASIC:
        return 20 + 9 * neq;
        
    case BDF_USER_FULL_JAC:
    case BDF_INTERNAL_FULL_JAC:
        return 22 + 9 * neq + neq * neq;
        
    case BDF_INTERNAL_DIAG_JAC:
        return 22 + 10 * neq;
    case BDF_USER_BAND_JAC:
    case BDF_INTERNAL_BAND_JAC:
        if (mu > 0 && mu > 0) {
            return 22 + 10 * neq + (2 *ml + mu)* neq;
        }
        else {
            return 0;
        }
    default :
        return 0;
    };
}

static int liw(CODEPACK_METHOD_FLAG mf, int neq)
{

    switch (mf) {
        
    case ADAMS_BASIC:
    case ADAMS_INTERNAL_DIAG_JAC:
    case BDF_BASIC:
    case BDF_INTERNAL_DIAG_JAC:
        return 20;
        
    case ADAMS_USER_FULL_JAC:
    case ADAMS_INTERNAL_FULL_JAC:
    case ADAMS_USER_BAND_JAC:
    case ADAMS_INTERNAL_BAND_JAC:
    case BDF_USER_FULL_JAC:
    case BDF_INTERNAL_FULL_JAC:
    case BDF_USER_BAND_JAC:
    case BDF_INTERNAL_BAND_JAC:
        return 20 + neq;
    default:
        return 0;
    }
}

static lsod_params* create_basic_lsode_params(int neq, codepack_ode_func f_func, CODEPACK_METHOD_FLAG mf)
{
    const double rtol = 0.0;
    const double atol = 1e-12;
    lsod_params* dsp = calloc(1, sizeof(lsod_params));
    if(dsp == NULL) {
        return NULL;
    }
    dsp->ode_func = DLSODE;
    dsp->f_func = f_func;
    dsp->neq    = neq;
    dsp->itol   = ALL_SCALAR;
    dsp->rtol.rtol_val = rtol;
    dsp->atol.atol_val = atol;
    dsp->itask = NORMAL;
    dsp->istate.istate_in = FIRST_CALL;
    dsp->iopt   = NO_OPTIONAL_INPUTS;

    dsp->lrw    = lrw(mf, neq, -1, -1);

    dsp->rwork     = (double*) calloc(dsp->lrw, sizeof(double));
    dsp->liw       = liw(mf,neq);
    dsp->iwork     = (int *)  calloc(dsp->liw, sizeof(int));
    dsp->jac.jac1  = NULL;
    dsp->mf        = mf;

    return dsp;
}

static void free_params(lsod_params* dsp)
{
    free(dsp->iwork);
    free(dsp->rwork);
    free(dsp);
}

static CODEPACK_ISTATE_OUT call_lsode(lsod_params *dlsodap, double tnext, double *t, double *q)
{

    dlsode_(dlsodap->f_func,
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
            dlsodap->jac.jac1,
            (const int*)&dlsodap->mf);

    return dlsodap->istate.istate_out;
}

CODEPACK_ODE_RETVAL lsode_basic(double* stack, double* q, codepack_ode_func f_func,int neq, double t0, double tf, double dt, CODEPACK_METHOD_FLAG mf)
{
    double t = 0.0;
    CODEPACK_ODE_RETVAL ode_ret = SUCCESS;
    int index = 0;
    lsod_params* dlsop = create_basic_lsode_params(neq, f_func, mf);
    stack = write_to_stack(stack, neq, &index, t, q);
    t = t0;
    while(t < tf){    
        CODEPACK_ISTATE_OUT ret = call_lsode(dlsop, t + dt, &t, q);
        ode_ret = istate(ret);
        if (ode_ret < 0) {
            break;
        }
        stack = write_to_stack(stack, neq, &index, t, q);
    }
    free_params(dlsop);
    return ode_ret;   
}

void dlsode(void (*f)(const int *neq, const double *t, const double *y, double *ydot),
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
            void (*jac)(const int *neq, const double *t, const double *y, const int *ml, const int *mu, double *pd, const int *nrowpd),
            const int *mf)
{
    dlsode_(f,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf);
}

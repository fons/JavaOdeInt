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

lsod_params* create_basic_lsode_params(int neq, codepack_ode_func f_func, CODEPACK_METHOD_FLAG mf)
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

CODEPACK_ISTATE_OUT lsode(double t, double *t0, double *q, lsod_params *dlsodap)
{

    //fprintf(stderr, "calling dlsoda_\n");
    dlsode_(dlsodap->f_func,
            &dlsodap->neq,
            q,
            t0,
            &t,
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
    double t;
    CODEPACK_ODE_RETVAL ode_ret;
    int index = 0;
    lsod_params* dlsop = create_basic_lsode_params(neq, f_func, mf);
    t = t0;
    while(t < tf){    
        CODEPACK_ISTATE_OUT ret = lsode(t + dt, &t, q, dlsop);
        ode_ret = istate(ret);
        if (ode_ret < 0) {
            return ode_ret;
        }
        stack = write_to_stack(stack, neq, &index, (t+dt), q);
    }
    
    return ode_ret;   
}

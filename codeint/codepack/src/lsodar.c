
#include <stdio.h>
#include <stdlib.h>


#include "codepack-internal.h"

/*
 * neq is static during the execution 
 */
static int lrw(int neq, int ng, CODEPACK_JAC_TYPE jc, int ml, int mu)
{
    int lrs = 0;
    int lrn = 20 + 16 * neq + 3 * ng;
    switch (jc){
    case  USER_PROVIDED:
    case INTERNAL:
        lrs = 22 + 9 * neq + neq * neq + 3 * ng;
        break;
    case USER_PROVIDED_BANDED:
    case INTERNAL_BANDED:
        lrs = 22 + 10 * neq + (2 * ml + mu)*neq + 3 * ng;
        break;
    default :
        lrs = 0;
    }
    return MAX(lrs, lrn);
}

lsod_params* create_basic_lsodar_params(int neq, codepack_ode_func f_func)
{
    const double rtol = 0.0;
    const double atol = 1e-12;
    lsod_params* dsp = calloc(1, sizeof(lsod_params));
    if(dsp == NULL) {
        return NULL;
    }
    dsp->ode_func = DLSODAR;
    dsp->f_func = f_func;
    dsp->neq    = neq;
    dsp->itol   = ALL_SCALAR;
    dsp->rtol.rtol_val = rtol;
    dsp->atol.atol_val = atol;
    dsp->itask = NORMAL;
    dsp->istate.istate_in = FIRST_CALL;
    dsp->iopt   = NO_OPTIONAL_INPUTS;
    dsp->lrw    = lrw(neq, 0, INTERNAL, 0, 0);
    dsp->rwork  = (double*) calloc(dsp->lrw, sizeof(double));
    dsp->liw    = neq + 20;
    dsp->iwork  = (int *)  calloc(dsp->liw, sizeof(int));
    dsp->jac.jac1    = NULL;
    dsp->jt     = INTERNAL;
    dsp->g      = NULL;
    dsp->ng     = 0;
    dsp->jroot  = (int *)  calloc(dsp->ng, sizeof(int));

    return dsp;
}

CODEPACK_ISTATE_OUT lsodar(double t, double *t0, double *q, lsod_params *dlsodap)
{

    //fprintf(stderr, "calling dlsoda_\n");
    dlsodar_(dlsodap->f_func,
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
             (const int*)&dlsodap->jt,
             dlsodap->g,
             (const int *)&dlsodap->ng,
             dlsodap->jroot
        );

    return dlsodap->istate.istate_out;
}

CODEPACK_ODE_RETVAL lsodar_basic(double* stack, double* q, codepack_ode_func f_func,int neq, double t0, double tf, double dt)
{
    double t;
    CODEPACK_ODE_RETVAL ode_ret;
    int index = 0;
    lsod_params* dlsop = create_basic_lsodar_params(neq, f_func);
    t = t0;
    while(t < tf){    
        CODEPACK_ISTATE_OUT ret = lsodar(t + dt, &t, q, dlsop);
        ode_ret = istate(ret);
        if (ode_ret < 0) {
            return ode_ret;
        }
        stack = write_to_stack(stack, neq, &index, (t+dt), q);
    }
    
    return ode_ret;
}

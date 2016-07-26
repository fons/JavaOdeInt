
#include <stdio.h>
#include <stdlib.h>


#include "codepack-internal.h"


lsod_params* create_basic_lsoda_params(int neq, codepack_ode_func f_func)
{
    const double rtol = 0.0;
    const double atol = 1e-12;
    lsod_params* dsp = calloc(1, sizeof(lsod_params));
    if(dsp == NULL) {
        return NULL;
    }
    dsp->ode_func = DLSODA;
    dsp->f_func = f_func;
    dsp->neq    = neq;
    dsp->itol   = ALL_SCALAR;
    dsp->rtol.rtol_val = rtol;
    dsp->atol.atol_val = atol;
    dsp->itask = NORMAL;
    dsp->istate.istate_in = FIRST_CALL;
    dsp->iopt   = NO_OPTIONAL_INPUTS;
    dsp->lrw    = 22 + neq * MAX(16, neq+9);
    dsp->rwork  = (double*) calloc(dsp->lrw, sizeof(double));
    dsp->liw    = neq + 20;
    dsp->iwork  = (int *)  calloc(dsp->liw, sizeof(int));
    dsp->jac.jac1    = NULL;
    dsp->jt     = INTERNAL;
    

    return dsp;
}

CODEPACK_ISTATE_OUT lsoda(double t, double *t0, double *q, lsod_params *dlsodap)
{

    //fprintf(stderr, "calling dlsoda_\n");
    dlsoda_(dlsodap->f_func,
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
            (const int*)&dlsodap->jt);

    return dlsodap->istate.istate_out;
}

CODEPACK_ODE_RETVAL lsoda_basic(double* stack, double* q, codepack_ode_func f_func,int neq, double t0, double tf, double dt)
{
    double t;
    CODEPACK_ODE_RETVAL ode_ret;
    int index = 0;
    lsod_params* dlsop = create_basic_lsoda_params(neq, f_func);
    t = t0;
    while(t < tf){    
        CODEPACK_ISTATE_OUT ret = lsoda(t + dt, &t, q, dlsop);
        ode_ret = istate(ret);
        if (ode_ret < 0) {
            return ode_ret;
        }
        stack = write_to_stack(stack, neq, &index, (t+dt), q);
    }
    
    return ode_ret;
}

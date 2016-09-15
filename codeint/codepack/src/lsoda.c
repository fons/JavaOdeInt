
#include <stdio.h>
#include <stdlib.h>


#include "codepack-internal.h"


static lsod_params* create_basic_lsoda_params(int neq, codepack_ode_func f_func)
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

static void free_params(lsod_params* dsp)
{
    free(dsp->iwork);
    free(dsp->rwork);
    free(dsp);
}

static CODEPACK_ISTATE_OUT lsoda(lsod_params *dlsodap, double tnext, double *t, double *q)
{
    dlsoda_(dlsodap->f_func,
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
           (const int*)&dlsodap->jt);

    return dlsodap->istate.istate_out;
}

CODEPACK_ODE_RETVAL lsoda_basic(double* stack, double* q, codepack_ode_func f_func,int neq, double t0, double tf, double dt)
{
    int max_retries = 5;
    int retry       = 0;
    double t;
    double tnext;
    CODEPACK_ODE_RETVAL ode_ret = SUCCESS;
    int index = 0;
    lsod_params* dlsop = create_basic_lsoda_params(neq, f_func);
    t = t0;
    stack = write_to_stack(stack, neq, &index, t, q);
    while(t < tf){
        retry = 0;
        tnext = t + dt;
        CODEPACK_ISTATE_OUT return_code = SUCCESS_DONE;
        do {
            return_code = lsoda(dlsop, tnext, &t, q);
            if (return_code == MAX_STEPS_EXCEEDED) {
                retry++;
                if (retry >= max_retries) {
                    break;
                }
                dlsop->iopt      = OPTIONAL_INPUTS;
                dlsop->iwork[5] += 2000;
                dlsop->istate.istate_in = NEXT_CALL_WITH_CHANGES;
                fprintf(stderr, "increased max steps to %d for retry %d \n", dlsop->iwork[5], retry);
                t = tnext - dt;
            }
            else {
                retry = 0;
                break;
            }
        } while (retry > 0);
        ode_ret = istate(return_code);
        if (ode_ret < 0) {
            break;
        }
        stack = write_to_stack(stack, neq, &index, t, q);
    }
    free_params(dlsop);
    return ode_ret;
}


void dlsoda(void (*f)(const int *neq, const double *t, const double *y, double *ydot),
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
            const int *jt)
{ 
    dlsoda_(f,neq,y,t,tout, itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac, jt);   
}

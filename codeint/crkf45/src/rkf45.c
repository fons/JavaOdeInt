#include <stdlib.h>
#include <stdio.h>
#include "rkf45-internal.h"

static rkf45_params* create_basic_rkf45_params(int neq, rkf45_ode_func f_func)
{
    double rtol = 0.0000000001;
    double atol = 0.0000000001;
    rkf45_params* dsp = calloc(1, sizeof(rkf45_params));
    if(dsp == NULL) {
        return NULL;
    }
    dsp->f_func = f_func;
    dsp->neq    = neq;
    dsp->rtol   = rtol;
    dsp->atol   = atol;
    dsp->iflag.iflag = NORMAL;
    dsp->lrw    = 6 * neq + 3;
    dsp->rwork  = (double*) calloc(dsp->lrw, sizeof(double));
    dsp->liw    = 5;
    dsp->iwork  = (int *)  calloc(dsp->liw, sizeof(int));
    return dsp;
}

static void free_params(rkf45_params* dsp)
{
    free(dsp->iwork);
    free(dsp->rwork);
    free(dsp);
}

static RKF45_RETVAL rkf45_w(rkf45_params* params, double tnext, double* t, double* y)
{
    params->iflag.iflag = NORMAL;
    r8_rkf45_(
        params->f_func,
        &params->neq,
        y,
        t,
        &tnext,
        &params->rtol,
        &params->atol,
        &params->iflag.iflag,
        params->rwork,
        params->iwork
        );
    return params->iflag.retval;
}

RKF45_RETVAL rkf45_basic(double* stack, double* q,  rkf45_ode_func f_func,int neq, double xstart, double xfinal, double deltax)
{

    RKF45_RETVAL retval;
    double x = xstart;
    int index = 0;
    rkf45_params* dlsop = create_basic_rkf45_params(neq, f_func);
    stack = write_to_stack(stack, neq, &index, x, q);
    while (x < xfinal) {
        retval = rkf45_w(dlsop, x + deltax, &x, q);
        if (retval != SUCCESS) {
            break;
        }
        stack = write_to_stack(stack, neq, &index, x, q);
    }
    free_params(dlsop);
    return retval;
}

void rkf45 ( void (*f)(double*, double* , double*),
             int* neqn,
             double* y,
             double* t,
             double* tout,
             double* relerr,
             double* abserr, 
             int* iflag,
             double* work,
             int* iwork )
{
    r8_rkf45_(f,neqn, y,t,tout,relerr,abserr,iflag,work,iwork);
}

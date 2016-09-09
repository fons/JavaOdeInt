#include <stdio.h>
#include <stdlib.h>

#include "dop853-internal.h"

static cdop853_params* create_basic_cdop853_params(int neq, dop853_ode_func f_func)
{
    double rtol = 0.000000001;
    double atol = 0.000000001;
    cdop853_params* dsp = calloc(1, sizeof(cdop853_params));
    if(dsp == NULL) {
        return NULL;
    }
    dsp->f_func = f_func;
    dsp->neq    = neq;
    dsp->itol   = ALL_SCALAR;
    dsp->rtol.rtol_val = rtol;
    dsp->atol.atol_val = atol;
    dsp->lrw    = 11 * neq + 21;
    dsp->rwork  = (double*) calloc(dsp->lrw, sizeof(double));
    /*
     * seems to work a bit better with 0 as an initial condition
     */
    dsp->rwork[6] = 0.0000001;
    dsp->liw    = 21;
    dsp->iwork  = (int *)  calloc(dsp->liw, sizeof(int));
    dsp->solout    = NULL;
    dsp->iout      = NEVER_CALLED;
    dsp->rpar      = NULL;
    dsp->ipar      = NULL;

    return dsp;
}

static void free_params(cdop853_params* dsp)
{
    free(dsp->iwork);
    free(dsp->rwork);
    free(dsp);
}

static DOP853_RETVAL dop853_w(cdop853_params *params, double tnext, double* t, double* y)
{
    dop853_(&params->neq,
            params->f_func,
            t,
            y,
            &tnext,
            &params->rtol.rtol_val,
            &params->atol.atol_val,
            (int*)&params->itol,
            params->solout,
            (int*)&params->iout,
            params->rwork,
            &params->lrw,
            params->iwork,
            &params->liw,
            params->rpar,
            params->ipar,
            (int *) &params->retval
        );
    return params->retval;
}

DOP853_RETVAL dop853_basic(double* stack, double* y, dop853_ode_func f_func, int n, double xstart, double xfinal, double deltax)
{
    
    DOP853_RETVAL retval = SUCCESS;
    double x = xstart;
    int index = 0;
    cdop853_params* dlsop = create_basic_cdop853_params(n, f_func);
    stack = write_to_stack(stack, n, &index, x, y);
    while (x < xfinal) {
        retval = dop853_w(dlsop, x + deltax, &x, y);
        if (retval != SUCCESS) {
            break;
        }
        stack = write_to_stack(stack, n, &index, x, y);
    }
    free_params(dlsop);
    return retval;
}


void dop853(
    int* n,
    void (*fcn) (const int* n, const double *x, const double *y, double *f, double* rpar, int* ipar),
    double* x,
    double* y,
    double* xend,
    double* rtol,
    double* atol,
    const int* itol,
    void (*solout)(int* nr, double* xold, double* x, double* y, int* n, double* con, int* icomp, int* nd, double* rpar, int* ipar, int* irtrn),
    int* iout,
    double* work,
    int* lwork,
    int* iwork,
    int* liwork,
    double* rpar,
    int* ipar,
    int* idid
    )
{
    dop853_(n, fcn, x,y,xend,rtol, atol, itol, solout, iout, work, lwork, iwork, liwork, rpar, ipar, idid);
}


#include <stdio.h>
#include <stdlib.h>

#include "dopri5-internal.h"

static cdopri5_params* create_basic_cdopri5_params(int neq, dopri5_ode_func f_func)
{
    double rtol = 0.000000001;
    double atol = 0.000000001;
    cdopri5_params* dsp = calloc(1, sizeof(cdopri5_params));
    if(dsp == NULL) {
        return NULL;
    }
    dsp->f_func = f_func;
    dsp->neq    = neq;
    dsp->itol   = ALL_SCALAR;
    dsp->rtol.rtol_val = rtol;
    dsp->atol.atol_val = atol;
    dsp->lrw    = 8 * neq + 21;
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

static void free_params(cdopri5_params* dsp)
{
    free(dsp->iwork);
    free(dsp->rwork);
    free(dsp);
}

static DOPRI5_RETVAL dopri5_w(cdopri5_params *params, double tnext, double* t, double* y)
{
    dopri5_(&params->neq,
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

DOPRI5_RETVAL dopri5_basic(double* stack, double* y, dopri5_ode_func f_func, int n, double xstart, double xfinal, double deltax)
{
    
    DOPRI5_RETVAL retval;
    double x = xstart;
    int index = 0;
    stack = write_to_stack(stack, n, &index, x, y);
    cdopri5_params* dlsop = create_basic_cdopri5_params(n, f_func);
    while (x < xfinal) {
        retval = dopri5_w(dlsop, x + deltax, &x, y);
        if (retval != SUCCESS) {
            break;
        }
        stack = write_to_stack(stack, n, &index, x, y);
    }
    free_params(dlsop);
    return retval;
}


void dopri5(
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
    dopri5_(n, fcn, x,y,xend,rtol, atol, itol, solout, iout, work, lwork, iwork, liwork, rpar, ipar, idid);
}


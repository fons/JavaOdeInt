#include <stdio.h>
#include <stdlib.h>

#include "radau5-internal.h"

static radau5_params* create_basic_radau5_params(int neq, radau5_ode_func f_func)
{
    double rtol = 0.0000000001;
    double atol = 0.0000000001;
    double initial_step = 0.0;
    radau5_params* dsp = calloc(1, sizeof(radau5_params));
    if(dsp == NULL) {
        return NULL;
    }
    dsp->neq            = neq;
    dsp->f_func         = f_func;
    dsp->jac_func       = NULL;
    dsp->ijac           = INTERNAL;
    dsp->mljac          = neq;
    dsp->mujac          = 0;
    dsp->mass_func      = NULL;
    dsp->imas           = IDENTITY_MATRIX;
    dsp->mlmas          = neq;
    dsp->mumas          = 0;
    dsp->solout         = NULL;
    dsp->iout           = NEVER_CALLED;
    dsp->itol           = ALL_SCALAR;
    dsp->rtol.rtol_val  = rtol;
    dsp->atol.atol_val  = atol;
    dsp->lrw = 10 *(4 * neq * neq + 12*neq + 20);/*default; full jacobian; no mass matrix*/
    dsp->rwork = (double *) calloc(dsp->lrw, sizeof(double));
    dsp->liw   = 10 * (3 * neq + 20);
    dsp->iwork = (int*) calloc(dsp->liw, sizeof(int));
    dsp->rpar  = NULL;
    dsp->ipar  = NULL;
    dsp->h     = initial_step;
    return dsp;
}

static void free_params(radau5_params* dsp)
{
    free(dsp->iwork);
    free(dsp->rwork);
}


static RADAU5_RETVAL radau5_w(radau5_params* dsp, double tnext, double* t, double* y)
{
    radau5_(&dsp->neq,
            dsp->f_func,
            t,
            y,
            &tnext,
            &dsp->h,
            &dsp->rtol.rtol_val,
            &dsp->atol.atol_val,
            (int *) &dsp->itol,
            dsp->jac_func,
            (int *)&dsp->ijac,
            &dsp->mljac,
            &dsp->mujac,
            dsp->mass_func,
            (int *)&dsp->imas,
            &dsp->mlmas,
            &dsp->mumas,
            dsp->solout,
            (int*) &dsp->iout,
            dsp->rwork,
            &dsp->lrw,
            dsp->iwork,
            &dsp->liw,
            dsp->rpar,
            dsp->ipar,
            &dsp->retval);
    
    return dsp->retval;
}

RADAU5_RETVAL radau5_basic(double* stack, double* y, radau5_ode_func f_func, int n, double xstart, double xfinal, double deltax)
{

    RADAU5_RETVAL retval;
    double x = xstart;
    int index = 0;
    radau5_params* dlsop = create_basic_radau5_params(n, f_func);
    stack = write_to_stack(stack, n, &index, x, y);

    while (x < xfinal) {
        retval = radau5_w(dlsop, x + deltax, &x, y);
        if (retval != SUCCESS) {
            break;
        }
        stack = write_to_stack(stack, n, &index, x, y);
    }
    free_params(dlsop);
    return retval;
}

void radau5(
    int* n,
    void (*fcn)(int* , double* , double* , double* , double* , int* ),
    double* x,
    double* y,
    double* xend,
    double* h,
    double* rtol,
    double* atol,
    int* itol,
    void (*jac)(int* , double* , double* , double** , int* , double* , int* ),
    int* ijac,
    int* mljac,
    int* mujac,
    void (*mas)(int* , double** ,int* ,double* , int* ),
    int* imas,
    int* mlmas,
    int* mumas,
    void (*solout)(int* , double* , double* , double* , double* , int*, int*, double*, int* , int* ),
    int* iout,
    double* work,
    int* lwork,
    int* iwork,
    int* liwork,
    double* rpar,
    int* ipar,
    int *idid
    )
{
    radau5_(n,fcn,x,y,xend,h,rtol,atol,itol,jac,ijac,mljac,mujac,mas,imas,mlmas,mumas,solout,iout,work,lwork,iwork,liwork,rpar,ipar, idid); 
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gnicodes-internal.h"


static gnicodes_params* create_basic_gni_lmm2_params(int neq, gnicodes_ode_func f_func, int nstep, GNI_LMM2_METH meth)
{
    
    gnicodes_params* dsp = calloc(1, sizeof(gnicodes_params));
    if(dsp == NULL) {
        return NULL;
    }
    dsp->neq            = neq;
    dsp->f_func         = f_func;
    dsp->solfix         = NULL;
    dsp->iout           = NEVER_CALLED;
    dsp->nstep          = nstep;
    dsp->meth.lmm2_meth = meth;
    dsp->lr             = 10;
    dsp->rpar           = (double *) calloc(dsp->lr, sizeof(double));
    dsp->li             = 10;
    dsp->ipar           = (int *) calloc(dsp->li, sizeof(int));

    return dsp;
}

static void gni_lmm2_w(gnicodes_params* dsp, double xend, double* x, double* p, double* q)
{
    gni_lmm2_(&dsp->neq,
              dsp->f_func,
              &dsp->nstep,  
              x,
              p,
              q,
              &xend,
              (int *) &dsp->meth.lmm2_meth,
              dsp->solfix,
              (int *)&dsp->iout,
              dsp->rpar,
              dsp->ipar);
}




void gni_lmm2_basic(double* stack, double* p, double* q, gnicodes_ode_func f_func, int n, double xstart, double xfinal, double deltax, GNI_LMM2_METH meth)
{


    double x = xstart;
    int index = 0;
    int steps = abs((int) ((xfinal - xstart)/deltax));
    
    gnicodes_params* dlsop = create_basic_gni_lmm2_params(n, f_func, 30, meth);
    stack = write_pq_to_stack(stack, n, &index, x, p, q);
    while (x < xfinal) {
        gni_lmm2_w(dlsop, x + deltax, &x, p, q);
        stack = write_pq_to_stack(stack, n, &index, x, q, p);
    }

}

void gni_lmm2(
    int* n,
    void(*f)(int* n, double* x, double* q, double* f, double* rpar, int* ipar),
    int* nstep,
    double* x,
    double* p,
    double* q,
    double* xend,
    int* meth,
    void (*s) (int* nr, double* xold, double* x, double* p, double* q, int* n, int* irtrn, double* rpar, int* ipar),
    int* iout,
    double* rpar,
    int* ipar
    )
{
    gni_lmm2_(n,f,nstep, x,p,q,xend,meth,s,iout,rpar,ipar);
}

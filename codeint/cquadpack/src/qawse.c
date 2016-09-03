#include <stdlib.h>
#include <stdio.h>

#include "../include/cquadpack.h"
#include "cquadpack-internal.h"

QUADPACK_ERRNO qawse_basic(cquadpack_ode_func func, double a, double b, double alfa, double beta, QUADPACK_LOG_WEIGHT_FUNCTION integr, double* result)
{
    QUADPACK_INPUT  inp;
    QUADPACK_OUTPUT out;

    default_quadpack_input(&inp);
    default_quadpack_output(inp, &out);

    dqawse_(func,
            &a,
            &b,
            &alfa,
            &beta,
            (int *)&integr,
            &inp.epsabs,
            &inp.epsrel,
            &inp.limit,
            result,
            &out.abserr,
            &out.neval,
            &out.ier,
            out.alist,
            out.blist,
            out.rlist,
            out.elist,
            out.iord,
            &out.last);

    free_quadpack_output(&out);

    return out.ier;
}

void dqawse(double (*f)(double *y),
            double* a,
            double* b,
            double* alfa,
            double* beta,
            int* integr,
            double* epsabs,
            double* epsrel,
            int* limit,
            double* result,
            double* abserr,
            int* neval,
            int* ier,
            double* alist,
            double* blist,
            double* rlist,
            double* elist,
            int* iord,
            int* last)

{
    dqawse_(f,a,b,alfa,beta,integr,epsabs,epsrel,limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last);
}


#include <stdlib.h>
#include <stdio.h>

#include "../include/cquadpack.h"
#include "cquadpack-internal.h"

QUADPACK_ERRNO qagse_basic(cquadpack_ode_func func, double start, double end, double* result)
{
    QUADPACK_INPUT  inp;
    QUADPACK_OUTPUT out;

    default_quadpack_input(&inp);
    default_quadpack_output(inp, &out);

    dqagse_(func,
            &start,
            &end,
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

void dqagse(double (*f)(double *y),
            double* a,
            double* b,
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
    dqagse_(f,a,b,epsabs,epsrel,limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last);
}

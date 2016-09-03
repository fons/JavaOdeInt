#include <stdlib.h>
#include <stdio.h>

#include "../include/cquadpack.h"
#include "cquadpack-internal.h"

QUADPACK_ERRNO qawce_basic(cquadpack_ode_func func, double start, double end, double weight, double* result)
{
    QUADPACK_INPUT  inp;
    QUADPACK_OUTPUT out;

    default_quadpack_input(&inp);
    default_quadpack_output(inp, &out);

    dqawce_(func,
            &start,
            &end,
            &weight,
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
void dqawce(double (*f)(double *y),
                    double* a,
                    double* b,
                    double* c,
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
    dqawce_(f,a,b,c,epsabs,epsrel,limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last);
}

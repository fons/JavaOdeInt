#include <stdlib.h>
#include <stdio.h>

#include "../include/cquadpack.h"
#include "cquadpack-internal.h"

QUADPACK_ERRNO qawoe_basic(cquadpack_ode_func func, double a, double b, double omega, QUADPACK_TRIG_WEIGHT_FUNCTION integr, double* result)
{
    QUADPACK_INPUT  inp;
    QUADPACK_OUTPUT out;

    default_quadpack_input(&inp);
    default_quadpack_output(inp, &out);


    dqawoe_(func,
            &a,
            &b,
            &omega,
            (int *)&integr,
            &inp.epsabs,
            &inp.epsrel,
            &inp.limit,
            &inp.icall,
            &inp.maxpl,
            result,
            &out.abserr,
            &out.neval,
            &out.ier,
            &out.last,
            out.alist,
            out.blist,
            out.rlist,
            out.elist,
            out.iord,
            out.nnlog,
            &out.momcom,
            out.chebmo
            );
    free_quadpack_output(&out);

    return out.ier;
}


void dqawoe(double (*f)(double *y),
            double* a,
            double* b,
            double* omega,
            int* integr,
            double* epsabs,
            double* epsrel,
            int* limit,
            int* icall,
            int* maxpl,
            double* result,
            double* abserr,
            int* neval,
            int* ier,
            int* last,
            double* alist,
            double* blist,
            double* rlist,
            double* elist,
            int* iord,
            int* nnlog,
            int* momcom,
            double* chebmo)
            {
            dqawoe_(f,a,b,omega,integr,epsabs,epsrel,limit,icall,maxpl,result,abserr,neval,ier,last,alist,blist,rlist,elist,iord,nnlog,momcom,chebmo);
            }

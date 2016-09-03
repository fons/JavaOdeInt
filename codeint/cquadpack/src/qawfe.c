#include <stdlib.h>
#include <stdio.h>

#include "../include/cquadpack.h"
#include "cquadpack-internal.h"

QUADPACK_ERRNO qawfe_basic(cquadpack_ode_func func, double a, double omega, QUADPACK_TRIG_WEIGHT_FUNCTION integr, double* result)
{
    QUADPACK_INPUT  inp;
    QUADPACK_OUTPUT out;

    default_quadpack_input(&inp);
    default_quadpack_output(inp, &out);
    
    dqawfe_(func,
            &a,
            &omega,
            (int *)&integr,
            &inp.epsabs,
            &inp.limlst,
            &inp.limit,
            &inp.maxpl,
            result,
            &out.abserr,
            &out.neval,
            &out.ier,
            out.rslst,
            out.erlst,
            out.ierlst,
            &out.lst,
            out.alist,
            out.blist,
            out.rlist,
            out.elist,
            out.iord,
            out.nnlog,
            out.chebmo
            );
    free_quadpack_output(&out);

    return out.ier;
}
void dqawfe(double (*f)(double *y),
            double* a,
            double* omega,
            int* integr,
            double* epsabs,
            int* limlst,
            int* limit,
            int* maxpl,
            double* result,
            double* abserr,
            int* neval,
            int* ier,
            double* rslst,
            double* erlst,
            int*    ierlst,
            int* lst,
            double* alist,
            double* blist,
            double* rlist,
            double* elist,
            int* iord,
            int* nnlog,
            double* chebmo
    )
{
    dqawfe_(f,a,omega,integr,epsabs,limlst,limit,maxpl,result,abserr,neval,ier,
            rslst,erlst,ierlst,lst,alist,blist,rlist,elist,iord,nnlog,chebmo);
}

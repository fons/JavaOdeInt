#include <stdlib.h>
#include <stdio.h>

#include "../include/cquadpack.h"
#include "cquadpack-internal.h"

QUADPACK_ERRNO qagpe_basic(cquadpack_ode_func func, double start, double end, int npts2, double* points, double* result)
{
    QUADPACK_INPUT  inp;
    QUADPACK_OUTPUT out;

    default_quadpack_input(&inp);
    default_quadpack_output(inp, &out);
    points_quadpack_output(npts2, &out);
    
    dqagpe_(func,
            &start,
            &end,
            &npts2,
            points,
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
            out.pts,
            out.level,
            out.ndin,
            out.iord,
            &out.last);

    free_quadpack_output(&out);

    return out.ier;
}

void dqagpe(double (*f)(double *y),
            double* a,
            double* b,
            int* npts2,
            double* points,
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
            double* pts,
            int* level,
            int*  ndin,
            int* iord,
            int* last)
{
    dqagpe_(f,a,b,npts2,points,epsabs,epsrel,limit,result,abserr,neval,ier,alist,blist,rlist,elist,pts,level,ndin,iord,last);
}

#include "daspk-internal.h"


void daspk(
           void(*res)(double*, double*, double*, double*, double* , int* , double*, int*),
           int* neq,
           double* t,
           double* y,
           double* yprime,
           double* tout,
           int* info,
           double* rtol,
           double* atol,
           int* idid,
           double* rwork,
           int* lrw,
           int* iwork,
           int* liw,
           double* rpar,
           int* ipar,
           void (*jac)(double*, double*, double*, double*, double*),
           void (*psol)(int*, double*, double*, double*, double*, double* , double*, double*, double*, int*, double* )
    )
{
    ddaspk_(res, neq, t,y,yprime,tout,info,rtol,atol,idid,rwork,lrw,iwork,liw,rpar,ipar,jac,psol);
}


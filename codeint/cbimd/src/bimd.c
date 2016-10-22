#include "bimd-internal.h"


void bimd(
    int* m,
    void (*fcn)(int* , double* , double* , double* , int* , double* , int* ),
    double* t0,
    double* tend,
    double* y0,
    double* h,
    double* rtol,
    double* atol,
    int* itol,
    void (*jac)(int* , double* , double* , double* , int* , int*, double* , int* ),
    int* ijac,
    int* mljac,
    int* mujac,
    void (*mas)(int* , double* ,int* , int*,double* , int* ),
    int* imas,
    int* mlmas,
    int* mumas,
    void (*solout)(int* , int*, int*, double* , double* , double* , double* , double*, double*, int* , int* ),
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
    bimd_(m,fcn,t0,tend, y0,h,rtol,atol,itol,jac,ijac,mljac,mujac,mas,imas,mlmas,mumas,solout,iout,work,lwork,iwork,liwork,rpar,ipar, idid);
}

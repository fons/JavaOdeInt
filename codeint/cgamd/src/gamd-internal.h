#ifndef H__GAMD_INTERNAL__H
#define H__GAMD_INTERNAL__H
#include "../include/cgamd.h"

void gamd_(
    int* r,
    void (*fcn)(int* , double* , double* , double* , int* , double* , int* ),
    double* t0,
    double* y0,
    double* tend,
    double* h,
    double* rtol,
    double* atol,
    int* itol,
    void (*jac)(int* , double* , double* , double* , int* , double* , int* ),
    int* ijac,
    int* mljac,
    int* mujac,
    void (*mas)(int* , double* ,int* ,double* , int* ),
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
    );

#endif

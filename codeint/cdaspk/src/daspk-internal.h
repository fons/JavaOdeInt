#ifndef H__DASPK__H
#define H__DASPK__H
#include "../include/cdaspk.h"
extern void ddaspk_(
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
    void (*psol)(int*, double*, double*, double*, double*, double* , double*, double*, double*, int*, double* ));

#endif

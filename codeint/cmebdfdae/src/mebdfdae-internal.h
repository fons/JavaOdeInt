#ifndef H__MEBDFDAE_H
#define H__MEBDFDAE_H
#include "../include/cmebdfdae.h"
extern void mebdf_(
    int* n,
    double* t0,
    double* h0,
    double* y0,
    double* tout,
    double* tend,
    int* mf,
    int* idid,
    int* lout,
    int* lwork,
    double* work,
    int* liwork,
    int* iwork,
    int* mbnd,
    int* maxder,
    int* itol,    
    double* rtol,
    double* atol,
    double* rpar,
    int* ipar,
    void (*f)(int*, double*, double*, double*, int*, double* ),
    void (*pderv)(double* , double* , double* , double* , int*, double* , int* ),
    void (*mas)(int*, double*, double* ,int*, double*),
    int *ierr
    );

#endif

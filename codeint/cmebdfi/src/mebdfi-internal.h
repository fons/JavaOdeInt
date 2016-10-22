#ifndef H__MEBDFI__H
#define H__MEBDFI__H
#include "../include/cmebdfi.h"
extern void mebdfi_(
    int* n,
    double* t0,
    double* h0,
    double* y0,
    double* yprime,
    double* tout,
    double* tend,
    int* mf,
    int* idid,
    int* lout,
    int* mbnd,
    int* maxder,
    int* itol,    
    double* rtol,
    double* atol,
    double* rpar,
    int* ipar,
    void (*pderv)(double* , double* , double* , double* , int*, double* , int* ),
    void (*resid)(int*, double*, double* , double*, int*, double* , int* ),
    int *ierr
    );

#endif

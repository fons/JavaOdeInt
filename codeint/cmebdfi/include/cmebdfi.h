#ifndef H__CMEBDFI__H
#define H__CMEBDFI__H
void mebdfi(
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

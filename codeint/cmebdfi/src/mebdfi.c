#include "mebdfi-internal.h"

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
            void (*pderv)(double* , /*t*/
                          double* , /*y*/
                          double* , /*pd*/
                          int*,    /*n*/
                          double* , /*yprime*/
                          int*,     /*mbnd*/
                          double* , /*con*/
                          int*    , /*ipar*/
                          double* , /*rpar*/
                          int* /*ier*/
                          ),
            void (*resid)(int*, double*, double* , double*, double*, int*, double* , int* ),
            int *ierr
            )
{
    mebdfi_(n,t0,h0, y0, yprime, tout,tend,mf,idid,lout,lwork, work, liwork, iwork, mbnd,maxder,itol,rtol, atol, rpar, ipar, pderv, resid,ierr);
}

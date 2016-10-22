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
            )
{
    mebdfi_(n,t0,h0, y0, yprime, tout,tend,mf,idid,lout,mbnd,maxder,itol,rtol, atol, rpar, ipar, pderv, resid,ierr);
}

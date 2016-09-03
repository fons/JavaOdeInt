#ifndef H__DOPRI5_INTERNAL__H
#define H__DOPRI5_INTERNAL__H
#include "../../codeintdeps/include/stack.h"
#include "../include/cdopri5.h"

typedef struct cdopri5_params_s {
    
    dopri5_ode_func f_func;
    int neq;
    DOPRI5_ITOLERANCE itol;
    union {
        double *rtol_vec;
        double rtol_val;
    } rtol;

    union {
        double *atol_vec;
        double atol_val;
    } atol;

    double* rwork;
    int*    iwork;
    int lrw;
    int liw;
    DOPRI5_IOUT iout;
    dopri5_solout solout;
    DOPRI5_RETVAL retval;
    double *rpar;
    int *ipar;
} cdopri5_params;


extern void dopri5_(
    int* n,
    void (*fcn) (const int* n, const double *x, const double *y, double *f, double* rpar, int* ipar),
    double* x,
    double* y,
    double* xend,
    double* rtol,
    double* atol,
    const int* itol,
    void (*solout)(int* nr, double* xold, double* x, double* y, int* n, double* con, int* icomp, int* nd, double* rpar, int* ipar, int* irtrn),
    int* iout,
    double* work,
    int* lwork,
    int* iwork,
    int* liwork,
    double* rpar,
    int* ipar,
    int* idid
    );

double* write_to_stack(double* stack, int neq, int* index, double t_new, double* q);

#endif

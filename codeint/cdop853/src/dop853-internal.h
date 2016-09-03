#ifndef H__DOP853_INTERNAL__H
#define H__DOP853_INTERNAL__H
#include "../../codeintdeps/include/stack.h"
#include "../include/cdop853.h"

typedef struct cdop853_params_s {
    
    dop853_ode_func f_func;
    int neq;
    DOP853_ITOLERANCE itol;
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
    DOP853_IOUT iout;
    dop853_solout solout;
    DOP853_RETVAL retval;
    double *rpar;
    int *ipar;
} cdop853_params;


extern void dop853_(
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

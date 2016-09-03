#ifndef H__RADAU5_INTERNAL__H
#define H__RADAU5_INTERNAL__H

#include "../../codeintdeps/include/stack.h"
#include "../include/cradau5.h"

typedef struct radau5_params_s {
    radau5_ode_func f_func;
    radau5_jacobian jac_func;
    radau5_mass     mass_func;
    radau5_solout   solout;
    int neq;
    RADAU5_ITOLERANCE itol;
    union {
        double *rtol_vec;
        double rtol_val;
    } rtol;

    union {
        double *atol_vec;
        double atol_val;
    } atol;
    int mljac;
    int mujac;
    double h;
    double* rwork;
    int*    iwork;
    int lrw;
    int liw;
    RADAU5_JACOBIAN ijac;
    RADAU5_MASS_MATRIX imas;
    int mlmas;
    int mumas;
    RADAU5_IOUT iout;
    RADAU5_RETVAL retval;
    double *rpar;
    int *ipar;
} radau5_params;

void radau5_(
    int* n,
    void (*fcn)(int* , double* , double* , double* , double* , int* ),
    double* x,
    double* y,
    double* xend,
    double* h,
    double* rtol,
    double* atol,
    int* itol,
    void (*jac)(int* , double* , double* , double** , int* , double* , int* ),
    int* ijac,
    int* mljac,
    int* mujac,
    void (*mas)(int* , double** ,int* ,double* , int* ),
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
    int* idid
    );


#endif

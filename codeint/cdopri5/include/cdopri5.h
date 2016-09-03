#ifndef H__DOPRI5__H
#define H__DOPRI5__H

/*
  n : diremsion
  x : independent variable
  y : vector
  f : derivative of y
 */

typedef void (*dopri5_ode_func)(const int* n, const double *x, const double *y, double *f, double* rpar, int* ipar);

typedef void (*dopri5_solout)(int* nr, double* xold, double* x, double* y, int* n, double* con, int* icomp, int* nd, double* rpar, int* ipar, int* irtrn);


typedef enum dopri5_idid_e {SUCCESS=1,
                            SUCCESS_INTR=2,
                            INPUT_INCONSISTENT=-1,
                            NMAX_TOO_SMALL=-2,
                            STEP_TOO_SMALL=-3,
                            STIFF_PROBLEM=-4} DOPRI5_RETVAL;

typedef enum dopri5_itol_e { ALL_SCALAR = 1, ALL_ARRAY=2} DOPRI5_ITOLERANCE;
typedef enum dopri5_iout_e { NEVER_CALLED=0, OUTPUT=2, DENSE_OUTPUT=3 } DOPRI5_IOUT; 

DOPRI5_RETVAL dopri5_basic(double* stack, double* y, dopri5_ode_func f_func, int n, double xstart, double xend, double deltax);

void dopri5(
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

#endif

#ifndef H__DOP853__H
#define H__DOP853__H

/*
  n : diremsion
  x : independent variable
  y : vector
  f : derivative of y
 */

typedef void (*dop853_ode_func)(const int* n, const double *x, const double *y, double *f, double* rpar, int* ipar);

typedef void (*dop853_solout)(int* nr, double* xold, double* x, double* y, int* n, double* con, int* icomp, int* nd, double* rpar, int* ipar, int* irtrn);


typedef enum dop853_idid_e {SUCCESS=1,
                            SUCCESS_INTR=2,
                            INPUT_INCONSISTENT=-1,
                            NMAX_TOO_SMALL=-2,
                            STEP_TOO_SMALL=-3,
                            STIFF_PROBLEM=-4} DOP853_RETVAL;

typedef enum dop853_itol_e { ALL_SCALAR = 1, ALL_ARRAY=2} DOP853_ITOLERANCE;
typedef enum dop853_iout_e { NEVER_CALLED=0, OUTPUT=2, DENSE_OUTPUT=3 } DOP853_IOUT; 

DOP853_RETVAL dop853_basic(double* stack, double* y, dop853_ode_func f_func, int n, double xstart, double xend, double deltax);

void dop853(
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

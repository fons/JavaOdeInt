#ifndef H__CRADUA5__H
#define H__CRADUA5__H

typedef void (*radau5_ode_func)(int* n, double* x, double* y, double* f, double* rpar, int* ipar);
typedef void (*radau5_jacobian)(int* n, double* x, double* y, double** dfy, int* ldfy, double* rpar, int* ipar);
typedef void (*radau5_mass)(int* n, double** am,int* lmas,double* rpar, int* ipar);
typedef void (*radau5_solout)(int* nr, double* xold, double* x, double* y, double* cont, int* lrc, int* n, double* rpar, int* ipar, int* irtrn); 

typedef enum radua5_itol_e { ALL_SCALAR = 0, ALL_ARRAY=1} RADAU5_ITOLERANCE;
typedef enum radau5_jacobian_e {INTERNAL=0, JAC_USER_PROVIDED=1 } RADAU5_JACOBIAN;
typedef enum radau5_mass_matrix_e {IDENTITY_MATRIX=0, MASS_USER_PROVIDED=1} RADAU5_MASS_MATRIX;
typedef enum radau5_iout_e { NEVER_CALLED=0, OUTPUT=1} RADAU5_IOUT;

typedef enum radau5_idid_e {SUCCESS=1,
                            SUCCESS_INTR=2,
                            INPUT_INCONSISTENT=-1,
                            NMAX_TOO_SMALL=-2,
                            STEP_TOO_SMALL=-3,
                            STIFF_PROBLEM=-4} RADAU5_RETVAL;

void radau5(
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

RADAU5_RETVAL radau5_basic(double* stack, double* y, radau5_ode_func f_func, int n, double xstart, double xfinal, double deltax);

#endif

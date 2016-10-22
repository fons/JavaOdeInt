#ifndef H__CGAMD__H
#define H__CGAMD__H


typedef enum gamd_itol_e { ALL_SCALAR = 0, ALL_ARRAY=1} GAMD_ITOLERANCE;
typedef enum gamd_jacobian_e {INTERNAL=0, JAC_USER_PROVIDED=1 } GAMD_JACOBIAN;
typedef enum gamd_mass_matrix_e {IDENTITY_MATRIX=0, MASS_USER_PROVIDED=1} GAMD_MASS_MATRIX;
typedef enum gamd_iout_e { NEVER_CALLED=0, OUTPUT=1} GAMD_IOUT;

typedef enum gamd_idid_e {SUCCESS=1,
                          INPUT_INCONSISTENT=-1,
                          NMAX_TOO_SMALL=-2,
                          STEP_TOO_SMALL=-3,
                          STIFF_PROBLEM=-4,
                          IERR_NEGATIVE=-5} GAMD_RETVAL;

void gamd(
    int* r,
    void (*fcn)(int* , double* , double* , double* , int* , double* , int* ),
    double* t0,
    double* y0,
    double* tend,
    double* h,
    double* rtol,
    double* atol,
    int* itol,
    void (*jac)(int* , double* , double* , double* , int* , double* , int* ),
    int* ijac,
    int* mljac,
    int* mujac,
    void (*mas)(int* , double* ,int* ,double* , int* ),
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
    int *idid
    );

#endif

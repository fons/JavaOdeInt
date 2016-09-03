#ifndef H__CDVODE__H
#define H__CDVODE__H

/*
 * function pointers compatible with the odepack 
 */


/*
 * ode function 
 */
typedef void (*cdvode_ode_func)(const int *neq, const double *t_, const double *y, double *ydot, double* rpar, int* ipar);

/*
 * jacobian 
 */
typedef void (*cdvode_jac_func) (const int *neq, const double *t_, const double *y, const int *ml, const int *mu, double *pd, const int *nrowpd, double* rpar, int* ipar);

typedef enum cdvode_itol_e { ALL_SCALAR = 1, ATOL_ARRAY=2, RTOL_ARRAY=3, ALL_ARRAY=4} CDVODE_ITOLERANCE;
typedef enum code_itask_e {NORMAL = 1, ONESTEP_ONLY = 2, STOP_AT_FIRST_MESH = 3, NORMAL_TCRIT = 4, ONESTEP_TCRIT = 5} CDVODE_ITASK;
typedef enum cdvode_istate_in_e {FIRST_CALL = 1, NEXT_CALL = 2, NEXT_CALL_WITH_CHANGES = 3} CDVODE_ISTATE_IN;
typedef enum cdvode_istate_out_e { NOTHING_DONE           = 1,
                                  SUCCESS_DONE           = 2,
                                  MXSTEPS_EXCEEDED       = -1,
                                  TO_MUCH_ACCURACY       = -2,
                                  ILLEGAL_INPUT          = -3 ,
                                  ERROR_TEST_FAILURES    = -4,
                                  CONVERGENCE_FAILURES   = -5,
                                  ZERO_ERR_TOLERANCE     = -6} CDVODE_ISTATE_OUT;

typedef enum cdvode_iopt_e {NO_OPTIONAL_INPUTS = 0, OPTIONAL_INPUTS = 1} CDVODE_OPTIONAL_INPUT_FLAG;

typedef enum cdvode_jacobian_saving_strategy_e {SAVE_COPY=1 , NO_COPY=-1} CDVODE_JSV;

typedef enum cdvode_method_e {ADAMS=1,
                             BDF=2 } CDVODE_METH;

typedef enum cdvode_iteration_method_e {NO_JACOBIAN=0,
                                       EXTERNAL_FULL_JACOBIAN = 1,
                                       INTERNAL_FULL_JACOBIAN = 2,
                                       INTERNAL_DIAG_JACOBIAN = 3,
                                       EXTERNAL_BAND_JACOBIAN = 4,
                                       INTERNAL_BAND_JACOBIAN = 5} CDVODE_MITER;


typedef enum cdvode_ode_err_e {UNCHANGED = 2, SUCCESS = 0, ERROR = -1, UNKNOWN_ERROR = -20} CDVODE_ODE_RETVAL;


CDVODE_ODE_RETVAL dvode_basic(double* stack, double* q, cdvode_ode_func f_func,int neq, double t0, double tf, double dt, CDVODE_METH meth);

void dvode(void (*f)(const int *neq, const double *t, const double *y, double *ydot, double* rpar, int* ipar),
           const int *neq,
           double *y,
           double *t,
           const double *tout,
           const int *itol,
           const double *rtol,
           const double *atol,
           const int *itask,
           int *istate,
           const int *iopt,
           double *rwork,
           const int *lrw,
           int *iwork,
           const int *liw,
           void (*jac)(const int *neq, const double *t, const double *y, const int *ml, const int *mu, double *pd, const int *nrowpd, double* rpar, int* ipar),
           const int *mf,
           double* rpar,
           int* ipar);


#endif

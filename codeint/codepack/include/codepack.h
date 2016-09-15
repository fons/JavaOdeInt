#ifndef H__CODEPACK__H
#define H__CODEPACK__H

/*
 * function pointers compatible with the odepack 
 */

/*
 * ode function 
 */
typedef void (*codepack_ode_func)(const int *neq, const double *t_, const double *y, double *ydot);

/*
 * jacobian 
 */
typedef void (*codepack_jac_func_1) (const int *neq, const double *t_, const double *y, const int *ml, const int *mu, double *pd, const int *nrowpd);
typedef void (*codepack_jac_func_2)(const int *neq, const double *t, const double *y, const int *j, const double *ian, double *jan, double *pdj);
/*
 * constraint function
 */
typedef void (*codepack_g_func) (const int *neq, const double *t_, const double *y, const int *ng, double *gout);

typedef enum codepack_itol_e { ALL_SCALAR = 1, ATOL_ARRAY=2, RTOL_ARRAY=3, ALL_ARRAY=4} CODEPACK_ITOLERANCE;
typedef enum codepack_itask_e {NORMAL = 1, ONESTEP_ONLY = 2, STOP_AT_FIRST_MESH = 3, NORMAL_TCRIT = 4, ONESTEP_TCRIT = 5} CODEPACK_ITASK;
typedef enum codepack_istate_in_e {FIRST_CALL = 1, NEXT_CALL = 2, NEXT_CALL_WITH_CHANGES = 3} CODEPACK_ISTATE_IN;
typedef enum codepack_istate_out_e { NOTHING_DONE           = 1,
                            SUCCESS_DONE           = 2,
                            MAX_STEPS_EXCEEDED     = -1,
                            TO_MUCH_ACCURACY       = -2,
                            ILLEGAL_INPUT          = -3 ,
                            ERROR_TEST_FAILURES    = -4,
                            CONVERGENCE_FAILURES   = -5,
                            ZERO_ERR_TOLERANCE     = -6,
                            TOO_SMALL_WORK_ARRAY   = -7 } CODEPACK_ISTATE_OUT;

typedef enum codepack_iopt_e {NO_OPTIONAL_INPUTS = 0, OPTIONAL_INPUTS = 1} CODEPACK_OPTIONAL_INPUT_FLAG;
typedef enum codepack_method_e {
    ADAMS_BASIC               = 10,
    ADAMS_USER_FULL_JAC       = 11,
    ADAMS_INTERNAL_FULL_JAC   = 12,
    ADAMS_INTERNAL_DIAG_JAC   = 13,
    ADAMS_USER_BAND_JAC       = 14,
    ADAMS_INTERNAL_BAND_JAC   = 15,
    BDF_BASIC                 = 20,
    BDF_USER_FULL_JAC         = 21,
    BDF_INTERNAL_FULL_JAC     = 22,
    BDF_INTERNAL_DIAG_JAC     = 23,
    BDF_USER_BAND_JAC         = 24,
    BDF_INTERNAL_BAND_JAC     = 25,
    BDF_USER_JAC_NO_IA_JA     = 121,
    BDF_INTERNAL_JAC_IA_JA    = 222
} CODEPACK_METHOD_FLAG;

typedef enum codepack_jac_type_e {USER_PROVIDED = 1, INTERNAL = 2, USER_PROVIDED_BANDED = 4, INTERNAL_BANDED = 5} CODEPACK_JAC_TYPE;

typedef enum codepack_ode_err_e {UNCHANGED = 2, SUCCESS = 0, ERROR = -1, UNKNOWN_ERROR = -20} CODEPACK_ODE_RETVAL;

CODEPACK_ODE_RETVAL lsode_basic(double* stack, double* q, codepack_ode_func f_func,int neq, double t0, double tf, double dt, CODEPACK_METHOD_FLAG mf);
CODEPACK_ODE_RETVAL lsodes_basic(double* stack, double* y, codepack_ode_func f_func,int neq, double t0, double tf, double dt, CODEPACK_METHOD_FLAG mf);
CODEPACK_ODE_RETVAL lsoda_basic(double* stack, double* q, codepack_ode_func f_func,int neq, double t0, double tf, double dt);
CODEPACK_ODE_RETVAL lsodar_basic(double* stack, double* q, codepack_ode_func f_func,int neq, double t0, double tf, double dt);

/*
 * full fortran interface here
 */


void dlsode(void (*f)(const int *neq, const double *t, const double *y, double *ydot),
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
            void (*jac)(const int *neq, const double *t, const double *y, const int *ml, const int *mu, double *pd, const int *nrowpd),
            const int *mf);

void dlsodes(void (*f)(const int *neq, const double *t, const double *y, double *ydot),
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
             void (*jac)(const int *neq, const double *t, const double *y, const int *j, const double *ian, double *jan, double *pdj),
             const int *mf);

void dlsoda(void (*f)(const int *neq, const double *t, const double *y, double *ydot),
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
            void (*jac)(const int *neq, const double *t, const double *y, const int *ml, const int *mu, double *pd, const int *nrowpd),
            const int *jt);


void dlsodar(void (*f)(const int *neq, const double *t, const double *y, double *ydot),
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
             void (*jac)(const int *neq, const double *t, const double *y, const int *ml, const int *mu, double *pd, const int *nrowpd), const int *jt, \
             void (*g)(const int *neq, const double *t, const double *y, const int *ng, double *gout), const int *ng,
             int *jroot);


#endif

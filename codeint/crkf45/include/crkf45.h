#ifndef H__CRKF45__H
#define H__CRKF45__H

typedef void (*rkf45_ode_func)(double* t, double* y, double* yp);

typedef enum rkf45_iflag_in_e  {ONE_STEP = -1, NORMAL=1} RKF45_IFLAG;

typedef enum rkf45_retval_e {ONE_STEP_SUCCES = -2,
                             SUCCESS=2,
                             RELATIVE_ERROR_TOO_SMALL=3,
                             TOO_MANY_EVAKIUATIONS=4,
                             SOLUTION_VANISHED=5,
                             REQUESTED_ACCURACY_TOO_SMALL=6,
                             USE_ONE_STEP_MODE=7,
                             INVALID_INPUT=8} RKF45_RETVAL;

RKF45_RETVAL rkf45_basic(double* stack, double* q,  rkf45_ode_func f_func,int neq, double t0, double tf, double dt);

void rkf45 ( void (*f)(double*, double* , double*),
             int* neqn,
             double* y,
             double* t,
             double* tout,
             double* relerr,
             double* abserr, 
             int* iflag,
             double* work,
             int* iwork );

#endif

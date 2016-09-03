#ifndef H__RKF45_INTERNAL__H
#define H__RKF45_INTERNAL__H
#include "../../codeintdeps/include/stack.h" 
#include "../include/crkf45.h"

typedef struct rkf45_params_s {

    rkf45_ode_func f_func;
    int neq;
    double rtol;
    double atol;
    union {
        RKF45_IFLAG  iflag;
        RKF45_RETVAL retval;
    } iflag;
    double* rwork;
    int*    iwork;
    int lrw;
    int liw;
} rkf45_params;

double* write_to_stack(double* stack, int neq, int* index, double t_new, double* q);
extern void r8_rkf45_ ( void (*f)(double*, double* , double*),
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

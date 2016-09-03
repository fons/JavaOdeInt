
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <czvode.h>


static void simple_pendulum(const int *neq, const double *t_, const ZVODE_COMPLEX *q, ZVODE_COMPLEX *qdot, ZVODE_COMPLEX* rpar, int* ipar)
{

    double alpha = 1;

    qdot[0].dr = q[1].dr;
    qdot[0].di = 0.0;
    qdot[1].dr = - alpha * sin(q[0].dr);
    qdot[1].di = 0.0;

    return;
}



static double* create_stack(double t0, double tf, double dt, int neq)
{
    double* s;
    int size = (tf - t0) / dt + 2;
    s  = (double*) calloc(size * (2*neq + 1) + 2, sizeof(double));
    *s = (double) size;
    *(s+1) = (double) neq;
    return s+2;
}

static void print_stack(FILE* fn, double* stack)
{
    double *f = (stack - 2);
    int size = (int) *f;
    int neq  = (int) *(f+1);
    int all  = size * (2 * neq + 1);
    for (int k = 0 ; k < all; k++) {
        fprintf(fn, "%.15lf,", stack[k]);
        
        if (((k + 1) % (2*neq + 1)) == 0) {
            fprintf(fn,"\n");
        }
    }
    fprintf(fn,"\n");
}

int main()
{
    double t0 = 0.0;
    double tf = 100.0;
    double dt = 0.1;
    int    neq = 2;
    ZVODE_COMPLEX q[2];

    q[0].dr = M_PI * 999/1000.0;
    q[0].di = 0.0;
    q[1].dr = 0.;
    q[1].di = 0.;

    double *stack = create_stack(t0, tf, dt, neq);
    CZVODE_ODE_RETVAL ret = zvode_basic(stack, q, &simple_pendulum, neq, t0, tf, dt, ADAMS); 
    if (ret != SUCCESS) {
        fprintf(stderr, "error encountered in vode_basic\n");
        return -1;
    }
    print_stack(stderr, stack);
    
    return ret;
}

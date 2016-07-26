
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <codepack.h>

static void simple_pendulum(const int *neq, const double *t_, const double *q, double *qdot)
{

    double alpha = 1;

    qdot[0] = q[1];
    qdot[1] = - alpha * sin(q[0]);

    return;
}



static double* create_stack(double t0, double tf, double dt, int neq)
{
    double* s;
    int size = (tf - t0) / dt + 2;
    s  = (double*) calloc(size * (neq + 1), sizeof(double));
    *s = (double) size;
    *(s+1) = (double) neq;
    return s+2;
}

static void print_stack(FILE* fn, double* stack)
{
    double *f = (stack - 2);
    int size = (int) *f;
    int neq  = (int) *(f+1);
    int all  = size * (neq + 1);
    //printf("%d %d  \n",size, neq);
    for (int k = 0 ; k < all - 1; k++) {
        fprintf(fn, "%.15lf,", stack[k]);
        if (((k + 1) % (neq + 1)) == 0) {
            fprintf(fn,"\n");
        }
    }
}

int main()
{
    double t0 = 0.0;
    double tf = 100.0;
    double dt = 0.1;
    int    neq = 2;
    double q[2];

    q[0] = M_PI * 999/1000.0;
    q[1] = 0.;

    double *stack = create_stack(t0, tf, dt, neq);
    CODEPACK_ODE_RETVAL ret = lsodar_basic(stack, q, &simple_pendulum, neq, t0, tf, dt); 
    if (ret != SUCCESS) {
        fprintf(stderr, "error encountered in lsodar_basic\n");
        return -1;
    }
    print_stack(stderr, stack);
    
    return ret;
}

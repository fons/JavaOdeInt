
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../codeintdeps/include/stack.h"
#include <codepack.h>

static void vanderpol(const int *neq, const double *t_, const double *q, double *qdot)
{

    double mu = 1000;

    qdot[0] = q[1];
    qdot[1] = mu * ( 1.0 - q[0]*q[0]) * q[1] - q[0];

    return;
}


int main()
{
    double t0 = 0.0;
    double tf = 2000.0;
    double dt = 0.01;
    int    neq = 2;
    double q[2];

    q[0] = 2.0;
    q[1] = 0.0;

    double *stack = create_stack(t0, tf, dt, neq);
    CODEPACK_ODE_RETVAL ret = lsodes_basic(stack, q, &vanderpol, neq, t0, tf, dt, ADAMS_BASIC); 
    if (ret != SUCCESS) {
        fprintf(stderr, "error encountered in lsodes_basic %d\n", ret);
        return -1;
    }
    print_stack(stdout, stack);
    
    return ret;
}

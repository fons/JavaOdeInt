
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../codeintdeps/include/stack.h"
#include <codepack.h>

static void simple_pendulum(const int *neq, const double *t_, const double *q, double *qdot)
{

    double alpha = 1;

    qdot[0] = q[1];
    qdot[1] = - alpha * sin(q[0]);

    return;
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
    CODEPACK_ODE_RETVAL ret = lsoda_basic(stack, q, &simple_pendulum, neq, t0, tf, dt); 
    if (ret != SUCCESS) {
        fprintf(stderr, "error encountered in lsodar_basic\n");
        return -1;
    }
    print_stack(stderr, stack);
    
    return ret;
}


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../codeintdeps/include/stack.h"
#include <cradau5.h>

static void simple_pendulum(int *neq, double *t_, double *q, double *qdot, double* rpar, int* ipar)
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
    q[1] = 0.00000000;

    double *stack = create_stack(t0, tf, dt, neq);
    RADAU5_RETVAL ret = radau5_basic(stack, q, &simple_pendulum, neq, t0, tf, dt); 
    if (ret != SUCCESS) {
        fprintf(stderr, "error encountered in radau5_basic\n");
        return -1;
    }
    print_stack(stderr, stack);
    fprintf(stderr, "\n");

    return ret;
}

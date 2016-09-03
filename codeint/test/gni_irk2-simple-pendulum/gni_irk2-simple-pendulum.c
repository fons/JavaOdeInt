
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../codeintdeps/include/stack.h"
#include <cgnicodes.h>

static void simple_pendulum(int *n, double* x, double* q, double* f, double* rpar, int* ipar)
{
    double alpha = 1;
    f[0] = - alpha * sin(q[0]);
    return;
}


int main()
{
    
    double t0 = 0.0;
    double tf = 100.0;
    double dt = 0.1;
    int    neq = 1;
    double q[1];
    double p[1];
    
    q[0] = M_PI * 999/1000.0;
    p[0] = 0.000000;

    double *stack = create_stack(t0, tf, dt, neq+1);

    gni_irk2_basic(stack, p,q, &simple_pendulum, neq, t0, tf, dt, IRK2_METH_6); 

    print_stack(stderr, stack);
    fprintf(stderr, "\n");

    return 0;
    
}

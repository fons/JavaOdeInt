
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../codeintdeps/include/stack.h"
#include <crkf45.h>

static void simple_pendulum(double *t_, double *q, double *qdot)
{

    double alpha = 1;
/*
    fprintf(stderr, "%lf %lf %lf\n", *t_, q[0], qdot[0]);    
*/
    qdot[0] = q[1];
    qdot[1] = - alpha * sin(q[0]);

    return;
}

#ifdef KLM
static void simple_1d(const double *x, const double *y, double *f)
{

    /*
      dy/dx = 2 * x => y = x^2
    */
    f[0] = 2 * (*x);
/*
  fprintf(stderr, "%d %lf %lf %lf\n", *n, *x, y[0], f[0]);
*/
    return;
}

static void simple_2d(const double *x, const double *y, double *f)
{
    double omega = 2.0 * M_PI;
    /*
      dy_0/dx = - y_1
      dy_1/dx =   y_0
    */
    f[0] = - omega * y[1];
    f[1] =   omega * y[0];
/*
    fprintf(stderr, "n : %d x : %lf y[0] %lf y[1] %lf f[0] : %lf f[1] : %lf\n", *n, *x, y[0],y[1], f[0], f[1]);
*/
    return;
}
#endif



int main()
{
/*
      double x0 = 0.0;
      double xf = 1.0;
      double dx = 0.1;
      int    neq = 1;
      double y[1];

      y[0] = 0;


      double *stack = create_stack(x0, xf, dx, neq);
      RKF45_RETVAL ret = rkf45_basic(stack, y, &simple_1d, neq, x0, xf, dx); 
      if (ret != SUCCESS) {
      fprintf(stderr, "error encountered in rkf45_basic\n");
      return -1;
      }
      print_stack(stderr, stack);

    
    double x0 = 0.0;
    double xf = 7.0;
    double dx = 0.1;
    int    neq = 2;
    double y[2];

    y[0] = 1.0;
    y[1] = -0.25;

    double *stack = create_stack(x0, xf, dx, neq);
    RKF45_RETVAL ret = rkf45_basic(stack, y, &simple_2d, neq, x0, xf, dx); 
    if (ret != SUCCESS) {
    fprintf(stderr, "error encountered in rkf45_basic\n");
        return -1;
    }
    print_stack(stderr, stack);
    fprintf(stderr, "\n");
*/
    double t0 = 0.0;
    double tf = 100.0;
    double dt = 0.1;
    int    neq = 2;
    double q[2];

    q[0] = M_PI * 999/1000.0;
    q[1] = 0.00000000;

    double *stack = create_stack(t0, tf, dt, neq);
    RKF45_RETVAL ret = rkf45_basic(stack, q, &simple_pendulum, neq, t0, tf, dt); 
    print_stack(stderr, stack);
    if (ret != SUCCESS) {
        fprintf(stderr, "error encountered in rkf45_basic %d : \n", ret);
        return -1;
    }
    fprintf(stderr, "\n");

    return ret;
}

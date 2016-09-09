#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cquadpack.h"

static double f1(double* x)
{
    return *x + 4.0 * (*x)*(*x);
}

static double f2(double* x)
{
    return exp(-0.4 * (*x) * (*x));
}

static double f3(double *x)
{
    if ((*x) > 0) {
        return log(*x);
    }
    return 0.0;
}

static double f4(double *x)
{
    return 1.0/(5.0 *(*x)*(*x)*(*x) + 6.0);
}

static double f5(double *x)
{
    if ((*x) > 0) {
        return 1.0/sqrt(*x);
    }
    return 0.0;
}

static double f6(double *x)
{
    if ((*x) > 0) {
        double v = 1.0 + log(*x)*log(*x);
        return 1.0 / (v*v);
    }
    return 0.0;
}

static double f7(double* x)
{
    double xsq = (*x)*(*x);
    return xsq * (*x) * log (fabs( (xsq - 1.0) * (xsq - 2.0)));
}

int main(int argc, char** argv)
{

    int ret;
    double result;
    ret = qagse_basic(f1,10.0,6.0,&result);
    fprintf(stderr, "\n 1 qagse_basic ==>%d %lf\n", ret, result);

    ret = qagie_basic(f2,1.0,MINUS_INFINITY ,&result);
    fprintf(stderr, "\n 2 qagie_basic ==>%d %lf\n", ret, result);

    ret = qagie_basic(f2,1.0,PLUS_INFINITY ,&result);
    fprintf(stderr, "\n 3 qagie_basic ==>%d %lf\n", ret, result);

    ret = qagie_basic(f2,10.0,MINUS_TO_PLUS_INFINITY ,&result);
    fprintf(stderr, "\n 4 qagie_basic ==>%d %lf\n", ret, result);
    
    double omega = 10.0 * M_PI;
    ret = qawoe_basic(f3,0.0001, 1.0 ,omega, SIN, &result);
    fprintf(stderr, "\n 5 qawoe_basic ==>%d %lf\n", ret, result);

    ret = qawce_basic(f4,-1.0, 5.0 ,0.0, &result);
    fprintf(stderr, "\n 6 qawce_basic ==>%d %lf\n", ret, result);


    omega = 0.50 * M_PI;
    ret = qawfe_basic(f5,0.00, omega ,COS, &result);
    fprintf(stderr, "\n 7 qawfe_basic ==>%d %lf\n", ret, result);

    fprintf(stderr, "%s\n", log_weight_to_string(LOGW_1));
    fprintf(stderr, "%s\n", log_weight_to_string(LOGW_2));
    fprintf(stderr, "%s\n", log_weight_to_string(LOGW_3));
    fprintf(stderr, "%s\n", log_weight_to_string(LOGW_4));

    ret = qawse_basic(f6, 0.0, 1.0 , 0.0,0.0, LOGW_2,  &result);
    fprintf(stderr, "\n 8 qawse_basic ==>%d %lf\n", ret, result);


    int npts2 = 4;

    double* points = (double*)calloc(npts2, sizeof(double));
    points[0] = 1.0;
    points[1] = sqrt(2.0);

    ret = qagpe_basic(f7, 0.0, 3.0 ,  npts2, points,  &result);

    fprintf(stderr, "\n 9 qagpe_basic ==>%d %lf\n", ret, result);

    free (points);
    return 0;
}

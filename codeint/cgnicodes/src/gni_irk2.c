#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "gnicodes-internal.h"


static gnicodes_params* create_basic_gnicodes_params(int neq, gnicodes_ode_func f_func, int nstep, GNI_IRK2_METH meth)
{
    
    gnicodes_params* dsp = calloc(1, sizeof(gnicodes_params));
    if(dsp == NULL) {
        return NULL;
    }
    dsp->neq            = neq;
    dsp->f_func         = f_func;
    dsp->solfix         = NULL;
    dsp->iout           = NEVER_CALLED;
    dsp->nstep          = nstep;
    dsp->meth.irk2_meth = meth;
    dsp->lr             = 10;
    dsp->rpar           = (double *) calloc(dsp->lr, sizeof(double));
    dsp->li             = 10;
    dsp->ipar           = (int *) calloc(dsp->li, sizeof(int));

    return dsp;
}

static void gni_irk2_w(gnicodes_params* dsp, double xend, double* x, double* p, double* q)
{
    gni_irk2_(&dsp->neq,
              dsp->f_func,
              &dsp->nstep,  
              x,
              p,
              q,
              &xend,
              (int *) &dsp->meth.irk2_meth,
              dsp->solfix,
              (int *)&dsp->iout,
              dsp->rpar,
              dsp->ipar);
}

static void simple_pendulum(int *n, double* x, double* q, double* f, double* rpar, int* ipar)
{
    double alpha = 1;
    f[0] = - alpha * sin(q[0]);
    return;
}
/*
static gnicodes_ode_func F;

static void wrapper(int *n, double* x, double* q, double* f, double* rpar, int* ipar)
{
    struct timeval t1, t2;
    double elapsedTime;

    gettimeofday(&t1, NULL);       

    F(n,x,q,f,rpar,ipar);
    gettimeofday(&t2, NULL);
    elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;
    fprintf(stderr, "\t\t\t\%lf \n", elapsedTime);

    return;
}

*/
void gni_irk2_basic(double* stack, double* p, double* q, gnicodes_ode_func f_func, int n, double xstart, double xfinal, double deltax, GNI_IRK2_METH meth)
{
/*
    struct timeval t1, t2;
    double elapsedTime;
*/  
    double x = xstart;
    int index = 0;
    int steps = abs((int) ((xfinal - xstart)/deltax));
    //fprintf(stderr, "start integration loop\n");
    gnicodes_params* dlsop = create_basic_gnicodes_params(n, f_func, steps, meth);
/*
    F = simple_pendulum;

    gnicodes_params* dlsop = create_basic_gnicodes_params(n, &wrapper, steps, meth);
*/  
    stack = write_pq_to_stack(stack, n, &index, x, p, q);
    while ( x < xfinal) {
        //gettimeofday(&t1, NULL);       
        gni_irk2_w(dlsop, x + deltax, &x, p, q);
        //gettimeofday(&t2, NULL);
        //elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
        //elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;
        //fprintf(stderr, "MAIN %lf \n", elapsedTime);
        stack = write_pq_to_stack(stack, n, &index, x, q, p);
    }
    //fprintf(stderr, "DONE integration loop\n");

}

void gni_irk2(
    int* n,
    void(*f)(int* n, double* x, double* q, double* f, double* rpar, int* ipar),
    int* nstep,
    double* x,
    double* p,
    double* q,
    double* xend,
    int* meth,
    void (*s) (int* nr, double* xold, double* x, double* p, double* q, int* n, int* irtrn, double* rpar, int* ipar),
    int* iout,
    double* rpar,
    int* ipar
    )
{
    gni_irk2_(n,f,nstep, x,p,q,xend,meth,s,iout,rpar,ipar);
}

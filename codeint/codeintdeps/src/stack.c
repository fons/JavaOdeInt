#include <stdlib.h>
#include <stdio.h>

double* write_pq_to_stack(double* stack, int neq, int* index, double t_new, double* pp, double* qq)
{
    double *p = (stack + *index);
    (*index)++;
    if (p != NULL) {
        *p = t_new;
        for (int k = 1; k < (neq+1);k++) {
            *(p + k) = pp[k-1];
            (*index)++;
        }
        for (int k = 1; k < (neq+1);k++) {
            *(p + neq + k) = qq[k-1];
            (*index)++;
        }
    }
    else {
        fprintf(stderr, "stack element is null \n");
    }
    return stack;
}
double* write_to_stack(double* stack, int neq, int* index, double t_new, double* q)
{
    double *p = (stack + *index);
    (*index)++;
    if (p != NULL) {
        *p = t_new;
        for (int k = 1; k < (neq+1);k++) {
            *(p + k) = q[k-1];
            (*index)++;
        }
    }
    else {
        fprintf(stderr, "stack element is null \n");
    }
    return stack;
}


double* create_stack(double t0, double tf, double dt, int neq)
{
    double* s;
    int size = (tf - t0) / dt + 2;
    s  = (double*) calloc(size * (neq + 1)+2, sizeof(double));
    *s = (double) size;
    *(s+1) = (double) neq;
    return s+2;
}

void print_stack(FILE* fn, double* stack)
{
    double *f = (stack - 2);
    int size = (int) *f;
    int neq  = (int) *(f+1);
    int all  = size * (neq + 1);
    //printf("%d %d  \n",size, neq);
    for (int k = 0 ; k < all; k++) {
        fprintf(fn, "%.15lf,", stack[k]);
        if (((k + 1) % (neq + 1)) == 0) {
            fprintf(fn,"\n");
        }
    }
}

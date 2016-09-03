#include <stdlib.h>
#include <stdio.h>
#include "czvode-internal.h"

double* write_to_stack(double* stack, int neq, int* index, double t_new, ZVODE_COMPLEX* q)
{
    double *p = (stack + *index);
    if (p != NULL) {
        *p = t_new;
        (*index)++;
        for (int k = 1; k < (2*neq+1);k=k+2) {
            int j = (k - 1)/2;
            *(p + k)     = q[j].dr;
            (*index)++;
            *(p + k + 1) = q[j].di;
            (*index)++;
        }
    }
    else {
        fprintf(stderr, "stack element is null \n");
    }
    return stack;
}

#include <stdlib.h>
#include <stdio.h>

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

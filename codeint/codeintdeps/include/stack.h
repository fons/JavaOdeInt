#ifndef H__STACK__H
#define H__STACK__H
#include <stdio.h>
double* write_pq_to_stack(double* stack, int neq, int* index, double t_new, double* pp, double* qq);
double* write_to_stack(double* stack, int neq, int* index, double t_new, double* q);
void print_stack(FILE* fn, double* stack);
double* create_stack(double t0, double tf, double dt, int neq);
#endif

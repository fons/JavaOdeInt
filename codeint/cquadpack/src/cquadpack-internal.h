#ifndef H__CQUADPACK_INTERNAL__H
#define H__CQUADPACK_INTERNAL__H

extern void dqagse_(double (*f)(double *y),
                    double* a,
                    double* b,
                    double* epsabs,
                    double* epsrel,
                    int* limit,
                    double* result,
                    double* abserr,
                    int* neval,
                    int* ier,
                    double* alist,
                    double* blist,
                    double* rlist,
                    double* elist,
                    int* iord,
                    int* last);

extern void dqagie_(double (*f)(double *y),
                    double* a,
                    int* inf,
                    double* epsabs,
                    double* epsrel,
                    int* limit,
                    double* result,
                    double* abserr,
                    int* neval,
                    int* ier,
                    double* alist,
                    double* blist,
                    double* rlist,
                    double* elist,
                    int* iord,
                    int* last);

extern void dqawoe_(double (*f)(double *y),
                    double* a,
                    double* b,
                    double* omega,
                    int* integr,
                    double* epsabs,
                    double* epsrel,
                    int* limit,
                    int* icall,
                    int* maxpl,
                    double* result,
                    double* abserr,
                    int* neval,
                    int* ier,
                    int* last,
                    double* alist,
                    double* blist,
                    double* rlist,
                    double* elist,
                    int* iord,
                    int* nnlog,
                    int* momcom,
                    double* chebmo
    );

extern void dqawce_(double (*f)(double *y),
                    double* a,
                    double* b,
                    double* c,
                    double* epsabs,
                    double* epsrel,
                    int* limit,
                    double* result,
                    double* abserr,
                    int* neval,
                    int* ier,
                    double* alist,
                    double* blist,
                    double* rlist,
                    double* elist,
                    int* iord,
                    int* last);

extern void dqawfe_(double (*f)(double *y),
                    double* a,
                    double* omega,
                    int* integr,
                    double* epsabs,
                    int* limlst,
                    int* limit,
                    int* maxpl,
                    double* result,
                    double* abserr,
                    int* neval,
                    int* ier,
                    double* rslst,
                    double* erlst,
                    int*    ierlst,
                    int* lst,
                    double* alist,
                    double* blist,
                    double* rlist,
                    double* elist,
                    int* iord,
                    int* nnlog,
                    double* chebmo
    );

extern void dqawse_(double (*f)(double *y),
                    double* a,
                    double* b,
                    double* alfa,
                    double* beta,
                    int* integr,
                    double* epsabs,
                    double* epsrel,
                    int* limit,
                    double* result,
                    double* abserr,
                    int* neval,
                    int* ier,
                    double* alist,
                    double* blist,
                    double* rlist,
                    double* elist,
                    int* iord,
                    int* last);

extern void dqagpe_(double (*f)(double *y),
                    double* a,
                    double* b,
                    int* npts2,
                    double* points,
                    double* epsabs,
                    double* epsrel,
                    int* limit,
                    double* result,
                    double* abserr,
                    int* neval,
                    int* ier,
                    double* alist,
                    double* blist,
                    double* rlist,
                    double* elist,
                    double* pts,
                    int* level,
                    int*  ndin,
                    int* iord,
                    int* last);

typedef struct quadpack_input_s {
    double epsabs;
    double epsrel;
    int    limlst;
    int    limit;
    int    icall;
    int    maxpl;
} QUADPACK_INPUT;

typedef struct quadpack_output_s {
    double  abserr;
    int     neval;
    int     ier;
    double* rslst;
    double* erlst;
    int*    ierlst;
    int     lst;
    double* alist;
    double* blist;
    double* rlist;
    double* elist;
    int*    iord;
    int     last;
    int*    nnlog;
    int     momcom;
    double* chebmo;
    double* pts;
    int*    level;
    int*    ndin;
} QUADPACK_OUTPUT;

typedef double QUADPACK_RESULT;


void default_quadpack_input(QUADPACK_INPUT* qpi);
void default_quadpack_output(QUADPACK_INPUT qpi, QUADPACK_OUTPUT* qpo);
void points_quadpack_output(int npts2, QUADPACK_OUTPUT* qpo);
void free_quadpack_output(QUADPACK_OUTPUT* qpo);
#endif

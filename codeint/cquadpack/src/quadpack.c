#include <stdlib.h>
#include <stdio.h>
#include "../include/cquadpack.h"
#include "cquadpack-internal.h"


void default_quadpack_input(QUADPACK_INPUT* qpi)
{
    qpi->epsabs = 1.49e-8;
    qpi->epsrel = 1.49e-8;
    qpi->limit  = 50;
    qpi->icall  = 1;
    qpi->maxpl  = 25;
    qpi->limlst = 20;
}

void default_quadpack_output(QUADPACK_INPUT qpi, QUADPACK_OUTPUT* qpo)
{
    
    qpo->alist   = (double*)  calloc(qpi.limit, sizeof(double));
    qpo->blist   = (double*)  calloc(qpi.limit, sizeof(double));
    qpo->rlist   = (double*)  calloc(qpi.limit, sizeof(double));
    qpo->elist   = (double*)  calloc(qpi.limit, sizeof(double));
    qpo->iord    = (int*)     calloc(qpi.limit, sizeof(int));
    qpo->nnlog   = (int*)     calloc(qpi.limit, sizeof(int));
    qpo->chebmo  = (double* ) calloc(qpi.maxpl * 25, sizeof(double));
    qpo->momcom  = 0;
    qpo->rslst   = (double*)  calloc(qpi.limlst, sizeof(double));
    qpo->erlst   = (double*)  calloc(qpi.limlst, sizeof(double));
    qpo->ierlst  = (int*)     calloc(qpi.limlst, sizeof(int));
    qpo->level   = (int*)     calloc(qpi.limit, sizeof(int));
    qpo->pts     = NULL;
    qpo->ndin    = NULL;
}
void points_quadpack_output(int npts2, QUADPACK_OUTPUT* qpo)
{
    qpo->pts     = (double*)  calloc(npts2, sizeof(double));
    qpo->ndin    = (int*)     calloc(npts2, sizeof(int));
    
}
void free_quadpack_output(QUADPACK_OUTPUT* qpo)
{
    free(qpo->alist);
    free(qpo->blist);
    free(qpo->rlist);
    free(qpo->elist);
    free(qpo->iord);
    free(qpo->nnlog);
    free(qpo->chebmo);
    free(qpo->rslst);
    free(qpo->erlst);
    free(qpo->ierlst);
    free(qpo->level);
    free(qpo->pts);
    free(qpo->ndin);
    
}

char* log_weight_to_string(QUADPACK_LOG_WEIGHT_FUNCTION w)
{
    return LOG_WEIGHT[w-1];
}

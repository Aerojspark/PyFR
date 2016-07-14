# -*- coding: utf-8 -*-
#include <math.h>
#include <omp.h>

typedef double fpdtype_t;

typedef void (*xsmm_gemm_t)(const fpdtype_t *, const fpdtype_t *, fpdtype_t *);

void par_xsmm_gemm(xsmm_gemm_t gemm, const int N,
                   const fpdtype_t * b, const fpdtype_t * a, fpdtype_t * c)
{
    #pragma omp parallel for
    for (int _i=0; _i < N; _i+=${nblock})
        gemm(b + _i, a, c + _i);
}
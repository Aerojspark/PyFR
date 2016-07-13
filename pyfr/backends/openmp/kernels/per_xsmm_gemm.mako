# -*- coding: utf-8 -*-
<%inherit file='base'/>

typedef void (*xsmm_gemm_t)(const fpdtype_t *, const fpdtype_t *, fpdtype_t *);

void per_xsmm_gemm(xsmm_gemm_t gemm,
                   const fpdtype_t * B, const fpdtype_t * A, fpdtype_t * C)
{
    #pragma omp parallel
    {
        int begin, end;
        loop_sched_1d(N, PYFR_ALIGN_BYTES / sizeof(fpdtype_t), &begin, &end);

        gemm(B + begin, A, C + begin);
    }
}
# -*- coding: utf-8 -*-
<%inherit file='base'/>

// XSMM GEMM prototype
typedef void (*xsmm_gemm_t)(const fpdtype_t *, const fpdtype_t *, fpdtype_t *);

void par_xsmm_gemm(xsmm_gemm_t gemm, const int n,
                   const fpdtype_t * b, const fpdtype_t * a, fpdtype_t * c)
{
    #pragma omp parallel
    {
        // TODO: Block is the same as align bytes (which value is optimal?)
        int nblock = PYFR_ALIGN_BYTES / sizeof(fpdtype_t);

        int begin, end;
        loop_sched_1d(n, nblock, &begin, &end);
        for (int _i=begin; _i < end; _i+=nblock)
            gemm(b + _i, a, c + _i);
    }
}
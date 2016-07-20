# -*- coding: utf-8 -*-

from ctypes import (CDLL, CFUNCTYPE, RTLD_GLOBAL, POINTER, cast,
                    c_int, c_void_p, c_float, c_double)

import numpy as np

from pyfr.backends.base import ComputeKernel, NotSuitableError
from pyfr.backends.openmp.provider import OpenMPKernelProvider
from pyfr.ctypesutil import load_library


class XSMMWrappers(object):
    def __init__(self, cblas):
        # TODO: libxsmm requires dgemm rourine
        cblas = CDLL(cblas, mode=RTLD_GLOBAL)

        # Load XSMM library
        self.lib = lib = load_library('xsmm')

        # GEMM function prototype
        gemm = CFUNCTYPE(c_void_p, c_void_p, c_void_p, c_void_p)

        # Dispatch
        smmdispatch = lib.libxsmm_smmdispatch
        smmdispatch.argtypes = [c_int, c_int, c_int,
                                POINTER(c_int), POINTER(c_int), POINTER(c_int),
                                POINTER(c_float), POINTER(c_float),
                                POINTER(c_int), POINTER(c_int)]
        smmdispatch.restype = gemm

        dmmdispatch = lib.libxsmm_dmmdispatch
        dmmdispatch.argtypes = [c_int, c_int, c_int,
                                POINTER(c_int), POINTER(c_int), POINTER(c_int),
                                POINTER(c_double), POINTER(c_double),
                                POINTER(c_int), POINTER(c_int)]
        dmmdispatch.restype = gemm

    def dmmdispatch(self, m, n, k, lda=None, ldb=None, ldc=None,
                    alpha=None, beta=None, flags=None, prefetch=None):
        return self.lib.libxsmm_dmmdispatch(m, n, k, lda, ldb, ldc,
                                            alpha, beta, flags, prefetch)

    def smmdispatch(self, m, n, k, lda=None, ldb=None, ldc=None,
                    alpha=None, beta=None, flags=None, prefetch=None):
        return self.lib.libxsmm_smmdispatch(m, n, k, lda, ldb, ldc,
                                            alpha, beta, flags, prefetch)


class OpenMPXSMMKernels(OpenMPKernelProvider):
    def __init__(self, backend):
        super().__init__(backend)

        cblas = backend.cfg.getpath('backend-openmp', 'cblas', abs=False)

        # Load Wrapper
        try:
            self._wrappers = XSMMWrappers(cblas)
        except OSError:
            pass

    def mul(self, a, b, out, alpha=1.0, beta=0.0):
        try:
            w = self._wrappers
        except:
            raise NotSuitableError('Failed to load xsmm library')

        # Ensure the matrices are compatible
        if a.nrow != out.nrow or a.ncol != b.nrow or b.ncol != out.ncol:
            raise ValueError('Incompatible matrices for out = a*b')
        
        # Check size (ldb < 65536)
        if b.leaddim > 65536:
            raise NotSuitableError('Leaddim is too big')

        # TODO: Block is the same as align bytes (which value is optimal?)
        nblock = self.backend.alignb // a.itemsize

        m, n, k = a.nrow, b.ncol, a.ncol
        lda, ldb, ldc = c_int(a.leaddim), c_int(b.leaddim), c_int(out.leaddim)

        # α and β factors for C = α*(A*op(B)) + β*C
        # Column major C^T = B^T*A^T
        if a.dtype == np.float64:
            alpha_ct, beta_ct = c_double(alpha), c_double(beta)
            gemm = w.dmmdispatch(nblock, m, k, ldb, lda, ldc, alpha_ct, beta_ct)
        else:
            alpha_ct, beta_ct = c_float(alpha), c_float(beta)
            gemm = w.smmdispatch(nblock, m, k, ldb, lda, ldc, alpha_ct, beta_ct)

        argt = [np.intp, np.int32, np.intp, np.intp, np.intp]

        # Render the kernel template
        tplargs = dict(nblock=nblock)
        src = self.backend.lookup.get_template('par_xsmm_gemm').render(**tplargs)

        # Build
        par_gemm = self._build_kernel('par_xsmm_gemm', src, argt)

        gemm_ptr = cast(gemm, c_void_p).value

        class MulKernel(ComputeKernel):
            def run(self, queue):
                par_gemm(gemm_ptr, n, b, a, out)

        return MulKernel()

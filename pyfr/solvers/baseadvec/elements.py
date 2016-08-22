# -*- coding: utf-8 -*-

from pyfr.backends.base.kernels import ComputeMetaKernel
from pyfr.solvers.base import BaseElements


class BaseAdvectionElements(BaseElements):
    @property
    def _scratch_bufs(self):
        if 'flux' in self.antialias:
            bufs = {'scal_fpts', 'matr_upts', 'scal_qpts', 'vect_qpts'}
        elif 'div-flux' in self.antialias:
            bufs = {'scal_fpts', 'vect_upts', 'matr_upts', 'scal_qpts'}
        else:
            bufs = {'scal_fpts', 'vect_upts', 'matr_upts'}

        if self._soln_in_src_exprs:
            if 'div-flux' in self.antialias:
                bufs |= {'scal_qpts_cpy'}
            else:
                bufs |= {'scal_upts_cpy'}

        return bufs

    def set_backend(self, backend, nscal_upts):
        super().set_backend(backend, nscal_upts)

        # Register pointwise kernels with the backend
        backend.pointwise.register(
            'pyfr.solvers.baseadvec.kernels.negdivconf'
        )

        # What anti-aliasing options we're running with
        fluxaa = 'flux' in self.antialias
        divfluxaa = 'div-flux' in self.antialias

        # What the source term expressions (if any) are a function of
        plocsrc = self._ploc_in_src_exprs
        solnsrc = self._soln_in_src_exprs

        # Source term kernel arguments
        srctplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'srcex': self._src_exprs
        }

        # Interpolation from elemental points
        if fluxaa or (divfluxaa and solnsrc):
            self.kernels['disu'] = lambda: backend.kernel(
                'mul', self.opmat('M8'), self.scal_upts_inb,
                out=self._scal_fqpts
            )
        else:
            self.kernels['disu'] = lambda: backend.kernel(
                'mul', self.opmat('M0'), self.scal_upts_inb,
                out=self._scal_fpts
            )

        # Interpolations and projections to/from quadrature points
        if divfluxaa:
            self.kernels['divf_qpts'] = lambda: backend.kernel(
                'mul', self.opmat('M7'), self.scal_upts_outb,
                out=self._scal_qpts
            )
            self.kernels['divf_upts'] = lambda: backend.kernel(
                'mul', self.opmat('M9'), self._scal_qpts,
                out=self.scal_upts_outb
            )

        # First flux correction kernel
        if fluxaa:
            nupts, nfpts, nqpts = self.nupts, self.nfpts, self.nqpts
            ndims = self.ndims

            def tdivdisf():
                vectu, matru = self._vect_qpts, self._matr_upts
                ndupts = ndims * nupts
                muls = [backend.kernel(
                    'mul', self.opmat('M4*M9'),
                    vectu.rslice(i * nqpts, (i + 1) * nqpts),
                    matru.rslice(i * ndupts, (i + 1) * ndupts)
                ) for i in range(ndims)]

                return ComputeMetaKernel(muls)

            def disf_fpts():
                vectu, vectf = self._vect_qpts, self._vect0_fpts
                muls = [backend.kernel(
                    'mul', self.opmat('M0*M9'),
                    vectu.rslice(i * nqpts, (i + 1) * nqpts),
                    vectf.rslice(i * nfpts, (i + 1) * nfpts)
                ) for i in range(ndims)]

                return ComputeMetaKernel(muls)
        else:
            nupts, nfpts, ndims = self.nupts, self.nfpts, self.ndims

            def tdivdisf():
                vectu, matru = self._vect_upts, self._matr_upts
                ndupts = ndims * nupts
                muls = [backend.kernel(
                    'mul', self.opmat('M4'),
                    vectu.rslice(i * nupts, (i + 1) * nupts),
                    matru.rslice(i * ndupts, (i + 1) * ndupts)
                ) for i in range(ndims)]

                return ComputeMetaKernel(muls)

            def disf_fpts():
                vectu, vectf = self._vect_upts, self._vect0_fpts
                muls = [backend.kernel(
                    'mul', self.opmat('M0'),
                    vectu.rslice(i * nupts, (i + 1) * nupts),
                    vectf.rslice(i * nfpts, (i + 1) * nfpts)
                ) for i in range(ndims)]

                return ComputeMetaKernel(muls)

        self.kernels['tdivdisf'] = tdivdisf
        self.kernels['disf_fpts'] = disf_fpts

        # Second flux correction kernel
        self.kernels['divconf'] = lambda: backend.kernel(
            'mul', self.opmat('M3'), self._scal_fpts, out=self.scal_upts_outb,
            beta=1.0
        )

        tplargs = {'nvars': self.nvars, 'ndims': self.ndims}
        backend.pointwise.register('pyfr.solvers.baseadvec.kernels.divdisf')
        self.kernels['divdisf'] = lambda: backend.kernel(
            'divdisf', dims=[self.nupts, self.neles], tplargs=tplargs,
            f=self._matr_upts, tsmats=self.smat_tr_at('upts'),
            u=self.scal_upts_outb)

        # Transformed to physical divergence kernel + source term
        if divfluxaa:
            plocqpts = self.ploc_at('qpts') if plocsrc else None
            solnqpts = self._scal_qpts_cpy if solnsrc else None

            if solnsrc:
                self.kernels['copy_soln'] = lambda: backend.kernel(
                    'copy', self._scal_qpts_cpy, self._scal_qpts
                )

            self.kernels['negdivconf'] = lambda: backend.kernel(
                'negdivconf', tplargs=srctplargs,
                dims=[self.nqpts, self.neles], tdivtconf=self._scal_qpts,
                rcpdjac=self.rcpdjac_at('qpts'), ploc=plocqpts, u=solnqpts
            )
        else:
            plocupts = self.ploc_at('upts') if plocsrc else None
            solnupts = self._scal_upts_cpy if solnsrc else None

            if solnsrc:
                self.kernels['copy_soln'] = lambda: backend.kernel(
                    'copy', self._scal_upts_cpy, self.scal_upts_inb
                )

            self.kernels['negdivconf'] = lambda: backend.kernel(
                'negdivconf', tplargs=srctplargs,
                dims=[self.nupts, self.neles], tdivtconf=self.scal_upts_outb,
                rcpdjac=self.rcpdjac_at('upts'), ploc=plocupts, u=solnupts
            )

        # In-place solution filter
        if self.cfg.getint('soln-filter', 'nsteps', '0'):
            def filter_soln():
                mul = backend.kernel(
                    'mul', self.opmat('M11'), self.scal_upts_inb,
                    out=self._scal_upts_temp
                )
                copy = backend.kernel(
                    'copy', self.scal_upts_inb, self._scal_upts_temp
                )

                return ComputeMetaKernel([mul, copy])

            self.kernels['filter_soln'] = filter_soln

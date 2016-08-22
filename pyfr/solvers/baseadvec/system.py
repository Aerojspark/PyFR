# -*- coding: utf-8 -*-

from pyfr.solvers.base import BaseSystem


class BaseAdvectionSystem(BaseSystem):
    _nqueues = 2

    def rhs(self, t, uinbank, foutbank):
        runall = self.backend.runall
        q1, q2 = self._queues
        kernels = self._kernels

        self.eles_scal_upts_inb.active = uinbank
        self.eles_scal_upts_outb.active = foutbank

        q1 << kernels['eles', 'disu']()
        q1 << kernels['mpiint', 'scal_fpts_pack']()
        runall([q1])

        if ('eles', 'copy_soln') in kernels:
            q1 << kernels['eles', 'copy_soln']()
        q1 << kernels['eles', 'disf_upts']()
        q1 << kernels['eles', 'disf_fpts']()
        q1 << kernels['iint', 'comfmdisfn']()
        q1 << kernels['bcint', 'comfmdisfn'](t=t)

        q2 << kernels['mpiint', 'scal_fpts_send']()
        q2 << kernels['mpiint', 'scal_fpts_recv']()
        q2 << kernels['mpiint', 'scal_fpts_unpack']()

        runall([q1, q2])

        q1 << kernels['mpiint', 'comfmdisfn']()
        q1 << kernels['eles', 'tdivdisf']()
        q1 << kernels['eles', 'divdisf']()
        q1 << kernels['eles', 'divconf']()
        if ('eles', 'divf_qpts') in kernels:
            q1 << kernels['eles', 'divf_qpts']()
            q1 << kernels['eles', 'negdivconf'](t=t)
            q1 << kernels['eles', 'divf_upts']()
        else:
            q1 << kernels['eles', 'negdivconf'](t=t)
        runall([q1])

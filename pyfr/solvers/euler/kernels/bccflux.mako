# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.euler.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.euler.kernels.bcs.${bctype}'/>

<%pyfr:kernel name='bccflux' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              fl='in view fpdtype_t[${str(ndims)}][${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              magnl='in fpdtype_t'
              ploc='in fpdtype_t[${str(ndims)}]'
              t='scalar fpdtype_t'>
    // Compute the RHS
    fpdtype_t ur[${nvars}];
    ${pyfr.expand('bc_rsolve_state', 'ul', 'nl', 'ur', 'ploc', 't')};

    // Perform the Riemann solve
    fpdtype_t fn[${nvars}];
    ${pyfr.expand('rsolve', 'ul', 'ur', 'nl', 'fn')};

    // Scale and write out the common normal fluxes
    fpdtype_t fln;
% for i in range(nvars):
    fln = ${'+'.join('nl[{0}]*fl[{0}][{1}]'.format(j, i)
                      for j in range(ndims))};
    ul[${i}] =  magnl*(fn[${i}] - fln);
% endfor
</%pyfr:kernel>

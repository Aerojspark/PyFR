# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.euler.kernels.rsolvers.${rsolver}'/>

<%pyfr:kernel name='intcflux' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              ur='inout view fpdtype_t[${str(nvars)}]'
              fl='in view fpdtype_t[${str(ndims)}][${str(nvars)}]'
              fr='in view fpdtype_t[${str(ndims)}][${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              magnl='in fpdtype_t'>
    // Perform the Riemann solve
    fpdtype_t fn[${nvars}];
    ${pyfr.expand('rsolve', 'ul', 'ur', 'nl', 'fn')};

    // Scale and write out the common normal fluxes
    fpdtype_t fln, frn;
% for i in range(nvars):
    fln = ${'+'.join('nl[{0}]*fl[{0}][{1}]'.format(j, i)
                      for j in range(ndims))};
    frn = ${'+'.join('nl[{0}]*fr[{0}][{1}]'.format(j, i)
                      for j in range(ndims))};
    ul[${i}] =  magnl*(fn[${i}] - fln);
    ur[${i}] = -magnl*(fn[${i}] - frn);
% endfor
</%pyfr:kernel>

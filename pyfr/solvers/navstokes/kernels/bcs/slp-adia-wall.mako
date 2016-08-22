# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.euler.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.common'/>
<%include file='pyfr.solvers.navstokes.kernels.flux'/>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur, ploc, t'>
    fpdtype_t nor = ${' + '.join('ul[{1}]*nl[{0}]'.format(i, i + 1)
                                 for i in range(ndims))};
    ur[0] = ul[0];
% for i in range(ndims):
    ur[${i + 1}] = ul[${i + 1}] - 2*nor*nl[${i}];
% endfor
    ur[${nvars - 1}] = ul[${nvars - 1}];
</%pyfr:macro>

<%pyfr:macro name='bc_common_flux_state' params='ul, gradul, fl, amul, nl, magnl, ploc, t'>
    // Ghost state r
    fpdtype_t ur[${nvars}];
    ${pyfr.expand('bc_ldg_state', 'ul', 'nl', 'ur', 'ploc', 't')};

    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}];
    ${pyfr.expand('rsolve', 'ul', 'ur', 'nl', 'ficomm')};

    fpdtype_t fln;
% for i in range(nvars):
    fln = ${'+'.join('nl[{0}]*fl[{0}][{1}]'.format(j, i)
                      for j in range(ndims))};

    ul[${i}] = magnl*(ficomm[${i}] - fln);
% endfor
</%pyfr:macro>

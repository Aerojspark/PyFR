# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>
<%include file='pyfr.solvers.navstokes.kernels.flux'/>

<%pyfr:kernel name='dflux' ndim='2'
              u='in fpdtype_t[${str(nvars)}]'
              amu='in fpdtype_t'
              f='inout fpdtype_t[${str(ndims)}][${str(nvars)}]'>
    // Compute the flux (F = Fi + Fv)
    fpdtype_t ftemp[${ndims}][${nvars}];
    fpdtype_t p, v[${ndims}];
    ${pyfr.expand('inviscid_flux', 'u', 'ftemp', 'p', 'v')};
    ${pyfr.expand('viscous_flux_add', 'u', 'f', 'amu', 'ftemp')};

    // Copy the fluxes
% for i, j in pyfr.ndrange(ndims, nvars):
    f[${i}][${j}] = ftemp[${i}][${j}];
% endfor
</%pyfr:kernel>

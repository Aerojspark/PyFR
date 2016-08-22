# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>

<%pyfr:kernel name='dflux' ndim='2'
              u='in fpdtype_t[${str(nvars)}]'
              f='out fpdtype_t[${str(ndims)}][${str(nvars)}]'>
    // Compute the flux
    fpdtype_t ftemp[${ndims}][${nvars}];
    fpdtype_t p, v[${ndims}];
    ${pyfr.expand('inviscid_flux', 'u', 'ftemp', 'p', 'v')};

    // Copy the fluxes
% for i, j in pyfr.ndrange(ndims, nvars):
    f[${i}][${j}] = ftemp[${i}][${j}];
% endfor
</%pyfr:kernel>

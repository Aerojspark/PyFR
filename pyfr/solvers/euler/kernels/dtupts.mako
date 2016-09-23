# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='dtupts' ndim='2'
              u='in fpdtype_t[${str(nvars)}]'
              wv='out fpdtype_t'>
    // Compute the flux
    fpdtype_t invrho = 1.0/u[0], E = u[${nvars - 1}];

    // Compute the velocities
    fpdtype_t rhov[${ndims}], v[${ndims}];
% for i in range(ndims):
    rhov[${i}] = u[${i + 1}];
    v[${i}] = invrho*rhov[${i}];
% endfor

    // Compute the pressure
    fpdtype_t p = ${c['gamma'] - 1}
                *(E - 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)});

    // Maximum characteristic wave speed (c+ |v|)
    wv = sqrt(${c['gamma']}*p*invrho) + sqrt(${pyfr.dot('v[{i}]', i=ndims)});
</%pyfr:kernel>

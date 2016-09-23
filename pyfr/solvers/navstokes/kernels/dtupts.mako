# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% cv = max(4/3, c['gamma'])/c['Pr'] %>

<%pyfr:kernel name='dtupts' ndim='2'
              u='in fpdtype_t[${str(nvars)}]'
              le='in broadcast fpdtype_t'
              artvisc='in broadcast fpdtype_t'
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

% if visc_corr == 'sutherland':
    // Compute the temperature and viscosity
    fpdtype_t cpT = ${c['gamma']}*(invrho*E - 0.5*${pyfr.dot('v[{i}]', i=ndims)});
    fpdtype_t Trat = ${1/c['cpTref']}*cpT;
    fpdtype_t mu_c = ${c['mu']*(c['cpTref'] + c['cpTs'])}*Trat*sqrt(Trat)
                   / (cpT + ${c['cpTs']});
% else:
    fpdtype_t mu_c = ${c['mu']};
% endif

% if shock_capturing == 'artificial-viscosity':
    mu_c += artvisc;
% endif

    // Maximum characteristic wave speed (c+ |v| + max(4/3, gamma)mu/rho/le/Pr)
    wv = sqrt(${c['gamma']}*p*invrho) + sqrt(${pyfr.dot('v[{i}]', i=ndims)})
       + ${cv}*mu_c*invrho/le;
</%pyfr:kernel>

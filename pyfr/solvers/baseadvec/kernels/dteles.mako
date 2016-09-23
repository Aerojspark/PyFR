# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='dteles' ndim='1'
              dt='inout fpdtype_t[${str(nupts)}]'
              le='in fpdtype_t'
              cfl='in scalar fpdtype_t'>
// Computea maximum wave speed
% for i in range(nupts - 1):
    dt[0] = max(dt[0], dt[${i+1}]);
% endfor

// Compute allowable time step based on CFL condition
dt[0] = cfl*le/dt[0];
</%pyfr:kernel>

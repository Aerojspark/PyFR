# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='divdisf' ndim='2'
              f='in fpdtype_t[${str((ndims)**2)}][${str(nvars)}]'
              u='out fpdtype_t[${str(nvars)}]'
              tsmats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'>
% for i in range(nvars):
    u[${i}] = ${'+'.join('tsmats[{j}][{k}]*f[{jk}][{i}]'
                         .format(i=i, j=j, k=k, jk=ndims*j+k)
                         for j, k in pyfr.ndrange(ndims, ndims))};
% endfor
</%pyfr:kernel>

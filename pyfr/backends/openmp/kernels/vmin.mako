# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

fpdtype_t
vmin(int n, fpdtype_t *__restrict__ x)
{
    fpdtype_t out = x[0];

    #pragma omp parallel for reduction(min:out)
    for (int i = 0; i < n; i++)
        out = min(out, x[i]);

    return out;
}

# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

void
vmin(long *n, double *out,
       fpdtype_t **xp)
{
    fpdtype_t *x = *xp;

    fpdtype_t val = x[0];

    #pragma omp parallel for reduction(min:val)
    for (int i = 0; i < *n; i++)
        val = max(val, x[i]);

    *out = err;
}

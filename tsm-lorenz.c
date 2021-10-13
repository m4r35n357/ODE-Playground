/*
 * Lorenz System
 *
 * Example: ./tsm-lorenz-dbg 9 32 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t sigma, rho, beta, _; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 13);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->sigma, &p->rho, &p->beta, &p->_);
    mpfr_div(p->beta, p->beta, p->_, RND);
    return p;
}

void ode (components *v, series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = S(y - x)
    mpfr_fmms(v->x, p->sigma, y[k], p->sigma, x[k], RND);
    //  y' = x(R - z) - y
    mpfr_fms(v->y, x[k], p->rho, *t_prod(x, z, k), RND);
    mpfr_sub(v->y, v->y, y[k], RND);
    //  z' = xy - Bz
    mpfr_fms(v->z, p->beta, z[k], *t_prod(x, y, k), RND);
    mpfr_neg(v->z, v->z, RND);
}

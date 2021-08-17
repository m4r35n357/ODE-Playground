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

void *get_p (int argc, char **argv, long n) {
    assert(argc == 13);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->sigma, &p->rho, &p->beta, &p->_);
    mpfr_div(p->beta, p->beta, p->_, RND);
    return p;
}

void ode (series x, series y, series z, components *c, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = S(y - x)
    mpfr_fmms(c->x, p->sigma, y[k], p->sigma, x[k], RND);
    //  y' = x(R - z) - y
    mpfr_fms(c->y, x[k], p->rho, *t_prod(x, z, k), RND);
    mpfr_sub(c->y, c->y, y[k], RND);
    //  z' = xy - Bz
    mpfr_fms(c->z, p->beta, z[k], *t_prod(x, y, k), RND);
    mpfr_neg(c->z, c->z, RND);
}

/*
 * Halvorsen Cyclic Attractor
 *
 * Example: ./tsm-halvorsen-dbg 9 32 10 .01 10000 1 0 0 1.4 4
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t a, b, w4x, w4y, w4z, _; } parameters;

void *get_p (int argc, char **argv, long n) {
    assert(argc == 11);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b);
    mpfr_inits(p->w4x, p->w4y, p->w4z, p->_, NULL);
    return p;
}

void ode (series x, series y, series z, components *c, void *params, int k) {
    parameters *p = (parameters *)params;
    mpfr_mul(p->w4x, x[k], p->b, RND);
    mpfr_mul(p->w4y, y[k], p->b, RND);
    mpfr_mul(p->w4z, z[k], p->b, RND);
    //  x' = - Ax - 4y - 4z - y^2
    mpfr_fma(p->_, p->a, x[k], *t_sqr(y, k), RND);
    mpfr_add(p->_, p->w4y, p->_, RND);
    mpfr_add(p->_, p->w4z, p->_, RND);
    mpfr_neg(c->x, p->_, RND);
    //  y' = - Ay - 4z - 4x - z^2
    mpfr_fma(p->_, p->a, y[k], *t_sqr(z, k), RND);
    mpfr_add(p->_, p->w4z, p->_, RND);
    mpfr_add(p->_, p->w4x, p->_, RND);
    mpfr_neg(c->y, p->_, RND);
    //  z' = - Az - 4x - 4y - x^2
    mpfr_fma(p->_, p->a, z[k], *t_sqr(x, k), RND);
    mpfr_add(p->_, p->w4x, p->_, RND);
    mpfr_add(p->_, p->w4y, p->_, RND);
    mpfr_neg(c->z, p->_, RND);
}

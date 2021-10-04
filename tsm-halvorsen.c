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

typedef struct { mpfr_t a, b, _; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 11);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b);
    mpfr_init(p->_);
    return p;
}

void ode (series x, series y, series z, components *c, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = - Ax - By - Bz - y^2
    mpfr_fmma(p->_, p->b, y[k], p->b, z[k], RND);
    mpfr_fma(p->_, p->a, x[k], p->_, RND);
    mpfr_add(p->_, *t_sqr(y, k), p->_, RND);
    mpfr_neg(c->x, p->_, RND);
    //  y' = - Ay - Bz - Bx - z^2
    mpfr_fmma(p->_, p->b, z[k], p->b, x[k], RND);
    mpfr_fma(p->_, p->a, y[k], p->_, RND);
    mpfr_add(p->_, *t_sqr(z, k), p->_, RND);
    mpfr_neg(c->y, p->_, RND);
    //  z' = - Az - Bx - By - x^2
    mpfr_fmma(p->_, p->b, x[k], p->b, y[k], RND);
    mpfr_fma(p->_, p->a, z[k], p->_, RND);
    mpfr_add(p->_, *t_sqr(x, k), p->_, RND);
    mpfr_neg(c->z, p->_, RND);
}

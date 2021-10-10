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

typedef struct { mpfr_t a, b; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 11);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b);
    return p;
}

void ode (series x, series y, series z, components *c, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = - Ax - By - Bz - y^2
    mpfr_fmma(c->x, p->b, y[k], p->b, z[k], RND);
    mpfr_fma(c->x, p->a, x[k], c->x, RND);
    mpfr_add(c->x, *t_sqr(y, k), c->x, RND);
    mpfr_neg(c->x, c->x, RND);
    //  y' = - Ay - Bz - Bx - z^2
    mpfr_fmma(c->y, p->b, z[k], p->b, x[k], RND);
    mpfr_fma(c->y, p->a, y[k], c->y, RND);
    mpfr_add(c->y, *t_sqr(z, k), c->y, RND);
    mpfr_neg(c->y, c->y, RND);
    //  z' = - Az - Bx - By - x^2
    mpfr_fmma(c->z, p->b, x[k], p->b, y[k], RND);
    mpfr_fma(c->z, p->a, z[k], c->z, RND);
    mpfr_add(c->z, *t_sqr(x, k), c->z, RND);
    mpfr_neg(c->z, c->z, RND);
}

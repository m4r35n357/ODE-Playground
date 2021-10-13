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

void ode (components *v, series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = - Ax - By - Bz - y^2
    mpfr_fmma(v->x, p->b, y[k], p->b, z[k], RND);
    mpfr_fma(v->x, p->a, x[k], v->x, RND);
    mpfr_add(v->x, *t_sqr(y, k), v->x, RND);
    mpfr_neg(v->x, v->x, RND);
    //  y' = - Ay - Bz - Bx - z^2
    mpfr_fmma(v->y, p->b, z[k], p->b, x[k], RND);
    mpfr_fma(v->y, p->a, y[k], v->y, RND);
    mpfr_add(v->y, *t_sqr(z, k), v->y, RND);
    mpfr_neg(v->y, v->y, RND);
    //  z' = - Az - Bx - By - x^2
    mpfr_fmma(v->z, p->b, x[k], p->b, y[k], RND);
    mpfr_fma(v->z, p->a, z[k], v->z, RND);
    mpfr_add(v->z, *t_sqr(x, k), v->z, RND);
    mpfr_neg(v->z, v->z, RND);
}

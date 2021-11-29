/*
 * Rossler System
 *
 * Example: ./tsm-rossler-dbg 9 32 10 0.01 50000 0.0 -6.78 0.02 .2 .2 5.7
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t a, b, c; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 12);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->c);
    return p;
}

void ode (components *v, series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = - y - z
    mpfr_add(v->x, y[k], z[k], RND);
    mpfr_neg(v->x, v->x, RND);
    //  y' = x + Ay
    mpfr_fma(v->y, p->a, y[k], x[k], RND);
    //  z' = B + z(x - C)
    mpfr_fms(v->z, p->c, z[k], *t_mul(z, x, k), RND);
    mpfr_sub(v->z, *t_const(p->b, k), v->z, RND);
}

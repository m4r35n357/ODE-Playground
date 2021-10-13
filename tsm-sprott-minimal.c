/*
 * Sprott Minimal System
 *
 * Example: ./tsm-sprott-minimal-dbg 9 32 4 0.1 10000 .02 0 0 2.017
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t a; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 10);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a);
    return p;
}

void ode (components *v, series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = y
    mpfr_set(v->x, y[k], RND);
    //  y' = z
    mpfr_set(v->y, z[k], RND);
    //  z' = - az + y^2 - x
    mpfr_fma(v->z, p->a, z[k], x[k], RND);
    mpfr_sub(v->z, *t_sqr(y, k), v->z, RND);
}

/*
 * Genesio-Tesi System
 *
 * Example: ./tsm-genesio-tesi-dbg 9 32 10 0.01 150000 .1 .1 .1 .44 1.1
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
    //  x' = y
    mpfr_set(c->x, y[k], RND);
    //  y' = z
    mpfr_set(c->y, z[k], RND);
    //  z' = - Az - By - x(1 + x)
    mpfr_fmma(c->z, p->a, z[k], p->b, y[k], RND);
    mpfr_add(c->z, c->z, x[k], RND);
    mpfr_add(c->z, c->z, *t_sqr(x, k), RND);
    mpfr_neg(c->z, c->z, RND);
}

/*
 * Rossler System
 *
 * Example: ./tsm-rossler-dbg 9 32 10 0.01 50000 0.0 -6.78 0.02 .2 .2 5.7
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t a, b, c, d0; } parameters;

void *get_p (int argc, char **argv, long n) {
    assert(argc == 12);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->c);
    mpfr_init_set_ui(p->d0, 0, RND);
    return p;
}

void ode (series x, series y, series z, components *c, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = - y - z
    mpfr_add(c->x, y[k], z[k], RND);
    mpfr_neg(c->x, c->x, RND);
    //  y' = x + Ay
    mpfr_fma(c->y, p->a, y[k], x[k], RND);
    //  z' = B + z(x - C)
    mpfr_fms(c->z, p->c, z[k], *t_prod(z, x, k), RND);
    mpfr_sub(c->z, *t_const(&p->b, &p->d0, k), c->z, RND);
}

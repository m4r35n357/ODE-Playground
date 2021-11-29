/*
 * Rucklidge Attractor
 *
 * Example: ./tsm-rucklidge-dbg 9 32 10 0.01 15000 1 0 0 6.7 2
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t alpha, kappa; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 11);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->alpha, &p->kappa);
    return p;
}

void ode (components *v, series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = ay - kx - yz
    mpfr_fmms(v->x, p->alpha, y[k], p->kappa, x[k], RND);
    mpfr_sub(v->x, v->x, *t_mul(y, z, k), RND);
    //  y' = x
    mpfr_set(v->y, x[k], RND);
    //  z' = y^2 - z
    mpfr_sub(v->z, *t_sqr(y, k), z[k], RND);
}

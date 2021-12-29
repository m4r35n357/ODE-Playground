/*
 * Sprott-Jafari System
 *
 * Example: ./tsm-sj-dbg 9 32 10 .01 10000 0 3.9 .7 8.888 4
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
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

void ode (components *vk, series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = y
    mpfr_set(vk->x, y[k], RND);
    //  y' = yz - x
    mpfr_sub(vk->y, *t_mul(y, z, k), x[k], RND);
    //  z' = z - ax^2 - y^2 - b
    mpfr_fma(vk->z, p->a, *t_sqr(x, k), *t_sqr(y, k), RND);
    mpfr_sub(vk->z, z[k], vk->z, RND);
    mpfr_sub(vk->z, vk->z, *t_const(p->b, k), RND);
}

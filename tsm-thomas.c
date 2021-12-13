/*
 * Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-thomas-dbg 9 32 10 0.1 30000 1 0 0 .19
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t b; series sx, sy, sz, cx, cy, cz; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 10);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->b);
    p->sx = t_jet(n); p->cx = t_jet(n);
    p->sy = t_jet(n); p->cy = t_jet(n);
    p->sz = t_jet(n); p->cz = t_jet(n);
    return p;
}

void ode (components *vk, series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = sin(y) - Bx
    mpfr_fms(vk->x, p->b, x[k], *t_sin_cos(p->sy, p->cy, y, k, TRIG).a, RND);
    mpfr_neg(vk->x, vk->x, RND);
    //  y' = sin(z) - By
    mpfr_fms(vk->y, p->b, y[k], *t_sin_cos(p->sz, p->cz, z, k, TRIG).a, RND);
    mpfr_neg(vk->y, vk->y, RND);
    //  z' = sin(x) - Bz
    mpfr_fms(vk->z, p->b, z[k], *t_sin_cos(p->sx, p->cx, x, k, TRIG).a, RND);
    mpfr_neg(vk->z, vk->z, RND);
}

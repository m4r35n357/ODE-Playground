/*
 * Yu-Wang System
 *
 * Example: ./tsm-yu-wang-dbg 9 32 10 .001 50000 1 0 0 10 40 2 2.5
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t a, b, c, d; series xy, e_xy; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 13);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->c, &p->d);
    p->xy = t_jet(n); p->e_xy = t_jet(n);
    return p;
}

void ode (components *vk, series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = A(y - x)
    mpfr_fmms(vk->x, p->a, y[k], p->a, x[k], RND);
    //  y' = Bx - cxz
    mpfr_fmms(vk->y, p->b, x[k], p->c, *t_mul(x, z, k), RND);
    //  z' = e^(xy) - Dz
    mpfr_swap(p->xy[k], *t_mul(x, y, k));
    mpfr_fms(vk->z, p->d, z[k], *t_exp(p->e_xy, p->xy, k), RND);
    mpfr_neg(vk->z, vk->z, RND);
}

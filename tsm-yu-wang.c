/*
 * Yu-Wang System
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { mpfr_t a, b, c, d; series xy, e_xy; };

void *get_p (int argc, char **argv, int n) {
    assert(argc == 13);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->c, &p->d);
    p->xy = t_jet(n); p->e_xy = t_jet(n);
    return p;
}

void ode (components *vk, series x, series y, series z, parameters *p, int k) {
    //  x' = A(y - x)
    mpfr_fmms(vk->x, p->a, y[k], p->a, x[k], RND);
    //  y' = Bx - cxz
    mpfr_fmms(vk->y, p->b, x[k], p->c, *t_mul(x, z, k), RND);
    //  z' = e^(xy) - Dz
    mpfr_swap(p->xy[k], *t_mul(x, y, k));
    mpfr_fms(vk->z, p->d, z[k], *t_exp(p->e_xy, p->xy, k), RND);
    mpfr_neg(vk->z, vk->z, RND);
}

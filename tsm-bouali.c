/*
 * Bouali Attractor
 *
 * Example: ./tsm-bouali-dbg 9 32 10 0.01 50000 1 1 0 3 2.2 1 .01
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { mpfr_t a, b, g, m; series gx2; };

void *get_p (int argc, char **argv, int n) {
    assert(argc == 13);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->g, &p->m);
    p->gx2 = t_jet(n);
    return p;
}

void ode (components *vk, series x, series y, series z, parameters *p, int k) {
    //  x' = Ax(1 - y) - Bz
    mpfr_fmms(vk->x, p->a, x[k], p->a, *t_mul(x, y, k), RND);
    mpfr_fms(vk->x, p->b, z[k], vk->x, RND);
    mpfr_neg(vk->x, vk->x, RND);
    //  y' = - Gy(1 - x^2)
    mpfr_mul(p->gx2[k], p->g, *t_sqr(x, k), RND);
    mpfr_fms(vk->y, p->g, y[k], *t_mul(y, p->gx2, k), RND);
    mpfr_neg(vk->y, vk->y, RND);
    //  z' = Mx
    mpfr_mul(vk->z, p->m, x[k], RND);
}

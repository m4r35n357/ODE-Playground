/*
 * Halvorsen Cyclic Attractor
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { mpfr_t a, b; };

parameters *get_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 11);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b);
    return p;
}

void ode (components *vk, series x, series y, series z, parameters *p, int k) {
    //  x' = - Ax - By - Bz - y^2
    mpfr_fmma(vk->x, p->b, y[k], p->b, z[k], RND);
    mpfr_fma(vk->x, p->a, x[k], vk->x, RND);
    mpfr_add(vk->x, *t_sqr(y, k), vk->x, RND);
    mpfr_neg(vk->x, vk->x, RND);
    //  y' = - Ay - Bz - Bx - z^2
    mpfr_fmma(vk->y, p->b, z[k], p->b, x[k], RND);
    mpfr_fma(vk->y, p->a, y[k], vk->y, RND);
    mpfr_add(vk->y, *t_sqr(z, k), vk->y, RND);
    mpfr_neg(vk->y, vk->y, RND);
    //  z' = - Az - Bx - By - x^2
    mpfr_fmma(vk->z, p->b, x[k], p->b, y[k], RND);
    mpfr_fma(vk->z, p->a, z[k], vk->z, RND);
    mpfr_add(vk->z, *t_sqr(x, k), vk->z, RND);
    mpfr_neg(vk->z, vk->z, RND);
}

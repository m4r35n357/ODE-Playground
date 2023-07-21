/*
 * Genesio-Tesi System
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { mpfr_t a, b; };

parameters *get_p (int argc, char **argv, int n) { (void)n;
    assert(argc == 11);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b);
    return p;
}

void ode (components *vk, series x, series y, series z, parameters *p, int k) {
    //  x' = y
    mpfr_set(vk->x, y[k], RND);
    //  y' = z
    mpfr_set(vk->y, z[k], RND);
    //  z' = - Az - By - x(1 + x)
    mpfr_fmma(vk->z, p->a, z[k], p->b, y[k], RND);
    mpfr_add(vk->z, vk->z, x[k], RND);
    mpfr_add(vk->z, vk->z, *t_sqr(x, k), RND);
    mpfr_neg(vk->z, vk->z, RND);
}

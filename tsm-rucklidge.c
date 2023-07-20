/*
 * Rucklidge Attractor
 *
 * Example: ./tsm-rucklidge-dbg 9 32 10 0.01 15000 1 0 0 6.7 2
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { mpfr_t alpha, kappa; };

void *get_p (int argc, char **argv, int n) { (void)n;
    assert(argc == 11);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->alpha, &p->kappa);
    return p;
}

void ode (components *vk, series x, series y, series z, parameters *p, int k) {
    //  x' = ay - kx - yz
    mpfr_fmms(vk->x, p->alpha, y[k], p->kappa, x[k], RND);
    mpfr_sub(vk->x, vk->x, *t_mul(y, z, k), RND);
    //  y' = x
    mpfr_set(vk->y, x[k], RND);
    //  z' = y^2 - z
    mpfr_sub(vk->z, *t_sqr(y, k), z[k], RND);
}

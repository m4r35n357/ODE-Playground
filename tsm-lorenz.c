/*
 * Lorenz System
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { mpfr_t sigma, rho, beta, _; };

parameters *get_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 13);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->sigma, &p->rho, &p->beta, &p->_);
    mpfr_div(p->beta, p->beta, p->_, RND);
    return p;
}

void ode (components *vk, series x, series y, series z, parameters *p, int k) {
    //  x' = S(y - x)
    mpfr_fmms(vk->x, p->sigma, y[k], p->sigma, x[k], RND);
    //  y' = x(R - z) - y
    mpfr_fms(vk->y, x[k], p->rho, *t_mul(x, z, k), RND);
    mpfr_sub(vk->y, vk->y, y[k], RND);
    //  z' = xy - Bz
    mpfr_fms(vk->z, p->beta, z[k], *t_mul(x, y, k), RND);
    mpfr_neg(vk->z, vk->z, RND);
}

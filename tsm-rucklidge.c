/*
 * Rucklidge Attractor
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { real alpha, kappa; };

model *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 11);
    model *_ = malloc(sizeof (model)); CHECK(_);
    tsm_get_p(argv, argc, &_->alpha, &_->kappa);
    return _;
}

void ode (triplet *v, series x, series y, series z, model *_, int k) {
    //  x' = ay - kx - yz
    mpfr_fmms(v->x, _->alpha, y[k], _->kappa, x[k], RND);
    mpfr_sub(v->x, v->x, *t_mul(y, z, k), RND);
    //  y' = x
    mpfr_set(v->y, x[k], RND);
    //  z' = y^2 - z
    mpfr_sub(v->z, *t_sqr(y, k), z[k], RND);
}

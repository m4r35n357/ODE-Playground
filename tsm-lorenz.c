/*
 * Lorenz System
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { real sigma, rho, beta, d; };

model *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 13);
    model *_ = malloc(sizeof (model)); CHECK(_);
    tsm_get_p(argv, argc, &_->sigma, &_->rho, &_->beta, &_->d);
    mpfr_div(_->beta, _->beta, _->d, RND);
    return _;
}

void ode (triplet *v, series x, series y, series z, model *_, int k) {
    //  x' = S(y - x)
    mpfr_fmms(v->x, _->sigma, y[k], _->sigma, x[k], RND);
    //  y' = x(R - z) - y
    mpfr_fms(v->y, x[k], _->rho, *t_mul(x, z, k), RND);
    mpfr_sub(v->y, v->y, y[k], RND);
    //  z' = xy - Bz
    mpfr_fms(v->z, _->beta, z[k], *t_mul(x, y, k), RND);
    mpfr_neg(v->z, v->z, RND);
}

/*
 * Bouali Attractor
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { real a, b, g, m; series gx2; };

model *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 13);
    model *_ = malloc(sizeof (model)); CHECK(_);
    tsm_get_p(argv, argc, &_->a, &_->b, &_->g, &_->m);
    _->gx2 = tsm_jet(n);
    return _;
}

void ode (triplet *v, series x, series y, series z, model *_, int k) {
    //  x' = Ax(1 - y) - Bz
    mpfr_fmms(v->x, _->a, x[k], _->a, *t_mul(x, y, k), RND);
    mpfr_fms(v->x, _->b, z[k], v->x, RND);
    mpfr_neg(v->x, v->x, RND);
    //  y' = - Gy(1 - x^2)
    mpfr_mul(_->gx2[k], _->g, *t_sqr(x, k), RND);
    mpfr_fms(v->y, _->g, y[k], *t_mul(y, _->gx2, k), RND);
    mpfr_neg(v->y, v->y, RND);
    //  z' = Mx
    mpfr_mul(v->z, _->m, x[k], RND);
}

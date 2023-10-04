/*
 * Halvorsen Cyclic Attractor
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { real a, D4; };

model *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 10);
    model *_ = malloc(sizeof (model)); CHECK(_);
    tsm_get_p(argv, argc, &_->a);
    mpfr_init_set_si(_->D4, 4, RND);
    return _;
}

void ode (triplet *v, series x, series y, series z, model *_, int k) {
    //  x' = - Ax - 4y - 4z - y^2
    mpfr_fmma(v->x, _->D4, y[k], _->D4, z[k], RND);
    mpfr_fma(v->x, _->a, x[k], v->x, RND);
    mpfr_add(v->x, *t_sqr(y, k), v->x, RND);
    mpfr_neg(v->x, v->x, RND);
    //  y' = - Ay - 4z - 4x - z^2
    mpfr_fmma(v->y, _->D4, z[k], _->D4, x[k], RND);
    mpfr_fma(v->y, _->a, y[k], v->y, RND);
    mpfr_add(v->y, *t_sqr(z, k), v->y, RND);
    mpfr_neg(v->y, v->y, RND);
    //  z' = - Az - 4x - 4y - x^2
    mpfr_fmma(v->z, _->D4, x[k], _->D4, y[k], RND);
    mpfr_fma(v->z, _->a, z[k], v->z, RND);
    mpfr_add(v->z, *t_sqr(x, k), v->z, RND);
    mpfr_neg(v->z, v->z, RND);
}

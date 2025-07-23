/*
 * Rossler System
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { real a, b, c; };

model *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 12);
    model *_ = malloc(sizeof (model)); CHECK(_);
    tsm_get_p(argv, argc, &_->a, &_->b, &_->c);
    return _;
}

void ode (triplet *v, series x, series y, series z, const model *_, const int k) {
    //  x' = - y - z
    mpfr_add(v->x, y[k], z[k], RND);
    mpfr_neg(v->x, v->x, RND);
    //  y' = x + Ay
    mpfr_fma(v->y, _->a, y[k], x[k], RND);
    //  z' = B + z(x - C)
    mpfr_fms(v->z, _->c, z[k], *t_mul(z, x, k), RND);
    mpfr_neg(v->z, v->z, RND);
    if (!k) mpfr_add(v->z, v->z, _->b, RND);
}

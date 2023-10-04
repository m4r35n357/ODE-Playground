/*
 * Genesio-Tesi System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { real a, b; };

model *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 11);
    model *_ = malloc(sizeof (model)); CHECK(_);
    tsm_get_p(argv, argc, &_->a, &_->b);
    return _;
}

void ode (triplet *v, series x, series y, series z, model *_, int k) {
    //  x' = y
    mpfr_set(v->x, y[k], RND);
    //  y' = z
    mpfr_set(v->y, z[k], RND);
    //  z' = - Az - By - x(1 + x)
    mpfr_fmma(v->z, _->a, z[k], _->b, y[k], RND);
    mpfr_add(v->z, v->z, x[k], RND);
    mpfr_add(v->z, v->z, *t_sqr(x, k), RND);
    mpfr_neg(v->z, v->z, RND);
}

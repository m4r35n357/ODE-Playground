/*
 * Burke & Shaw System - http://www.atomosyd.net/spip.php?article33
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { real s, v; };

model *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 11);
    model *_ = malloc(sizeof (model)); CHECK(_);
    tsm_get_p(argv, argc, &_->s, &_->v);
    return _;
}

void ode (triplet *v, series x, series y, series z, const model *_, const int k) {
    //  x' = - S(x + y)
    mpfr_fmma(v->x, _->s, x[k], _->s, y[k], RND);
    mpfr_neg(v->x, v->x, RND);
    //  y' = - (Sxz + y)
    mpfr_fma(v->y, _->s, *t_mul(x, z, k), y[k], RND);
    mpfr_neg(v->y, v->y, RND);
    //  z' = Sxy + V
    mpfr_mul(v->z, _->s, *t_mul(x, y, k), RND);
    if (!k) mpfr_add(v->z, v->z, _->v, RND);
}

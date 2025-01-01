/*
 * Wimol-Banlue System
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { real a; series tx, sx; };

model *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 10);
    model *_ = malloc(sizeof (model)); CHECK(_);
    tsm_get_p(argv, argc, &_->a);
    _->tx = tsm_jet(n);
    _->sx = tsm_jet(n);
    return _;
}

void ode (triplet *v, series x, series y, series z, const model *_, int k) {
    //  x' = y - x
    mpfr_sub(v->x, y[k], x[k], RND);
    //  y' = - z * tan(x)
    t_tan_sec2(_->tx, _->sx, x, k, false);
    mpfr_neg(v->y, *t_mul(z, _->tx, k), RND);
    //  z' = - A + xy + |y|
    mpfr_add(v->z, *t_mul(x, y, k), *t_abs(y, k), RND);
    if (!k) mpfr_sub(v->z, v->z, _->a, RND);
}

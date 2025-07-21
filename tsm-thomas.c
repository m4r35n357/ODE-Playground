/*
 * Thomas' cyclically symmetric attractor
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { real b; series sx, sy, sz, cx, cy, cz; };

model *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 10);
    model *_ = malloc(sizeof (model)); CHECK(_);
    tsm_get_p(argv, argc, &_->b);
    _->sx = tsm_jet(n); _->cx = tsm_jet(n);
    _->sy = tsm_jet(n); _->cy = tsm_jet(n);
    _->sz = tsm_jet(n); _->cz = tsm_jet(n);
    return _;
}

void ode (triplet *v, series x, series y, series z, const model *_, int k) {
    //  x' = sin(y) - Bx
    mpfr_fms(v->x, _->b, x[k], *t_sin_cos(_->sy, _->cy, y, k, true).a, RND);
    mpfr_neg(v->x, v->x, RND);
    //  y' = sin(z) - By
    mpfr_fms(v->y, _->b, y[k], *t_sin_cos(_->sz, _->cz, z, k, true).a, RND);
    mpfr_neg(v->y, v->y, RND);
    //  z' = sin(x) - Bz
    mpfr_fms(v->z, _->b, z[k], *t_sin_cos(_->sx, _->cx, x, k, true).a, RND);
    mpfr_neg(v->z, v->z, RND);
}

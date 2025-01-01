/*
 * Yu-Wang System
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { real a, b, c, d; series xy, e_xy; };

model *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 13);
    model *_ = malloc(sizeof (model)); CHECK(_);
    tsm_get_p(argv, argc, &_->a, &_->b, &_->c, &_->d);
    _->xy = tsm_jet(n);
    _->e_xy = tsm_jet(n);
    return _;
}

void ode (triplet *v, series x, series y, series z, const model *_, int k) {
    //  x' = A(y - x)
    mpfr_fmms(v->x, _->a, y[k], _->a, x[k], RND);
    //  y' = Bx - cxz
    mpfr_fmms(v->y, _->b, x[k], _->c, *t_mul(x, z, k), RND);
    //  z' = e^(xy) - Dz
    mpfr_swap(_->xy[k], *t_mul(x, y, k));
    t_exp(_->e_xy, _->xy, k);
    mpfr_fms(v->z, _->d, z[k], _->e_xy[k], RND);
    mpfr_neg(v->z, v->z, RND);
}

/*
 * Rabinovich–Fabrikant System
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { real alpha, gamma, D1, D4; series a, b, c; };

model *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 11);
    model *_ = malloc(sizeof (model)); CHECK(_);
    tsm_get_p(argv, argc, &_->alpha, &_->gamma);
    mpfr_init_set_si(_->D1, 1, RND);
    mpfr_init_set_si(_->D4, 4, RND);
    _->a = tsm_jet(n);
    _->b = tsm_jet(n);
    _->c = tsm_jet(n);
    return _;
}

void ode (triplet *v, series x, series y, series z, const model *_, int k) {
    //  x' = y(z - 1 + x^2) + Gx
    mpfr_set(_->a[k], *t_sqr(x, k), RND);
    if (!k) mpfr_sub(_->a[k], *t_sqr(x, k), _->D1, RND);
    mpfr_add(_->a[k], z[k], _->a[k], RND);
    mpfr_fma(v->x, _->gamma, x[k], *t_mul(y, _->a, k), RND);
    //  y' = x(3z + 1 - x^2) + Gy
    mpfr_fms(_->b[k], _->D4, z[k], _->a[k], RND);
    mpfr_fma(v->y, _->gamma, y[k], *t_mul(x, _->b, k), RND);
    //  z' = -2z(A + xy)
    mpfr_set(_->c[k], *t_mul(x, y, k), RND);
    if (!k) mpfr_add(_->c[k], *t_mul(x, y, k), _->alpha, RND);
    mpfr_mul_2si(v->z, *t_mul(z, _->c, k), 1, RND);
    mpfr_neg(v->z, v->z, RND);
}

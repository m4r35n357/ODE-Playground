/*
 * Rabinovich–Fabrikant System
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real alpha, gamma; series _a, _b, _c; };

model *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 10);
    model *_ = malloc(sizeof (model)); CHECK(_);
    _->_a = tsm_jet(n);
    _->_b = tsm_jet(n);
    _->_c = tsm_jet(n);
    tsm_get_p(argv, argc, &_->alpha, &_->gamma);
    return _;
}

triplet ode (series x, series y, series z, const model *_, const int k) {
    _->_a[k] = z[k] + t_sqr(x, k) - t_const(1.0L, k);
    _->_b[k] = 4.0L * z[k] - _->_a[k];
    _->_c[k] = t_const(_->alpha, k) + t_mul(x, y, k);
    return (triplet) {
        .x = t_mul(y, _->_a, k) + _->gamma * x[k],
        .y = t_mul(x, _->_b, k) + _->gamma * y[k],
        .z = - 2.0L * t_mul(z, _->_c, k)
    };
}

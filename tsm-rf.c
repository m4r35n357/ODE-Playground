/*
 * Rabinovich–Fabrikant System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real alpha, gamma; series a, b, c, _ALPHA, _1; };

parameters *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 10);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    tsm_get_p(argv, argc, &_->alpha, &_->gamma);
    _->a = tsm_var(n);
    _->b = tsm_var(n);
    _->c = tsm_var(n);
    _->_ALPHA = tsm_const(n, _->alpha);
    _->_1 = tsm_const(n, 1.0L);
    return _;
}

triplet ode (series x, series y, series z, parameters *_, int k) {
    _->a[k] = z[k] + t_sqr(x, k) - _->_1[k];
    _->b[k] = 4.0L * z[k] - _->a[k];
    _->c[k] = _->_ALPHA[k] + t_mul(x, y, k);
    return (triplet) {
        .x = t_mul(y, _->a, k) + _->gamma * x[k],
        .y = t_mul(x, _->b, k) + _->gamma * y[k],
        .z = - 2.0L * t_mul(z, _->c, k)
    };
}

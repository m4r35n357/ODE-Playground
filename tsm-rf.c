/*
 * Rabinovich–Fabrikant System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a, g; series sa, sb, sc, _A, _1; };

model *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 10);
    model *_ = malloc(sizeof (model)); CHECK(_);
    tsm_get_p(argv, argc, &_->a, &_->g);
    _->sa = tsm_jet(n);
    _->sb = tsm_jet(n);
    _->sc = tsm_jet(n);
    _->_A = tsm_const(n, _->a);
    _->_1 = tsm_const(n, 1.0L);
    return _;
}

triplet ode (series x, series y, series z, model *_, int k) {
    _->sa[k] = z[k] + t_sqr(x, k) - _->_1[k];
    _->sb[k] = 4.0L * z[k] - _->sa[k];
    _->sc[k] = _->_A[k] + t_mul(x, y, k);
    return (triplet) {
        .x = t_mul(y, _->sa, k) + _->g * x[k],
        .y = t_mul(x, _->sb, k) + _->g * y[k],
        .z = - 2.0L * t_mul(z, _->sc, k)
    };
}

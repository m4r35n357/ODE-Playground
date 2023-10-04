/*
 * Wimol-Banlue System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a; series tx, sx, _A; };

model *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 9);
    model *_ = malloc(sizeof (model)); CHECK(_);
    tsm_get_p(argv, argc, &_->a);
    _->tx = tsm_jet(n);
    _->sx = tsm_jet(n);
    _->_A = tsm_const(n, _->a);
    return _;
}

triplet ode (series x, series y, series z, model *_, int k) {
    t_tan_sec2(_->tx, _->sx, x, k, false);
    return (triplet) {
        .x = y[k] - x[k],
        .y = - t_mul(z, _->tx, k),
        .z = - _->_A[k] + t_mul(x, y, k) + t_abs(y, k)
    };
}

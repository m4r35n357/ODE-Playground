/*
 * Wimol-Banlue System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a; series tx, s2x, _A; };

parameters *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 9);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    tsm_get_p(argv, argc, &_->a);
    _->tx = t_jet(n);
    _->s2x = t_jet(n);
    _->_A = t_const(n, _->a);
    return _;
}

triplet ode (series x, series y, series z, parameters *p, int k) {
    t_tan_sec2(p->tx, p->s2x, x, k, false);
    return (triplet) {
        .x = y[k] - x[k],
        .y = - t_mul(z, p->tx, k),
        .z = - p->_A[k] + t_mul(x, y, k) + t_abs(y, k)
    };
}

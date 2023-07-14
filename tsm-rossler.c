/*
 * Rossler System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a, b, c; series _B; };

parameters *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 11);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    tsm_get_p(argv, argc, &_->a, &_->b, &_->c);
    _->_B = t_const(n, _->b);
    return _;
}

triplet ode (series x, series y, series z, parameters *p, int k) {
    return (triplet) {
        .x = - y[k] - z[k],
        .y = x[k] + p->a * y[k],
        .z = p->_B[k] + t_mul(x, z, k) - p->c * z[k]
    };
}

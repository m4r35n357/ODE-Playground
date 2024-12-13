/*
 * Lorenz System
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real sigma, rho, beta; };

model *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 12);
    model *_ = malloc(sizeof (model)); CHECK(_);
    real d;
    tsm_get_p(argv, argc, &_->sigma, &_->rho, &_->beta, &d);
    _->beta /= d;
    return _;
}

triplet ode (series x, series y, series z, const model *_, int k) {
    return (triplet) {
        .x = _->sigma * (y[k] - x[k]),
        .y = _->rho * x[k] - y[k] - t_mul(x, z, k),
        .z = t_mul(x, y, k) - _->beta * z[k]
    };
}

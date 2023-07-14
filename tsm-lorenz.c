/*
 * Lorenz System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real sigma, rho, beta; };

parameters *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 12);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    real d;
    tsm_get_p(argv, argc, &_->sigma, &_->rho, &_->beta, &d);
    _->beta /= d;
    return _;
}

triplet ode (series x, series y, series z, parameters *p, int k) {
    return (triplet) {
        .x = p->sigma * (y[k] - x[k]),
        .y = p->rho * x[k] - y[k] - t_mul(x, z, k),
        .z = t_mul(x, y, k) - p->beta * z[k]
    };
}

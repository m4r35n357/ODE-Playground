/*
 * Rucklidge Attractor
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real alpha, kappa; };

parameters *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 10);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    tsm_get_p(argv, argc, &_->alpha, &_->kappa);
    return _;
}

triplet ode (series x, series y, series z, parameters *_,  int k) {
    return (triplet) {
        .x = _->alpha * y[k] - _->kappa * x[k] - t_mul(y, z, k),
        .y = x[k],
        .z = t_sqr(y, k) - z[k]
    };
}

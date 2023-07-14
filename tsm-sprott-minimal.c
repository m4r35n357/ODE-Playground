/*
 * Sprott Minimal System
 *
 * Example: ./tsm-sprott-minimal-std 15 10 0.01 10000 .02 0 0 2.017
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a; };

parameters *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 9);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    tsm_get_p(argv, argc, &_->a);
    return _;
}

triplet ode (series x, series y, series z, parameters *_, int k) {
    return (triplet) {
        .x = y[k],
        .y = z[k],
        .z = - _->a * z[k] + t_sqr(y, k) - x[k]
    };
}

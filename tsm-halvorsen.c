/*
 * Halvorsen Cyclic Attractor
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a; };

model *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 9);
    model *_ = malloc(sizeof (model)); CHECK(_);
    tsm_get_p(argv, argc, &_->a);
    return _;
}

triplet ode (series x, series y, series z, model *_, int k) {
    return (triplet) {
        .x = - _->a * x[k] - 4.0L * (y[k] + z[k]) - t_sqr(y, k),
        .y = - _->a * y[k] - 4.0L * (z[k] + x[k]) - t_sqr(z, k),
        .z = - _->a * z[k] - 4.0L * (x[k] + y[k]) - t_sqr(x, k)
    };
}

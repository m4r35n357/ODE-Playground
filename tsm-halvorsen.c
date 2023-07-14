/*
 * Halvorsen Cyclic Attractor
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

triplet ode (series x, series y, series z, parameters *p, int k) {
    return (triplet) {
        .x = - p->a * x[k] - 4.0L * (y[k] + z[k]) - t_sqr(y, k),
        .y = - p->a * y[k] - 4.0L * (z[k] + x[k]) - t_sqr(z, k),
        .z = - p->a * z[k] - 4.0L * (x[k] + y[k]) - t_sqr(x, k)
    };
}

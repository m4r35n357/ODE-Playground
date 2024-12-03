/*
 * Rossler System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a, b, c; };

model *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 11);
    model *_ = malloc(sizeof (model)); CHECK(_);
    tsm_get_p(argv, argc, &_->a, &_->b, &_->c);
    return _;
}

triplet ode (series x, series y, series z, const model *_, int k) {
    return (triplet) {
        .x = - y[k] - z[k],
        .y = x[k] + _->a * y[k],
        .z = t_const(_->b, k) + t_mul(x, z, k) - _->c * z[k]
    };
}

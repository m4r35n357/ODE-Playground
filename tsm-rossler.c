/*
 * Rossler System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a, c; series B; };

model *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 11);
    model *_ = malloc(sizeof (model)); CHECK(_);
    _->B = tsm_jet(n);
    tsm_get_p(argv, argc, &_->a, _->B, &_->c);
    return _;
}

triplet ode (series x, series y, series z, model *_, int k) {
    return (triplet) {
        .x = - y[k] - z[k],
        .y = x[k] + _->a * y[k],
        .z = _->B[k] + t_mul(x, z, k) - _->c * z[k]
    };
}

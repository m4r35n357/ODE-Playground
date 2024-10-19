/*
 * Rucklidge Attractor
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a, k; };

model *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 10);
    model *_ = malloc(sizeof (model)); CHECK(_);
    tsm_get_p(argv, argc, &_->a, &_->k);
    return _;
}

triplet ode (series x, series y, series z, const model *_,  int k) {
    return (triplet) {
        .x = _->a * y[k] - _->k * x[k] - t_mul(y, z, k),
        .y = x[k],
        .z = t_sqr(y, k) - z[k]
    };
}

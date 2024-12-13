/*
 * Genesio-Tesi System - http://www.atomosyd.net/spip.php?article153
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a, b; };

model *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 10);
    model *_ = malloc(sizeof (model)); CHECK(_);
    tsm_get_p(argv, argc, &_->a, &_->b);
    return _;
}

triplet ode (series x, series y, series z, const model *_, int k) {
    return (triplet) {
        .x = y[k],
        .y = z[k],
        .z = - t_sqr(x, k) - x[k] - _->b * y[k] - _->a * z[k]
    };
}

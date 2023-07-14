/*
 * Genesio-Tesi System - http://www.atomosyd.net/spip.php?article153
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a, b; };

parameters *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 10);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    tsm_get_p(argv, argc, &_->a, &_->b);
    return _;
}

triplet ode (series x, series y, series z, parameters *p, int k) {
    return (triplet) {
        .x = y[k],
        .y = z[k],
        .z = - t_sqr(x, k) - x[k] - p->b * y[k] - p->a * z[k]
    };
}

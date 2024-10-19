/*
 * Burke & Shaw System - http://www.atomosyd.net/spip.php?article33
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real s; series V; };

model *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 10);
    model *_ = malloc(sizeof (model)); CHECK(_);
    _->V = tsm_jet(n);
    tsm_get_p(argv, argc, &_->s, _->V);
    return _;
}

triplet ode (series x, series y, series z, const model *_, int k) {
    return (triplet) {
        .x = - _->s * (x[k] + y[k]),
        .y = - _->s * t_mul(x, z, k) - y[k],
        .z =   _->s * t_mul(x, y, k) + _->V[k]
    };
}

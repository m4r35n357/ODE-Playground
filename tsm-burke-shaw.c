/*
 * Burke & Shaw System - http://www.atomosyd.net/spip.php?article33
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real s, v; series _V; };

parameters *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 10);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    tsm_get_p(argv, argc, &_->s, &_->v);
    _->_V = tsm_const(n, _->v);
    return _;
}

triplet ode (series x, series y, series z, parameters *_, int k) {
    return (triplet) {
        .x = - _->s * (x[k] + y[k]),
        .y = - _->s * t_mul(x, z, k) - y[k],
        .z =   _->s * t_mul(x, y, k) + _->_V[k]
    };
}

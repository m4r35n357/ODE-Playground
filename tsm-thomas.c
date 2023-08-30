/*
 * Thomas' cyclically symmetric attractor
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real b; series sx, sy, sz, cx, cy, cz; };

parameters *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 9);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    tsm_get_p(argv, argc, &_->b);
    _->sx = tsm_var(n); _->cx = tsm_var(n);
    _->sy = tsm_var(n); _->cy = tsm_var(n);
    _->sz = tsm_var(n); _->cz = tsm_var(n);
    return _;
}

triplet ode (series x, series y, series z, parameters *_, int k) {
    return (triplet) {
        .x = t_sin_cos(_->sy, _->cy, y, k, true).a - _->b * x[k],
        .y = t_sin_cos(_->sz, _->cz, z, k, true).a - _->b * y[k],
        .z = t_sin_cos(_->sx, _->cx, x, k, true).a - _->b * z[k]
    };
}

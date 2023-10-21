/*
 * Thomas' cyclically symmetric attractor
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real b; series sx, sy, sz, cx, cy, cz; };

model *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 9);
    model *_ = malloc(sizeof (model)); CHECK(_);
    _->sx = tsm_jet(n); _->cx = tsm_jet(n);
    _->sy = tsm_jet(n); _->cy = tsm_jet(n);
    _->sz = tsm_jet(n); _->cz = tsm_jet(n);
    tsm_get_p(argv, argc, &_->b);
    return _;
}

triplet ode (series x, series y, series z, model *_, int k) {
    return (triplet) {
        .x = t_sin_cos(_->sy, _->cy, y, k, true).a - _->b * x[k],
        .y = t_sin_cos(_->sz, _->cz, z, k, true).a - _->b * y[k],
        .z = t_sin_cos(_->sx, _->cx, x, k, true).a - _->b * z[k]
    };
}

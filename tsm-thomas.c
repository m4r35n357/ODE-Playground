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
    _->sx = t_jet(n); _->cx = t_jet(n);
    _->sy = t_jet(n); _->cy = t_jet(n);
    _->sz = t_jet(n); _->cz = t_jet(n);
    return _;
}

triplet ode (series x, series y, series z, parameters *p, int k) {
    return (triplet) {
        .x = t_sin_cos(p->sy, p->cy, y, k, true).a - p->b * x[k],
        .y = t_sin_cos(p->sz, p->cz, z, k, true).a - p->b * y[k],
        .z = t_sin_cos(p->sx, p->cx, x, k, true).a - p->b * z[k]
    };
}

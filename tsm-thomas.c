/*
 * Thomas' cyclically symmetric attractor
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real b; series sx, sy, sz, cx, cy, cz; };

void *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 9);
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
    tsm_get_p(argv, argc, &p->b);
    p->sx = t_jet(n); p->cx = t_jet(n);
    p->sy = t_jet(n); p->cy = t_jet(n);
    p->sz = t_jet(n); p->cz = t_jet(n);
    return p;
}

triplet ode (series x, series y, series z, parameters *p, int k) {
    return (triplet) {
        .x = t_sin_cos(p->sy, p->cy, y, k, true).a - p->b * x[k],
        .y = t_sin_cos(p->sz, p->cz, z, k, true).a - p->b * y[k],
        .z = t_sin_cos(p->sx, p->cx, x, k, true).a - p->b * z[k]
    };
}

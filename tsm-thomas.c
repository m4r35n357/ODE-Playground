/*
 * Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-thomas-std  6 8  0.1 30000  1 0 0  .185
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

typedef struct Parameters { real b; series sx, sy, sz, cx, cy, cz; } parameters;

void *get_p (int argc, char **argv, int n) {
    CHECK(argc == 9);
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
    t_params(argv, argc, &p->b);
    p->sx = t_jet(n); p->cx = t_jet(n);
    p->sy = t_jet(n); p->cy = t_jet(n);
    p->sz = t_jet(n); p->cz = t_jet(n);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = t_sin_cos(p->sy, p->cy, y, k, true).a - p->b * x[k],
        .y = t_sin_cos(p->sz, p->cz, z, k, true).a - p->b * y[k],
        .z = t_sin_cos(p->sx, p->cx, x, k, true).a - p->b * z[k]
    };
}

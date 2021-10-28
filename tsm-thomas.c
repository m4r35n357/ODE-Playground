/*
 * Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-thomas-dbg 15 10 0.1 30000 1 0 0 .19
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct { real b; series sx, sy, sz, cx, cy, cz; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 9);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->b);
    p->sx = t_jet(n); p->cx = t_jet(n);
    p->sy = t_jet(n); p->cy = t_jet(n);
    p->sz = t_jet(n); p->cz = t_jet(n);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = t_sin_cos(p->sy, p->cy, y, k, TRIG).a - p->b * x[k],
        .y = t_sin_cos(p->sz, p->cz, z, k, TRIG).a - p->b * y[k],
        .z = t_sin_cos(p->sx, p->cx, x, k, TRIG).a - p->b * z[k]
    };
}

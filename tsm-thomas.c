/*
 * Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-thomas-dbg NA NA 10 0.1 30000 1 0 0 .19
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real b;
} parameters;

typedef struct {
    series sx;
    series sy;
    series sz;
    series cx;
    series cy;
    series cz;
} intermediates;

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    intermediates *i = (intermediates *)inters;
    return (components) {
        .x = t_sin_cos(i->sy, i->cy, y, k, TRIG).a - p->b * x[k],
        .y = t_sin_cos(i->sz, i->cz, z, k, TRIG).a - p->b * y[k],
        .z = t_sin_cos(i->sx, i->cx, x, k, TRIG).a - p->b * z[k]
    };
}

int main (int argc, char **argv) {
    long order, steps;
    real x0, y0, z0, stepsize;

    assert(argc == 10);
    t_stepper(argv, &order, &stepsize, &steps);
    parameters p;
    t_args(argv, argc, &x0, &y0, &z0, &p.b);
    intermediates i = (intermediates) {
        .sx = t_jet(order),
        .sy = t_jet(order),
        .sz = t_jet(order),
        .cx = t_jet(order),
        .cy = t_jet(order),
        .cz = t_jet(order)
    };

    taylor(order, steps, stepsize, x0, y0, z0, &p, &i, ode);
    return 0;
}

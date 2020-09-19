/*
 * Yu-Wang System
 *
 * Example: ./tsm-yu-wang-dbg NA NA 10 .001 50000 1 0 0 10 40 2 2.5
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real a;
    real b;
    real c;
    real d;
} parameters;

typedef struct {
    series xy;
    series e_xy;
} intermediates;

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    intermediates *i = (intermediates *)inters;
    i->xy[k] = t_prod(x, y, k);
    return (components) {
        .x = p->a * (y[k] - x[k]),
        .y = p->b * x[k] - p->c * t_prod(x, z, k),
        .z = t_exp(i->e_xy, i->xy, k) - p->d * z[k]
    };
}

int main (int argc, char **argv) {
    long order, steps;
    real x0, y0, z0, stepsize;

    assert(argc == 13);
    t_stepper(argv, &order, &stepsize, &steps);
    parameters p;
    t_args(argv, argc, &x0, &y0, &z0, &p.a, &p.b, &p.c, &p.d);
    intermediates i = (intermediates) {
        .xy = t_jet(order),
        .e_xy = t_jet(order)
    };

    taylor(order, steps, stepsize, x0, y0, z0, &p, &i, ode);
    return 0;
}

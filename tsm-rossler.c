/*
 * Rossler System
 *
 * Example: ./tsm-rossler-dbg NA NA 10 0.01 50000 0.0 -6.78 0.02 .2 .2 5.7
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real a;
    series b;
    real c;
} parameters;

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = - y[k] - z[k],
        .y = x[k] + p->a * y[k],
        .z = p->b[k] + t_prod(x, z, k) - p->c * z[k]
    };
}

int main (int argc, char **argv) {
    long order, steps;
    real x0, y0, z0, stepsize;

    assert(argc == 12);
    t_stepper(argv, &order, &stepsize, &steps);
    parameters p = (parameters) {
        .b = t_jet(order)
    };
    t_args(argv, argc, &x0, &y0, &z0, &p.a, p.b, &p.c);

    tsm(order, steps, stepsize, x0, y0, z0, &p, NULL, ode);
    return 0;
}

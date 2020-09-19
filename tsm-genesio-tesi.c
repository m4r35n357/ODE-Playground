/*
 * Genesio-Tesi System
 *
 * Example: ./tsm-genesio-tesi-dbg NA NA 10 0.01 50000 .1 .1 .1 .44 1.1 1
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
} parameters;

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = y[k],
        .y = z[k],
        .z = t_sqr(x, k) - p->c * x[k] - p->b * y[k] - p->a * z[k]
    };
}

int main (int argc, char **argv) {
    long order, steps;
    real stepsize, x0, y0, z0;

    assert(argc == 12);
    t_stepper(argv, &order, &stepsize, &steps);
    parameters p;
    t_args(argv, argc, &x0, &y0, &z0, &p.a);

    taylor(order, steps, stepsize, x0, y0, z0, &p, NULL, ode);
    return 0;
}

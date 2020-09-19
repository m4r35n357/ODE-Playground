/*
 * Sprott Minimal System
 *
 * Example: ./tsm-sprott-minimal-dbg NA NA 10 0.01 10000 .02 0 0 2.017
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real a;
} parameters;

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = y[k],
        .y = z[k],
        .z = - p->a * z[k] + t_sqr(y, k) - x[k]
    };
}

int main (int argc, char **argv) {
    long order, steps;
    real x0, y0, z0, stepsize;

    assert(argc == 10);
    t_stepper(argv, &order, &stepsize, &steps);
    parameters p;
    t_args(argv, argc, &x0, &y0, &z0, &p.a);

    taylor(order, steps, stepsize, x0, y0, z0, &p, NULL, ode);
    return 0;
}

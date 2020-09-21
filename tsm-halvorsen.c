/*
 * Halvorsen Cyclic Attractor
 *
 * Example: ./tsm-halvorsen-dbg NA NA 10 .01 10000 1 0 0 1.4
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
    (void)inters;
    return (components) {
        .x = - p->a * x[k] - 4.0 * (y[k] + z[k]) - t_sqr(y, k),
        .y = - p->a * y[k] - 4.0 * (z[k] + x[k]) - t_sqr(z, k),
        .z = - p->a * z[k] - 4.0 * (x[k] + y[k]) - t_sqr(x, k)
    };
}

int main (int argc, char **argv) {
    long order, steps;
    real stepsize, x0, y0, z0;

    assert(argc == 10);
    t_stepper(argv, &order, &stepsize, &steps);
    parameters p;
    t_args(argv, argc, &x0, &y0, &z0, &p.a);

    tsm(order, steps, stepsize, x0, y0, z0, &p, NULL, ode);
    return 0;
}

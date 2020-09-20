/*
 * Sprott-Jafari System
 *
 * Example: ./tsm-sj-dbg NA NA 10 .01 10000 0 3.9 .7 8.888 4
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real a;
    series b;
} parameters;

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = y[k],
        .y = - x[k] + t_prod(y, z, k),
        .z = z[k] + p->a * t_sqr(x, k) - t_sqr(y, k) - p->b[k]
    };
}

int main (int argc, char **argv) {
    long order, steps;
    real stepsize, x0, y0, z0;

    assert(argc == 11);
    t_stepper(argv, &order, &stepsize, &steps);
    parameters p = (parameters) {
        .b = t_jet(order)
    };
    t_args(argv, argc, &x0, &y0, &z0, &p.a, p.b);

    tsm(order, steps, stepsize, x0, y0, z0, &p, NULL, ode);
    return 0;
}

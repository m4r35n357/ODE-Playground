/*
 * Lorenz System
 *
 * Example: ./tsm-lorenz-dbg NA NA 16 .01 10000 -15.8 -17.48 35.64 10 28 8 3
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real sigma;
    real rho;
    real beta;
} parameters;

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = p->sigma * (y[k] - x[k]),
        .y = p->rho * x[k] - y[k] - t_prod(x, z, k),
        .z = t_prod(x, y, k) - p->beta * z[k]
    };
}

int main (int argc, char **argv) {
    long order, steps;
    real x0, y0, z0, stepsize, _;

    assert(argc == 13);
    t_stepper(argv, &order, &stepsize, &steps);
    parameters p;
    t_args(argv, argc, &x0, &y0, &z0, &p.sigma, &p.rho, &p.beta, &_);
    p.beta /= _;

    tsm(order, steps, stepsize, x0, y0, z0, &p, NULL, ode);
    return 0;
}

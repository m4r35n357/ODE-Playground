/*
 * Lorenz System
 *
 * Example: ./rk4-lorenz-dbg NA NA 1 .01 10000 -15.8 -17.48 35.64 10 28 8 3
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

static components ode (real x, real y, real z, void *params) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = p->sigma * (y - x),
        .y = p->rho * x - y - x * z,
        .z = x * y - p->beta * z
    };
}

int main (int argc, char **argv) {
    long nsteps, interval;
    real x0, y0, z0, h, _;
    parameters p;

    assert(argc == 13);
    t_stepper(argv, &interval, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &p.sigma, &p.rho, &p.beta, &_);
    p.beta /= _;

    rk4(interval, nsteps, h, x0, y0, z0, &p, ode);
    return 0;
}

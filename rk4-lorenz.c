/*
 * Lorenz System
 *
 * Example: ./rk4-lorenz-dbg 15 1 .01 10000 -15.8 -17.48 35.64 10 28 8 3
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "runge-kutta-ode.h"

typedef struct { real sigma, rho, beta; } parameters;

void *get_p (int argc, char **argv) {
    assert(argc == 12);
    parameters *p = malloc(sizeof (parameters));
    real _;
    t_params(argv, argc, &p->sigma, &p->rho, &p->beta, &_);
    p->beta /= _;
    return p;
}

components ode (real x, real y, real z, void *params) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = p->sigma * (y - x),
        .y = p->rho * x - y - x * z,
        .z = x * y - p->beta * z
    };
}

/*
 * Thomas' cyclically symmetric attractor
 *
 * Example: ./rk4-thomas-dbg 15 1 0.1 30000 1 0 0 .19
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "runge-kutta-ode.h"

typedef struct { real b; } parameters;

void *get_p (int argc, char **argv) {
    assert(argc == 9);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->b);
    return p;
}

components ode (real x, real y, real z, void *params) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = sinl(y) - p->b * x,
        .y = sinl(z) - p->b * y,
        .z = sinl(x) - p->b * z
    };
}

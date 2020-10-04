/*
 * Sprott-Thomas' cyclically symmetric attractor
 *
 * Example: ./rk4-sprott-thomas-dbg 15 _ 1 0.01 30000 1 0 0 4.75 1
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "taylor-ode.h"

typedef struct {
    real a;
    real b;
} parameters;

static void *get_p (int argc, char **argv) {
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b);
    return p;
}

static components ode (real x, real y, real z, void *params) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = sin(p->a * y) - p->b * tan(x),
        .y = sin(p->a * z) - p->b * tan(y),
        .z = sin(p->a * x) - p->b * tan(z)
    };
}

int main (int argc, char **argv) {
    assert(argc == 11);
    rk4(argc, argv, ode, get_p);
    return 0;
}
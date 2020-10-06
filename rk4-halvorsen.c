/*
 * Halvorsen Cyclic Attractor
 *
 * Example: ./rk4-halvorsen-dbg 15 _ 1 .01 10000 1 0 0 1.4 4
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
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
        .x = - p->a * x - p->b * (y + z) - y * y,
        .y = - p->a * y - p->b * (z + x) - z * z,
        .z = - p->a * z - p->b * (x + y) - x * x
    };
}

int main (int argc, char **argv) {
    assert(argc == 11);
    rk4(argc, argv, ode, get_p);
    return 0;
}

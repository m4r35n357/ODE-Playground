/*
 * Halvorsen Cyclic Attractor
 *
 * Example: ./rk4-halvorsen-dbg 15 NA 1 .01 10000 1 0 0 1.4
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real a;
} parameters;

static void *get_p (int argc, char **argv) {
    parameters *p = malloc(sizeof (parameters));
    t_args(argv, argc, &p->a);
    return p;
}

static components ode (real x, real y, real z, void *params) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = - p->a * x - 4.0 * (y + z) - y * y,
        .y = - p->a * y - 4.0 * (z + x) - z * z,
        .z = - p->a * z - 4.0 * (x + y) - x * x
    };
}

int main (int argc, char **argv) {
    assert(argc == 10);
    rk4(argc, argv, ode, get_p);
    return 0;
}

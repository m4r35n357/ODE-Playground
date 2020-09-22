/*
 * Sprott Minimal System
 *
 * Example: ./rk4-sprott-minimal-dbg NA NA 1 0.01 10000 .02 0 0 2.017
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real a;
} parameters;

static components ode (real x, real y, real z, void *params) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = y,
        .y = z,
        .z = - p->a * z + y * y - x
    };
}

static void *get_p (int argc, char **argv, long order) {
    (void)order;
    parameters *p = malloc(sizeof (parameters));
    t_args(argv, argc, &p->a);
    return p;
}

int main (int argc, char **argv) {
    assert(argc == 10);
    rk4(argc, argv, ode, get_p);
    return 0;
}

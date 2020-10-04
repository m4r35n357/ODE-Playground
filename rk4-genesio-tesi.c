/*
 * Genesio-Tesi System
 *
 * Example: ./rk4-genesio-tesi-dbg 15 NA 1 0.01 50000 .1 .1 .1 .44 1.1 1
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real a;
    real b;
    real c;
} parameters;

static void *get_p (int argc, char **argv) {
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->c);
    return p;
}

static components ode (real x, real y, real z, void *params) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = y,
        .y = z,
        .z = x * x - p->c * x - p->b * y - p->a * z
    };
}

int main (int argc, char **argv) {
    assert(argc == 12);
    rk4(argc, argv, ode, get_p);
    return 0;
}

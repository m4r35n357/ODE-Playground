/*
 * Bouali Attractor
 *
 * Example: ./rk4-bouali-dbg 15 NA 1 0.01 50000 1 1 0 3 2.2 1 .01
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
    real d;
} parameters;

static void *get_p (int argc, char **argv) {
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->c, &p->d);
    return p;
}

static components ode (real x, real y, real z, void *params) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = p->a * x * (1.0 - y) - p->b * z,
        .y = - p->c * y * (1.0 - x * x),
        .z = p->d * x
    };
}

int main (int argc, char **argv) {
    assert(argc == 13);
    rk4(argc, argv, ode, get_p);
    return 0;
}

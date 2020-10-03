/*
 * Yu-Wang System
 *
 * Example: ./rk4-yu-wang-dbg 15 NA 1 .001 50000 1 0 0 10 40 2 2.5
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
    real c;
    real d;
} parameters;

static void *get_p (int argc, char **argv) {
    parameters *p = malloc(sizeof (parameters));
    t_args(argv, argc, &p->a, &p->b, &p->c, &p->d);
    return p;
}

static components ode (real x, real y, real z, void *params) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = p->a * (y - x),
        .y = p->b * x - p->c * x * z,
        .z = exp(x * y) - p->d * z
    };
}

int main (int argc, char **argv) {
    assert(argc == 13);
    rk4(argc, argv, ode, get_p);
    return 0;
}

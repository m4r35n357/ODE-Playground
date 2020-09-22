/*
 * Wimol-Banlue System
 *
 * Example: ./rk4-wimol-banlue-dbg NA NA 1 0.1 10000 1 0 0 2.0
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "taylor-ode.h"

typedef struct {
    real a;
} parameters;

static components ode (real x, real y, real z, void *params) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = y - x,
        .y = - z * tan(x),
        .z = - p->a + x * y + abs(y)
    };
}

static void *get_p (int argc, char **argv, long order) {
    (void)order;
    parameters *p = malloc(sizeof (parameters));
    t_args(argv, argc, p->a);
    return p;
}

int main (int argc, char **argv) {
    assert(argc == 10);
    rk4(argc, argv, ode, get_p);
    return 0;
}

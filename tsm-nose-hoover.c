/*
 * Nose-Hoover System
 *
 * Example: ./tsm-nose-hoover-dbg 15 10 0.01 10000 1 0 0 6.0
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct { real a; } parameters;

void *get_p (int argc, char **argv, int n) { (void)n;
    assert(argc == 9);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = y[k],
        .y = t_mul(y, z, k) - x[k],
        .z = p->a - t_sqr(y, k)
    };
}

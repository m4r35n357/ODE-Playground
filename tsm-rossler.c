/*
 * Rossler System
 *
 * Example: ./tsm-rossler-dbg 15 10 0.01 50000 0.0 -6.78 0.02 .2 .2 5.7
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct { real a, b, c; } parameters;

void *get_p (int argc, char **argv, int n) { (void)n;
    assert(argc == 11);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->c);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = - y[k] - z[k],
        .y = x[k] + p->a * y[k],
        .z = t_const(p->b, k) + t_mul(x, z, k) - p->c * z[k]
    };
}

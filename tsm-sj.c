/*
 * Sprott-Jafari System
 *
 * Example: ./tsm-sj-dbg 9 10 .01 10000 0 3.9 .7 8.888 4
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct { real a, b; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 10);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = y[k],
        .y = t_prod(y, z, k) - x[k],
        .z = z[k] - p->a * t_sqr(x, k) - t_sqr(y, k) - t_const(p->b, k)
    };
}

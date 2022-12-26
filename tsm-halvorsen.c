/*
 * Halvorsen Cyclic Attractor
 *
 * Example: ./tsm-halvorsen-std  6 8  .01 10000  1 0 0  1.4 4
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct Parameters { real a, b; } parameters;

void *get_p (int argc, char **argv, int n) { (void)n;
    assert(argc == 10);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = - p->a * x[k] - p->b * (y[k] + z[k]) - t_sqr(y, k),
        .y = - p->a * y[k] - p->b * (z[k] + x[k]) - t_sqr(z, k),
        .z = - p->a * z[k] - p->b * (x[k] + y[k]) - t_sqr(x, k)
    };
}

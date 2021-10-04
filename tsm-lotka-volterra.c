/*
 * Lotka-Volterra (Predator-Prey) System
 *
 * Example: ./tsm-lotka-volterra-dbg 15 10 .01 10000 10 10 0 1 .5 .05 .02
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct { real a, b, c, d; } parameters;

void *get_p (int argc, char **argv, long order) {
    assert(argc == 12);
    (void)order;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->c, &p->d);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    (void)z;
    real xy = t_prod(x, y, k);
    return (components) {
        .x = p->a * x[k] - p->c * xy,
        .y = p->d * xy - p->b * y[k],
        .z = 0.0L
    };
}

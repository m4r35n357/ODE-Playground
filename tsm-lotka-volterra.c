/*
 * Lotka-Volterra (Predator-Prey) System
 *
 * Example: ./tsm-lotka-volterra-std 15 10 .01 10000 10 10 0 1 .5 .05 .02
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct Parameters { real a, b, c, d; } parameters;

void *get_p (int argc, char **argv, int n) { (void)n;
    assert(argc == 12);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->c, &p->d);
    return p;
}

components ode (series x, series y, series z, void *params, int k) { (void)z;
    parameters *p = (parameters *)params;
    real xy = t_mul(x, y, k);
    return (components) {
        .x = p->a * x[k] - p->c * xy,
        .y = p->d * xy - p->b * y[k]
    };
}

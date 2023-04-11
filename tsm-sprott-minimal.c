/*
 * Sprott Minimal System
 *
 * Example: ./tsm-sprott-minimal-std 15 10 0.01 10000 .02 0 0 2.017
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

typedef struct Parameters { real a; } parameters;

void *get_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 9);
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
    t_params(argv, argc, &p->a);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = y[k],
        .y = z[k],
        .z = - p->a * z[k] + t_sqr(y, k) - x[k]
    };
}

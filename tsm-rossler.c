/*
 * Rossler System
 *
 * Example: ./tsm-rossler-std 15 10 0.01 50000 0.0 -6.78 0.02 .2 .2 5.7
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

typedef struct Parameters { real a, b, c; series wb; } parameters;

void *get_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 11);
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
    t_params(argv, argc, &p->a, &p->b, &p->c);
    p->wb = t_const(n, p->b);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = - y[k] - z[k],
        .y = x[k] + p->a * y[k],
        .z = p->wb[k] + t_mul(x, z, k) - p->c * z[k]
    };
}

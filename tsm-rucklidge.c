/*
 * Rucklidge Attractor
 *
 * Example: ./tsm-rucklidge-std 15 10 0.01 15000 1 0 0 6.7 2
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

typedef struct Parameters { real alpha, kappa; } parameters;

void *get_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 10);
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
    t_params(argv, argc, &p->alpha, &p->kappa);
    return p;
}

components ode (series x, series y, series z, void *params,  int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = p->alpha * y[k] - p->kappa * x[k] - t_mul(y, z, k),
        .y = x[k],
        .z = t_sqr(y, k) - z[k]
    };
}

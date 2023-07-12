/*
 * Rucklidge Attractor
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

typedef struct Parameters { real alpha, kappa; } parameters;

void *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 10);
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
    tsm_get_p(argv, argc, &p->alpha, &p->kappa);
    return p;
}

triplet ode (series x, series y, series z, void *params,  int k) {
    parameters *p = (parameters *)params;
    return (triplet) {
        .x = p->alpha * y[k] - p->kappa * x[k] - t_mul(y, z, k),
        .y = x[k],
        .z = t_sqr(y, k) - z[k]
    };
}

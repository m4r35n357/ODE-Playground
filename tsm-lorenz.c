/*
 * Lorenz System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

typedef struct Parameters { real sigma, rho, beta; } parameters;

void *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 12);
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
    real _;
    tsm_get_p(argv, argc, &p->sigma, &p->rho, &p->beta, &_);
    p->beta /= _;
    return p;
}

triplet ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    return (triplet) {
        .x = p->sigma * (y[k] - x[k]),
        .y = p->rho * x[k] - y[k] - t_mul(x, z, k),
        .z = t_mul(x, y, k) - p->beta * z[k]
    };
}

/*
 * Lorenz System
 *
 * Example: ./tsm-lorenz-std  6 8 .01 10000  -15.8 -17.48 35.64  10 28 8 3
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

typedef struct Parameters { real sigma, rho, beta; } parameters;

void *get_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 12);
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
    real _;
    t_params(argv, argc, &p->sigma, &p->rho, &p->beta, &_);
    p->beta /= _;
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = p->sigma * (y[k] - x[k]),
        .y = p->rho * x[k] - y[k] - t_mul(x, z, k),
        .z = t_mul(x, y, k) - p->beta * z[k]
    };
}

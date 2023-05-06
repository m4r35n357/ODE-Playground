/*
 * Rabinovich–Fabrikant System
 *
 * Example: ./tsm-rf-std 15 10 .01 50000 .05 -.05 .3 .2873 .1
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

typedef struct Parameters { real alpha, gamma; series a, b, c, walpha, w1; } parameters;

void *get_p (int argc, char **argv, int n) {
    CHECK(argc == 10);
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
    t_params(argv, argc, &p->alpha, &p->gamma);
    p->a = t_jet(n);
    p->b = t_jet(n);
    p->c = t_jet(n);
    p->walpha = t_const(n, p->alpha);
    p->w1 = t_const(n, 1.0L);
    return p;
}

triplet ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    p->a[k] = z[k] + t_sqr(x, k) - p->w1[k];
    p->b[k] = 4.0L * z[k] - p->a[k];
    p->c[k] = p->walpha[k] + t_mul(x, y, k);
    return (triplet) {
        .x = t_mul(y, p->a, k) + p->gamma * x[k],
        .y = t_mul(x, p->b, k) + p->gamma * y[k],
        .z = - 2.0L * t_mul(z, p->c, k)
    };
}

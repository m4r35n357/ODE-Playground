/*
 * Rabinovich–Fabrikant System
 *
 * Example: ./tsm-rf-std 15 10 .01 50000 .05 -.05 .3 .2873 .1
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct Parameters { real alpha, gamma; series a, b, c; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 10);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->alpha, &p->gamma);
    p->a = t_jet(n);
    p->b = t_jet(n);
    p->c = t_jet(n);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    p->a[k] = z[k] + t_sqr(x, k) - (!k ? 1.0L : 0.0L);
    p->b[k] = 4.0L * z[k] - p->a[k];
    p->c[k] = (!k ? p->alpha : 0.0L) + t_mul(x, y, k);
    return (components) {
        .x = t_mul(y, p->a, k) + p->gamma * x[k],
        .y = t_mul(x, p->b, k) + p->gamma * y[k],
        .z = - 2.0L * t_mul(z, p->c, k)
    };
}

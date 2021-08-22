/*
 * Rabinovich–Fabrikant System
 *
 * Example: ./tsm-rf-dbg 15 10 .01 50000 .05 -.05 .3 .28713 .1
 *          ./tsm-rf-dbg 15 16 .01 50000 .05 -.05 .3 .105 .1 | ./plotPi3d.py
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real alpha;
    real gamma;
    series a;
    series b;
    series c;
} parameters;

void *get_p (int argc, char **argv, long order) {
    assert(argc == 10);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->alpha, &p->gamma);
    p->a = t_jet(order);
    p->b = t_jet(order);
    p->c = t_jet(order);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    p->a[k] = z[k] + t_prod(x, x, k) - (k == 0 ? 1.0L : 0.0L);
    p->b[k] = 4.0L * z[k] - p->a[k];
    p->c[k] = (k == 0 ? p->alpha : 0.0L) + t_prod(x, y, k);
    return (components) {
        .x = t_prod(y, p->a, k) + p->gamma * x[k],
        .y = t_prod(x, p->b, k) + p->gamma * y[k],
        .z = - 2.0L * t_prod(z, p->c, k)
    };
}

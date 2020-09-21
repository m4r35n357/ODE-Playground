/*
 * Rucklidge Attractor
 *
 * Example: ./tsm-rucklidge-dbg NA NA 10 0.01 15000 1 0 0 6.7 2
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real alpha;
    real kappa;
} parameters;

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    (void)inters;
    return (components) {
        .x = p->alpha * y[k] - p->kappa * x[k] - t_prod(y, z, k),
        .y = x[k],
        .z = t_sqr(y, k) - z[k]
    };
}

int main (int argc, char **argv) {
    long order, steps;
    real x0, y0, z0, stepsize;

    assert(argc == 11);
    t_stepper(argv, &order, &stepsize, &steps);
    parameters p;
    t_args(argv, argc, &x0, &y0, &z0, &p.alpha, &p.kappa);

    tsm(order, steps, stepsize, x0, y0, z0, &p, NULL, ode);
    return 0;
}

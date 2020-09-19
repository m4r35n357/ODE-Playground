/*
 * Wimol-Banlue System
 *
 * Example: ./tsm-wimol-banlue-dbg NA NA 8 0.1 10000 1 0 0 2.0
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    series a;
} parameters;

typedef struct {
    series tx;
    series s2x;
} intermediates;

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    intermediates *i = (intermediates *)inters;
    i->tx[k] = t_tan_sec2(i->tx, i->s2x, x, k, HYP).a;
    return (components) {
        .x = y[k] - x[k],
        .y = - t_prod(z, i->tx, k),
        .z = - p->a[k] + t_prod(x, y, k) + t_abs(y, k)
    };
}

int main (int argc, char **argv) {
    long order, steps;
    real stepsize, x0, y0, z0;

    assert(argc == 10);
    t_stepper(argv, &order, &stepsize, &steps);
    parameters p = (parameters) {
        .a = t_jet(order)
    };
    t_args(argv, argc, &x0, &y0, &z0, p.a);
    intermediates i = (intermediates) {
        .tx = t_jet(order),
        .s2x = t_jet(order)
    };

    taylor(order, steps, stepsize, x0, y0, z0, &p, &i, ode);
    return 0;
}

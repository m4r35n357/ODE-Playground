/*
 * Bouali Attractor
 *
 * Example: ./tsm-bouali-dbg NA NA 10 0.01 50000 1 1 0 3 2.2 1 .01
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real a;
    real b;
    real c;
    real d;
} parameters;

typedef struct {
    series wa;
    series wb;
    series w1;
} intermediates;

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    intermediates *i = (intermediates *)inters;
    i->wa[k] = i->w1[k] - y[k];
    i->wb[k] = i->w1[k] - t_sqr(x, k);
    return (components) {
        .x = p->a * t_prod(x, i->wa, k) - p->b * z[k],
        .y = - p->c * t_prod(y, i->wb, k),
        .z = p->d * x[k]
    };
}

int main (int argc, char **argv) {
    long order, steps;
    real x0, y0, z0, stepsize;

    assert(argc == 13);
    t_stepper(argv, &order, &stepsize, &steps);
    parameters p;
    t_args(argv, argc, &x0, &y0, &z0, &p.a, &p.b, &p.c, &p.d);
    intermediates i = (intermediates) {
        .wa = t_jet(order),
        .wb = t_jet(order),
        .w1 = t_jet_c(order, 1.0)
    };

    tsm(order, steps, stepsize, x0, y0, z0, &p, &i, ode);
    return 0;
}

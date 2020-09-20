/*
 * Rabinovich–Fabrikant System
 *
 * Example: ./tsm-rf-dbg NA NA 10 .01 50000 .05 -.05 .3 .28713 .1
 *          ./tsm-rf-dbg NA NA 16 .01 50000 .05 -.05 .3 .105 .1 | ./plotPi3d.py
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    series alpha;
    real gamma;
} parameters;

typedef struct {
    series a;
    series b;
    series c;
    series w1;
    real x2_1;
} intermediates;

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    intermediates *i = (intermediates *)inters;
    i->x2_1 = t_sqr(x, k) - i->w1[k];
    i->a[k] = z[k] + i->x2_1;
    i->b[k] = 3.0 * z[k] - i->x2_1;
    i->c[k] = p->alpha[k] + t_prod(x, y, k);
    return (components) {
        .x = t_prod(y, i->a, k) + p->gamma * x[k],
        .y = t_prod(x, i->b, k) + p->gamma * y[k],
        .z = - 2.0 * t_prod(z, i->c, k)
    };
}

int main (int argc, char **argv) {
    long order, steps;
    real stepsize, x0, y0, z0;

    assert(argc == 11);
    t_stepper(argv, &order, &stepsize, &steps);
    parameters p = (parameters) {
        .alpha = t_jet(order)
    };
    t_args(argv, argc, &x0, &y0, &z0, p.alpha, &p.gamma);
    intermediates i = (intermediates) {
        .a = t_jet(order),
        .b = t_jet(order),
        .c = t_jet(order),
        .w1 = t_jet_c(order, 1.0)
    };

    tsm(order, steps, stepsize, x0, y0, z0, &p, &i, ode);
    return 0;
}

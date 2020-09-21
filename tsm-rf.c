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

static void *get_p (int argc, char **argv, long order) {
    parameters *p = malloc(sizeof (parameters));
    p->alpha = t_jet(order);
    t_args(argv, argc, p->alpha, &p->gamma);
    return p;
}

static void *get_i (long order) {
    intermediates *i = malloc(sizeof (intermediates));
    i->a = t_jet(order);
    i->b = t_jet(order);
    i->c = t_jet(order);
    i->w1 = t_jet_c(order, 1.0);
    return i;
}

int main (int argc, char **argv) {
    assert(argc == 11);
    tsm(argc, argv, ode, get_p, get_i);
    return 0;
}

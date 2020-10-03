/*
 * Rabinovich–Fabrikant System
 *
 * Example: ./tsm-rf-dbg 15 NA 10 .01 50000 .05 -.05 .3 .28713 .1
 *          ./tsm-rf-dbg 15 NA 16 .01 50000 .05 -.05 .3 .105 .1 | ./plotPi3d.py
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

static void *get_p (int argc, char **argv, long order) {
    parameters *p = malloc(sizeof (parameters));
    p->alpha = t_jet(order);
    t_args(argv, argc, p->alpha, &p->gamma);
    return p;
}

typedef struct {
    series a;
    series b;
    series c;
    series w1;
} intermediates;

static void *get_i (long order) {
    intermediates *i = malloc(sizeof (intermediates));
    i->a = t_jet(order);
    i->b = t_jet(order);
    i->c = t_jet(order);
    i->w1 = t_jet_c(order, 1.0);
    return i;
}

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    intermediates *i = (intermediates *)inters;
    i->a[k] = z[k] + t_sqr(x, k) - i->w1[k];
    i->b[k] = 4.0 * z[k] - i->a[k];
    i->c[k] = p->alpha[k] + t_prod(x, y, k);
    return (components) {
        .x = t_prod(y, i->a, k) + p->gamma * x[k],
        .y = t_prod(x, i->b, k) + p->gamma * y[k],
        .z = - 2.0 * t_prod(z, i->c, k)
    };
}

int main (int argc, char **argv) {
    assert(argc == 11);
    tsm(argc, argv, ode, get_p, get_i);
    return 0;
}

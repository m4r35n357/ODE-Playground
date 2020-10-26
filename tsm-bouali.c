/*
 * Bouali Attractor
 *
 * Example: ./tsm-bouali-dbg 15 10 0.01 50000 1 1 0 3 2.2 1 .01
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

static void *get_p (int argc, char **argv, long order) {
    (void)order;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->c, &p->d);
    return p;
}

typedef struct {
    series wa;
    series wb;
    series w1;
} intermediates;

static void *get_i (long order) {
    intermediates *i = malloc(sizeof (intermediates));
    i->wa = t_jet(order);
    i->wb = t_jet(order);
    i->w1 = t_jet_c(order, 1.0);
    return i;
}

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
    assert(argc == 12);
    tsm(argc, argv, ode, get_p, get_i);
    return 0;
}

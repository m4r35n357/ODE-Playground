/*
 * Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-thomas-dbg NA NA 10 0.1 30000 1 0 0 .19
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real b;
} parameters;

typedef struct {
    series sx; series cx;
    series sy; series cy;
    series sz; series cz;
} intermediates;

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    intermediates *i = (intermediates *)inters;
    return (components) {
        .x = t_sin_cos(i->sy, i->cy, y, k, TRIG).a - p->b * x[k],
        .y = t_sin_cos(i->sz, i->cz, z, k, TRIG).a - p->b * y[k],
        .z = t_sin_cos(i->sx, i->cx, x, k, TRIG).a - p->b * z[k]
    };
}

static void *get_p (int argc, char **argv, long order) {
    (void)order;
    parameters *p = malloc(sizeof (parameters));
    t_args(argv, argc, &p->b);
    return p;
}

static void *get_i (long order) {
    intermediates *i = malloc(sizeof (intermediates));
    i->sx = t_jet(order); i->cx = t_jet(order);
    i->sy = t_jet(order); i->cy = t_jet(order);
    i->sz = t_jet(order); i->cz = t_jet(order);
    return i;
}

int main (int argc, char **argv) {
    assert(argc == 10);
    tsm(argc, argv, ode, get_p, get_i);
    return 0;
}

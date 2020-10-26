/*
 * Yu-Wang System
 *
 * Example: ./tsm-yu-wang-dbg 15 10 .001 50000 1 0 0 10 40 2 2.5
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
    series xy;
    series e_xy;
} intermediates;

static void *get_i (long order) {
    intermediates *i = malloc(sizeof (intermediates));
    i->xy = t_jet(order);
    i->e_xy = t_jet(order);
    return i;
}

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    intermediates *i = (intermediates *)inters;
    i->xy[k] = t_prod(x, y, k);
    return (components) {
        .x = p->a * (y[k] - x[k]),
        .y = p->b * x[k] - p->c * t_prod(x, z, k),
        .z = t_exp(i->e_xy, i->xy, k) - p->d * z[k]
    };
}

int main (int argc, char **argv) {
    assert(argc == 12);
    tsm(argc, argv, ode, get_p, get_i);
    return 0;
}

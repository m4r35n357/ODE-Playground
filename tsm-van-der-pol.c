/*
 * Van der Pol oscillator
 *
 * Example:  ./tsm-van-der-pol-dbg 15 10 .01 10000 1 0 0 5
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real mu;
} parameters;

static void *get_p (int argc, char **argv, long order) {
    (void)order;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->mu);
    return p;
}

typedef struct {
    series x2;
} intermediates;

static void *get_i (long order) {
    intermediates *i = malloc(sizeof (intermediates));
    i->x2 = t_jet(order);
    return i;
}

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    intermediates *i = (intermediates *)inters;
    (void)z;
    i->x2[k] = t_prod(x, x, k);
    return (components) {
        .x = y[k],
        .y = p->mu * (y[k] - t_prod(i->x2, y, k)) - x[k]
    };
}

int main (int argc, char **argv) {
    assert(argc == 9);
    tsm(argc, argv, ode, get_p, get_i);
    return 0;
}

/*
 * Wimol-Banlue System
 *
 * Example: ./tsm-wimol-banlue-dbg 15 8 0.1 10000 1 0 0 2.0
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    series a;
} parameters;

static void *get_p (int argc, char **argv, long order) {
    parameters *p = malloc(sizeof (parameters));
    p->a = t_jet(order);
    t_params(argv, argc, p->a);
    return p;
}

typedef struct {
    series tx;
    series s2x;
} intermediates;

static void *get_i (long order) {
    intermediates *i = malloc(sizeof (intermediates));
    i->tx = t_jet(order);
    i->s2x = t_jet(order);
    return i;
}

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
    assert(argc == 9);
    tsm(argc, argv, ode, get_p, get_i);
    return 0;
}

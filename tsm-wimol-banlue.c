/*
 * Wimol-Banlue System
 *
 * Example: ./tsm-wimol-banlue-dbg 15 8 0.1 10000 1 0 0 2.0
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct { real a; series tx, s2x; } parameters;

void *get_p (int argc, char **argv, long order) {
    assert(argc == 9);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a);
    p->tx = t_jet(order);
    p->s2x = t_jet(order);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    t_tan_sec2(p->tx, p->s2x, x, k, HYP);
    return (components) {
        .x = y[k] - x[k],
        .y = - t_prod(z, p->tx, k),
        .z = - (k == 0 ? p->a : 0.0L) + t_prod(x, y, k) + t_abs(y, k)
    };
}

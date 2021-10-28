/*
 * Van der Pol oscillator
 *
 * Example:  ./tsm-van-der-pol-dbg 15 10 .01 10000 1 0 0 5
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct { real mu; series x2; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 9);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->mu);
    p->x2 = t_jet(n);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    (void)z;
    p->x2[k] = t_sqr(x, k);
    return (components) {
        .x = y[k],
        .y = p->mu * (y[k] - t_prod(p->x2, y, k)) - x[k],
        .z = 0.0L
    };
}

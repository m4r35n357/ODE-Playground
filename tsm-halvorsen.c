/*
 * Halvorsen Cyclic Attractor
 *
 * Example: ./tsm-halvorsen-dbg 15 NA 10 .01 10000 1 0 0 1.4
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real a;
} parameters;

static void *get_p (int argc, char **argv, long order) {
    (void)order;
    parameters *p = malloc(sizeof (parameters));
    t_args(argv, argc, &p->a);
    return p;
}

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    (void)inters;
    return (components) {
        .x = - p->a * x[k] - 4.0 * (y[k] + z[k]) - t_sqr(y, k),
        .y = - p->a * y[k] - 4.0 * (z[k] + x[k]) - t_sqr(z, k),
        .z = - p->a * z[k] - 4.0 * (x[k] + y[k]) - t_sqr(x, k)
    };
}

int main (int argc, char **argv) {
    assert(argc == 10);
    tsm(argc, argv, ode, get_p, NULL);
    return 0;
}

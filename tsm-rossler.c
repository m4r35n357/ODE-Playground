/*
 * Rossler System
 *
 * Example: ./tsm-rossler-dbg 15 10 0.01 50000 0.0 -6.78 0.02 .2 .2 5.7
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real a;
    series b;
    real c;
} parameters;

static void *get_p (int argc, char **argv, long order) {
    parameters *p = malloc(sizeof (parameters));
    p->b = t_jet(order);
    t_params(argv, argc, &p->a, p->b, &p->c);
    return p;
}

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    (void)inters;
    return (components) {
        .x = - y[k] - z[k],
        .y = x[k] + p->a * y[k],
        .z = p->b[k] + t_prod(x, z, k) - p->c * z[k]
    };
}

int main (int argc, char **argv) {
    assert(argc == 11);
    tsm(argc, argv, ode, get_p, NULL);
    return 0;
}

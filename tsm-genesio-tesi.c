/*
 * Genesio-Tesi System
 *
 * Example: ./tsm-genesio-tesi-dbg 15 NA 10 0.01 50000 .1 .1 .1 .44 1.1 1
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
} parameters;

static void *get_p (int argc, char **argv, long order) {
    (void)order;
    parameters *p = malloc(sizeof (parameters));
    t_args(argv, argc, &p->a, &p->b, &p->c);
    return p;
}

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    (void)inters;
    return (components) {
        .x = y[k],
        .y = z[k],
        .z = t_sqr(x, k) - p->c * x[k] - p->b * y[k] - p->a * z[k]
    };
}

int main (int argc, char **argv) {
    assert(argc == 12);
    tsm(argc, argv, ode, get_p, NULL);
    return 0;
}

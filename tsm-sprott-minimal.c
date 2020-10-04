/*
 * Sprott Minimal System
 *
 * Example: ./tsm-sprott-minimal-dbg 15 _ 10 0.01 10000 .02 0 0 2.017
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
    t_params(argv, argc, &p->a);
    return p;
}

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    (void)inters;
    return (components) {
        .x = y[k],
        .y = z[k],
        .z = - p->a * z[k] + t_sqr(y, k) - x[k]
    };
}

int main (int argc, char **argv) {
    assert(argc == 10);
    tsm(argc, argv, ode, get_p, NULL);
    return 0;
}

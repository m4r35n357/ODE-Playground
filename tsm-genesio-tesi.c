/*
 * Genesio-Tesi System
 *
 * Example: ./tsm-genesio-tesi-dbg 15 10 0.01 50000 .1 .1 .1 .44 1.1 1
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

void *get_p (int argc, char **argv, long order) {
    assert(argc == 11);
    (void)order;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->c);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = y[k],
        .y = z[k],
        .z = t_sqr(x, k) - p->c * x[k] - p->b * y[k] - p->a * z[k]
    };
}

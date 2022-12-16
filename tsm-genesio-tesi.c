/*
 * Genesio-Tesi System - http://www.atomosyd.net/spip.php?article153
 *
 * Example: ./tsm-genesio-tesi-std 15 10 0.01 50000 .1 .1 .1 .44 1.1
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct Parameters { real a, b; } parameters;

void *get_p (int argc, char **argv, int n) { (void)n;
    assert(argc == 10);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = y[k],
        .y = z[k],
        .z = - t_sqr(x, k) - x[k] - p->b * y[k] - p->a * z[k]
    };
}

/*
 * Burke & Shaw System - http://www.atomosyd.net/spip.php?article33
 *
 * Example: ./tsm-burke-shaw-std 9 10 .01 10000 1 1 1 10 4.272
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct Parameters { real s, v; } parameters;

void *get_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 10);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->s, &p->v);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = - p->s * (x[k] + y[k]),
        .y = - (p->s * t_mul(x, z, k) + y[k]),
        .z = p->s * t_mul(x, y, k) + p->v
    };
}

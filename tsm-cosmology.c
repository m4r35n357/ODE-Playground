/*
 * No-frills cosmology (https://www.youtube.com/watch?v=vcGDqRm7ZK4&list=PLaNkJORnlhZkgIyPFNxhJPIVewGckJCGr&index=7)
 *
 * Example:  ./tsm-cosmology-dbg 15 10 .01 10000  1.0 1.0 0.0  0.0
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct { real w; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 11);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->w);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    (void)z;
    return (components) {  // .x maps to rho, .y to theta
        .x = - (1.0L + p->w) * t_mul(x, y, k),
        .y = - t_sqr(y, k) / 3.0L - 4.0L * MY_PI * (1.0L + 3.0L * p->w) * x[k]
    };
}

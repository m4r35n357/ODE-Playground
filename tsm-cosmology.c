/*
 * No-frills cosmology (https://www.youtube.com/watch?v=vcGDqRm7ZK4&list=PLaNkJORnlhZkgIyPFNxhJPIVewGckJCGr&index=7)
 *
 * Example:  ./tsm-cosmology-dbg 15 10 .01 10000  1.0 1.0 0.0  0.0
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

static const real MY_PI = 3.1415926535897932384626433832795029L;

typedef struct {
    real w;
} parameters;

static void *get_p (int argc, char **argv, long order) {
    (void)order;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->w);
    return p;
}

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    (void)inters;
    (void)z;
    return (components) {
        .x = - (1.0L + p->w) * t_prod(x, y, k),
        .y = - t_prod(y, y, k) / 3.0L - 4.0L * MY_PI * (1.0 + 3.0 * p->w) * x[k]
    };
}

int main (int argc, char **argv) {
    assert(argc == 9);
    tsm(argc, argv, ode, get_p, NULL);
    return 0;
}

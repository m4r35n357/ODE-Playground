/*
 * Lotka-Volterra (Predator-Prey) System
 *
 * Example: ./tsm-lotka-volterra-dbg 15 10 .01 10000 10 10 0 1 .5 .05 .02
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
    real d;
} parameters;

static void *get_p (int argc, char **argv, long order) {
    (void)order;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->c, &p->d);
    return p;
}

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    (void)inters;
    (void)z;
    real xy = t_prod(x, y, k);
    return (components) {
        .x = p->a * x[k] - p->c * xy,
        .y = p->d * xy - p->b * y[k]
    };
}

int main (int argc, char **argv) {
    assert(argc == 12);
    tsm(argc, argv, ode, get_p, NULL);
    return 0;
}

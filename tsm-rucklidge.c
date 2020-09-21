/*
 * Rucklidge Attractor
 *
 * Example: ./tsm-rucklidge-dbg NA NA 10 0.01 15000 1 0 0 6.7 2
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real alpha;
    real kappa;
} parameters;

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    (void)inters;
    return (components) {
        .x = p->alpha * y[k] - p->kappa * x[k] - t_prod(y, z, k),
        .y = x[k],
        .z = t_sqr(y, k) - z[k]
    };
}

static void *get_p (int argc, char **argv, long order) {
    (void)order;
    parameters *p = malloc(sizeof (parameters));
    t_args(argv, argc, &p->alpha, &p->kappa);
    return p;
}

int main (int argc, char **argv) {
    assert(argc == 11);
    tsm(argc, argv, ode, get_p, NULL);
    return 0;
}

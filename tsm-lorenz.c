/*
 * Lorenz System
 *
 * Example: ./tsm-lorenz-dbg 15 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real sigma;
    real rho;
    real beta;
} parameters;

static void *get_p (int argc, char **argv, long order) {
    (void)order;
    parameters *p = malloc(sizeof (parameters));
    real _;
    t_params(argv, argc, &p->sigma, &p->rho, &p->beta, &_);
    p->beta /= _;
    return p;
}

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    (void)inters;
    return (components) {
        .x = p->sigma * (y[k] - x[k]),
        .y = p->rho * x[k] - y[k] - t_prod(x, z, k),
        .z = t_prod(x, y, k) - p->beta * z[k]
    };
}

int main (int argc, char **argv) {
    assert(argc == 12);
    tsm(argc, argv, ode, get_p, NULL);
    return 0;
}

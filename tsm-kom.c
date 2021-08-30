/*
 * Kom, Kengne et. al.
 *
 * Example: ./tsm-kom-dbg 15 10 0.01 10000 1 0 0 1 1 1.1 5.4291e-4
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real alpha;
    real rho;
    real mu;
    series epsilon;
    series e_y;
} parameters;

void *get_p (int argc, char **argv, long order) {
    assert(argc == 12);
    parameters *p = malloc(sizeof (parameters));
    p->epsilon = t_jet(order);
    t_params(argv, argc, &p->alpha, &p->rho, &p->mu, p->epsilon);
    p->e_y = t_jet(order);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    p->e_y[k] = t_exp(p->e_y, y, k);
    return (components) {
        .x = y[k],
        .y = p->alpha * z[k],
        .z = - p->rho * x[k] - p->mu * z[k] - p->epsilon[0] * p->e_y[k] + p->epsilon[k]
    };
}

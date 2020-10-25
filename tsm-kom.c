/*
 * Kom, Kengne et. al.
 *
 * Example: ./tsm-kom-dbg 15 _ 10 0.01 10000 1 0 0 1 1 1.1 5.4291e-4
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real alpha;
    real rho;
    real mu;
    series epsilon;
} parameters;

static void *get_p (int argc, char **argv, long order) {
    parameters *p = malloc(sizeof (parameters));
    p->epsilon = t_jet(order);
    t_params(argv, argc, &p->alpha, &p->rho, &p->mu, p->epsilon);
    return p;
}

typedef struct {
    series e_y;
} intermediates;

static void *get_i (long order) {
    intermediates *i = malloc(sizeof (intermediates));
    i->e_y = t_jet(order);
    return i;
}

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    intermediates *i = (intermediates *)inters;
    i->e_y[k] = t_exp(i->e_y, y, k);
    return (components) {
        .x = y[k],
        .y = p->alpha * z[k],
        .z = - p->rho * x[k] - p->mu * z[k] - p->epsilon[0] * i->e_y[k] + p->epsilon[k]
    };
}

int main (int argc, char **argv) {
    assert(argc == 13);
    tsm(argc, argv, ode, get_p, get_i);
    return 0;
}


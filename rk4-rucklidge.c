/*
 * Rucklidge Attractor
 *
 * Example: ./rk4-rucklidge-dbg 15 NA 1 0.01 15000 1 0 0 6.7 2
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

static void *get_p (int argc, char **argv) {
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->alpha, &p->kappa);
    return p;
}

static components ode (real x, real y, real z, void *params) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = p->alpha * y - p->kappa * x - y * z,
        .y = x,
        .z = y * y - z
    };
}

int main (int argc, char **argv) {
    assert(argc == 11);
    rk4(argc, argv, ode, get_p);
    return 0;
}

/*
 * Sprott-Jafari System
 *
 * Example: ./tsm-sj-dbg 15 _ 10 .01 10000 0 3.9 .7 8.888 4
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real a;
    series b;
} parameters;

static void *get_p (int argc, char **argv, long order) {
    parameters *p = malloc(sizeof (parameters));
    p->b = t_jet(order);
    t_params(argv, argc, &p->a, p->b);
    return p;
}

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    (void)inters;
    return (components) {
        .x = y[k],
        .y = - x[k] + t_prod(y, z, k),
        .z = z[k] + p->a * t_sqr(x, k) - t_sqr(y, k) - p->b[k]
    };
}

int main (int argc, char **argv) {
    assert(argc == 11);
    tsm(argc, argv, ode, get_p, NULL);
    return 0;
}

/*
 * Bouali Attractor
 *
 * Example: ./tsm-bouali-dbg 15 10 0.01 50000 1 1 0 3 2.2 1 .01
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real a;
    real b;
    real c;
    real d;
    series wa;
    series wb;
} parameters;

void *get_p (int argc, char **argv, long order) {
    assert(argc == 12);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->c, &p->d);
    p->wa = t_jet(order);
    p->wb = t_jet(order);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    p->wa[k] = (k == 0 ? 1.0L : 0.0L) - y[k];
    p->wb[k] = (k == 0 ? 1.0L : 0.0L) - t_prod(x, x, k);
    return (components) {
        .x = p->a * t_prod(x, p->wa, k) - p->b * z[k],
        .y = - p->c * t_prod(y, p->wb, k),
        .z = p->d * x[k]
    };
}

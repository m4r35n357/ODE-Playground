/*
 * Bouali Attractor
 *
 * Example: ./tsm-bouali-std 15 10 0.01 50000 1 1 0 3 2.2 1 .01
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

typedef struct Parameters { real a, b, c, d; series wa, wb, w1; } parameters;

void *get_p (int argc, char **argv, int n) {
    CHECK(argc == 12);
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
    t_params(argv, argc, &p->a, &p->b, &p->c, &p->d);
    p->wa = t_jet(n);
    p->wb = t_jet(n);
    p->w1 = t_const(n, 1.0L);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    p->wa[k] = p->w1[k] - y[k];
    p->wb[k] = p->w1[k] - t_sqr(x, k);
    return (components) {
        .x = p->a * t_mul(x, p->wa, k) - p->b * z[k],
        .y = - p->c * t_mul(y, p->wb, k),
        .z = p->d * x[k]
    };
}

/*
 * Wimol-Banlue System
 *
 * Example: ./tsm-wimol-banlue-std 15 8 0.1 10000 1 0 0 2.0
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

typedef struct Parameters { real a; series tx, s2x, wa; } parameters;

void *get_p (int argc, char **argv, int n) {
    CHECK(argc == 9);
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
    t_params(argv, argc, &p->a);
    p->tx = t_jet(n);
    p->s2x = t_jet(n);
    p->wa = t_const(n, p->a);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    t_tan_sec2(p->tx, p->s2x, x, k, false);
    return (components) {
        .x = y[k] - x[k],
        .y = - t_mul(z, p->tx, k),
        .z = - p->wa[k] + t_mul(x, y, k) + t_abs(y, k)
    };
}

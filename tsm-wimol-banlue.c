/*
 * Wimol-Banlue System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

typedef struct Parameters { real a; series tx, s2x, _A; } parameters;

void *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 9);
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
    tsm_get_p(argv, argc, &p->a);
    p->tx = t_jet(n);
    p->s2x = t_jet(n);
    p->_A = t_const(n, p->a);
    return p;
}

triplet ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    t_tan_sec2(p->tx, p->s2x, x, k, false);
    return (triplet) {
        .x = y[k] - x[k],
        .y = - t_mul(z, p->tx, k),
        .z = - p->_A[k] + t_mul(x, y, k) + t_abs(y, k)
    };
}

/*
 * Rabinovich–Fabrikant System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

typedef struct Parameters { real alpha, gamma; series sa, sb, sc, _ALPHA, _1; } parameters;

void *get_p (int argc, char **argv, int n) {
    CHECK(argc == 10);
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
    t_params(argv, argc, &p->alpha, &p->gamma);
    p->sa = t_jet(n);
    p->sb = t_jet(n);
    p->sc = t_jet(n);
    p->_ALPHA = t_const(n, p->alpha);
    p->_1 = t_const(n, 1.0L);
    return p;
}

triplet ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    p->sa[k] = z[k] + t_sqr(x, k) - p->_1[k];
    p->sb[k] = 4.0L * z[k] - p->sa[k];
    p->sc[k] = p->_ALPHA[k] + t_mul(x, y, k);
    return (triplet) {
        .x = t_mul(y, p->sa, k) + p->gamma * x[k],
        .y = t_mul(x, p->sb, k) + p->gamma * y[k],
        .z = - 2.0L * t_mul(z, p->sc, k)
    };
}

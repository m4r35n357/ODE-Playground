/*
 * Bouali Attractor
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

typedef struct Parameters { real a, b, c, d; series sa, sb, _1; } parameters;

void *get_p (int argc, char **argv, int n) {
    CHECK(argc == 12);
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
    t_params(argv, argc, &p->a, &p->b, &p->c, &p->d);
    p->sa = t_jet(n);
    p->sb = t_jet(n);
    p->_1 = t_const(n, 1.0L);
    return p;
}

triplet ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    p->sa[k] = p->_1[k] - y[k];
    p->sb[k] = p->_1[k] - t_sqr(x, k);
    return (triplet) {
        .x = p->a * t_mul(x, p->sa, k) - p->b * z[k],
        .y = - p->c * t_mul(y, p->sb, k),
        .z = p->d * x[k]
    };
}

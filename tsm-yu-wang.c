/*
 * Yu-Wang System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

typedef struct Parameters { real a, b, c, d; series xy, e_xy; } parameters;

void *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 12);
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
    tsm_get_p(argv, argc, &p->a, &p->b, &p->c, &p->d);
    p->xy = t_jet(n);
    p->e_xy = t_jet(n);
    return p;
}

triplet ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    p->xy[k] = t_mul(x, y, k);
    return (triplet) {
        .x = p->a * (y[k] - x[k]),
        .y = p->b * x[k] - p->c * t_mul(x, z, k),
        .z = t_exp(p->e_xy, p->xy, k) - p->d * z[k]
    };
}

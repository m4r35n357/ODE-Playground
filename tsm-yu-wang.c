/*
 * Yu-Wang System
 *
 * Example: ./tsm-yu-wang-std 15 10 .001 50000 1 0 0 10 40 2 2.5
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

typedef struct Parameters { real a, b, c, d; series xy, e_xy; } parameters;

void *get_p (int argc, char **argv, int n) {
    CHECK(argc == 12);
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
    t_params(argv, argc, &p->a, &p->b, &p->c, &p->d);
    p->xy = t_jet(n);
    p->e_xy = t_jet(n);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    p->xy[k] = t_mul(x, y, k);
    return (components) {
        .x = p->a * (y[k] - x[k]),
        .y = p->b * x[k] - p->c * t_mul(x, z, k),
        .z = t_exp(p->e_xy, p->xy, k) - p->d * z[k]
    };
}

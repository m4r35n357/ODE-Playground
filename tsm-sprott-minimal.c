/*
 * Sprott Minimal System
 *
 * Example: ./tsm-sprott-minimal-std 15 10 0.01 10000 .02 0 0 2.017
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a; };

void *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 9);
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
    tsm_get_p(argv, argc, &p->a);
    return p;
}

triplet ode (series x, series y, series z, parameters *p, int k) {
    return (triplet) {
        .x = y[k],
        .y = z[k],
        .z = - p->a * z[k] + t_sqr(y, k) - x[k]
    };
}

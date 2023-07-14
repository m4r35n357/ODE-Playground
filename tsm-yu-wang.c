/*
 * Yu-Wang System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a, b, c, d; series xy, e_xy; };

parameters *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 12);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    tsm_get_p(argv, argc, &_->a, &_->b, &_->c, &_->d);
    _->xy = t_jet(n);
    _->e_xy = t_jet(n);
    return _;
}

triplet ode (series x, series y, series z, parameters *p, int k) {
    p->xy[k] = t_mul(x, y, k);
    return (triplet) {
        .x = p->a * (y[k] - x[k]),
        .y = p->b * x[k] - p->c * t_mul(x, z, k),
        .z = t_exp(p->e_xy, p->xy, k) - p->d * z[k]
    };
}

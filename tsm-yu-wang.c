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
    _->xy = tsm_jet(n);
    _->e_xy = tsm_jet(n);
    return _;
}

triplet ode (series x, series y, series z, parameters *_, int k) {
    _->xy[k] = t_mul(x, y, k);
    return (triplet) {
        .x = _->a * (y[k] - x[k]),
        .y = _->b * x[k] - _->c * t_mul(x, z, k),
        .z = t_exp(_->e_xy, _->xy, k) - _->d * z[k]
    };
}

/*
 * Bouali Attractor
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a, b, c, d; series sa, sb, _1; };

parameters *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 12);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    tsm_get_p(argv, argc, &_->a, &_->b, &_->c, &_->d);
    _->sa = t_jet(n);
    _->sb = t_jet(n);
    _->_1 = t_const(n, 1.0L);
    return _;
}

triplet ode (series x, series y, series z, parameters *_, int k) {
    _->sa[k] = _->_1[k] - y[k];
    _->sb[k] = _->_1[k] - t_sqr(x, k);
    return (triplet) {
        .x = _->a * t_mul(x, _->sa, k) - _->b * z[k],
        .y = - _->c * t_mul(y, _->sb, k),
        .z = _->d * x[k]
    };
}

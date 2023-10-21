/*
 * Bouali Attractor
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a, b, c, d; series sa, sb, _1; };

model *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 12);
    model *_ = malloc(sizeof (model)); CHECK(_);
    _->sa = tsm_jet(n);
    _->sb = tsm_jet(n);
    _->_1 = tsm_jet(n); _->_1[0] = 1.0L;
    tsm_get_p(argv, argc, &_->a, &_->b, &_->c, &_->d);
    return _;
}

triplet ode (series x, series y, series z, model *_, int k) {
    _->sa[k] = _->_1[k] - y[k];
    _->sb[k] = _->_1[k] - t_sqr(x, k);
    return (triplet) {
        .x = _->a * t_mul(x, _->sa, k) - _->b * z[k],
        .y = - _->c * t_mul(y, _->sb, k),
        .z = _->d * x[k]
    };
}

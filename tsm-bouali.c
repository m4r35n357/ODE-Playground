/*
 * Bouali Attractor
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a, b, c, d; series sa, sb, _1; };

void *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 12);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    tsm_get_p(argv, argc, &_->a, &_->b, &_->c, &_->d);
    _->sa = t_jet(n);
    _->sb = t_jet(n);
    _->_1 = t_const(n, 1.0L);
    return _;
}

triplet ode (series x, series y, series z, parameters *p, int k) {
    p->sa[k] = p->_1[k] - y[k];
    p->sb[k] = p->_1[k] - t_sqr(x, k);
    return (triplet) {
        .x = p->a * t_mul(x, p->sa, k) - p->b * z[k],
        .y = - p->c * t_mul(y, p->sb, k),
        .z = p->d * x[k]
    };
}

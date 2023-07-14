/*
 * Rabinovich–Fabrikant System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real alpha, gamma; series sa, sb, sc, _ALPHA, _1; };

void *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 10);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    tsm_get_p(argv, argc, &_->alpha, &_->gamma);
    _->sa = t_jet(n);
    _->sb = t_jet(n);
    _->sc = t_jet(n);
    _->_ALPHA = t_const(n, _->alpha);
    _->_1 = t_const(n, 1.0L);
    return _;
}

triplet ode (series x, series y, series z, parameters *p, int k) {
    p->sa[k] = z[k] + t_sqr(x, k) - p->_1[k];
    p->sb[k] = 4.0L * z[k] - p->sa[k];
    p->sc[k] = p->_ALPHA[k] + t_mul(x, y, k);
    return (triplet) {
        .x = t_mul(y, p->sa, k) + p->gamma * x[k],
        .y = t_mul(x, p->sb, k) + p->gamma * y[k],
        .z = - 2.0L * t_mul(z, p->sc, k)
    };
}

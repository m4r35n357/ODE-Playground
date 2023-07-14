/*
 * Sprott-Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-sprott-thomas-std 15 10 0.01 30000 1 0 0 4.75 1
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a, b; series ax, ay, az, sax, say, saz, cax, cay, caz, tx, ty, tz, s2x, s2y, s2z; };

parameters *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 10);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    tsm_get_p(argv, argc, &_->a, &_->b);
    _->ax = t_jet(n); _->sax = t_jet(n); _->cax = t_jet(n); _->tx = t_jet(n); _->s2x = t_jet(n);
    _->ay = t_jet(n); _->say = t_jet(n); _->cay = t_jet(n); _->ty = t_jet(n); _->s2y = t_jet(n);
    _->az = t_jet(n); _->saz = t_jet(n); _->caz = t_jet(n); _->tz = t_jet(n); _->s2z = t_jet(n);
    return _;
}

triplet ode (series x, series y, series z, parameters *p, int k) {
    p->ax[k] = p->a * x[k];
    p->ay[k] = p->a * y[k];
    p->az[k] = p->a * z[k];
    return (triplet) {
        .x = t_sin_cos(p->say, p->cay, p->ay, k, true).a - p->b * t_tan_sec2(p->tx, p->s2x, x, k, true).a,
        .y = t_sin_cos(p->saz, p->caz, p->az, k, true).a - p->b * t_tan_sec2(p->ty, p->s2y, y, k, true).a,
        .z = t_sin_cos(p->sax, p->cax, p->ax, k, true).a - p->b * t_tan_sec2(p->tz, p->s2z, z, k, true).a
    };
}

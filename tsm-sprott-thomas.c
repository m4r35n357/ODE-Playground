/*
 * Sprott-Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-sprott-thomas-std 15 10 0.01 30000 1 0 0 4.75 1
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct Parameters { real a, b; series ax, ay, az, sax, say, saz, cax, cay, caz, tx, ty, tz, s2x, s2y, s2z; } parameters;

void *get_p (int argc, char **argv, int n) {
    CHECK(argc == 10);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b);
    p->ax = t_jet(n); p->sax = t_jet(n); p->cax = t_jet(n); p->tx = t_jet(n); p->s2x = t_jet(n);
    p->ay = t_jet(n); p->say = t_jet(n); p->cay = t_jet(n); p->ty = t_jet(n); p->s2y = t_jet(n);
    p->az = t_jet(n); p->saz = t_jet(n); p->caz = t_jet(n); p->tz = t_jet(n); p->s2z = t_jet(n);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    p->ax[k] = p->a * x[k];
    p->ay[k] = p->a * y[k];
    p->az[k] = p->a * z[k];
    return (components) {
        .x = t_sin_cos(p->say, p->cay, p->ay, k, true).a - p->b * t_tan_sec2(p->tx, p->s2x, x, k, true).a,
        .y = t_sin_cos(p->saz, p->caz, p->az, k, true).a - p->b * t_tan_sec2(p->ty, p->s2y, y, k, true).a,
        .z = t_sin_cos(p->sax, p->cax, p->ax, k, true).a - p->b * t_tan_sec2(p->tz, p->s2z, z, k, true).a
    };
}

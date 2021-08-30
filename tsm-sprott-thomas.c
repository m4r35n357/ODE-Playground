/*
 * Sprott-Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-sprott-thomas-dbg 15 10 0.01 30000 1 0 0 4.75 1
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real a;
    real b;
    series ax; series sax; series cax; series tx; series s2x;
    series ay; series say; series cay; series ty; series s2y;
    series az; series saz; series caz; series tz; series s2z;
} parameters;

void *get_p (int argc, char **argv, long order) {
    assert(argc == 10);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b);
    p->ax = t_jet(order); p->sax = t_jet(order); p->cax = t_jet(order); p->tx = t_jet(order); p->s2x = t_jet(order);
    p->ay = t_jet(order); p->say = t_jet(order); p->cay = t_jet(order); p->ty = t_jet(order); p->s2y = t_jet(order);
    p->az = t_jet(order); p->saz = t_jet(order); p->caz = t_jet(order); p->tz = t_jet(order); p->s2z = t_jet(order);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    p->ax[k] = p->a * x[k];
    p->ay[k] = p->a * y[k];
    p->az[k] = p->a * z[k];
    return (components) {
        .x = t_sin_cos(p->say, p->cay, p->ay, k, TRIG).a - p->b * t_tan_sec2(p->tx, p->s2x, x, k, TRIG).a,
        .y = t_sin_cos(p->saz, p->caz, p->az, k, TRIG).a - p->b * t_tan_sec2(p->ty, p->s2y, y, k, TRIG).a,
        .z = t_sin_cos(p->sax, p->cax, p->ax, k, TRIG).a - p->b * t_tan_sec2(p->tz, p->s2z, z, k, TRIG).a
    };
}

/*
 * Sprott-Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-sprott-thomas-dbg 15 NA 10 0.01 30000 1 0 0 4.75 1
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real a;
    real b;
} parameters;

static void *get_p (int argc, char **argv, long order) {
    (void)order;
    parameters *p = malloc(sizeof (parameters));
    t_args(argv, argc, &p->a, &p->b);
    return p;
}

typedef struct {
    series ax; series sax; series cax; series tx; series s2x;
    series ay; series say; series cay; series ty; series s2y;
    series az; series saz; series caz; series tz; series s2z;
} intermediates;

static void *get_i (long order) {
    intermediates *i = malloc(sizeof (intermediates));
    i->ax = t_jet(order); i->sax = t_jet(order); i->cax = t_jet(order); i->tx = t_jet(order); i->s2x = t_jet(order);
    i->ay = t_jet(order); i->say = t_jet(order); i->cay = t_jet(order); i->ty = t_jet(order); i->s2y = t_jet(order);
    i->az = t_jet(order); i->saz = t_jet(order); i->caz = t_jet(order); i->tz = t_jet(order); i->s2z = t_jet(order);
    return i;
}

static components ode (series x, series y, series z, void *params, void *inters, int k) {
    parameters *p = (parameters *)params;
    intermediates *i = (intermediates *)inters;
    i->ax[k] = p->a * x[k];
    i->ay[k] = p->a * y[k];
    i->az[k] = p->a * z[k];
    return (components) {
        .x = t_sin_cos(i->say, i->cay, i->ay, k, TRIG).a - p->b * t_tan_sec2(i->tx, i->s2x, x, k, TRIG).a,
        .y = t_sin_cos(i->saz, i->caz, i->az, k, TRIG).a - p->b * t_tan_sec2(i->ty, i->s2y, y, k, TRIG).a,
        .z = t_sin_cos(i->sax, i->cax, i->ax, k, TRIG).a - p->b * t_tan_sec2(i->tz, i->s2z, z, k, TRIG).a
    };
}

int main (int argc, char **argv) {
    assert(argc == 11);
    tsm(argc, argv, ode, get_p, get_i);
    return 0;
}

/*
 * Sprott-Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-sprott-thomas-dbg NA NA 10 0.01 30000 1 0 0 4.75 1
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

typedef struct {
    series ax;
    series ay;
    series az;
    series sax;
    series say;
    series saz;
    series cax;
    series cay;
    series caz;
    series tx;
    series ty;
    series tz;
    series s2x;
    series s2y;
    series s2z;
} intermediates;

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
    long order, steps;
    real x0, y0, z0, stepsize;

    assert(argc == 11);
    t_stepper(argv, &order, &stepsize, &steps);
    parameters p;
    t_args(argv, argc, &x0, &y0, &z0, &p.a, &p.b);
    intermediates i = (intermediates) {
        .ax = t_jet(order),
        .ay = t_jet(order),
        .az = t_jet(order),
        .sax = t_jet(order),
        .say = t_jet(order),
        .saz = t_jet(order),
        .cax = t_jet(order),
        .cay = t_jet(order),
        .caz = t_jet(order),
        .tx = t_jet(order),
        .ty = t_jet(order),
        .tz = t_jet(order),
        .s2x = t_jet(order),
        .s2y = t_jet(order),
        .s2z = t_jet(order)
    };

    taylor(order, steps, stepsize, x0, y0, z0, &p, &i, ode);
    return 0;
}

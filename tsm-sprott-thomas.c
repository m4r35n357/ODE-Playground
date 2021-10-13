/*
 * Sprott-Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-sprott-thomas-dbg 9 32 10 0.01 30000 1 0 0 4.75 1
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t a, b; series sax, say, saz, cax, cay, caz, ax, ay, az, tx, ty, tz, sx, sy, sz; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 11);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b);
    p->sax = t_jet(n); p->cax = t_jet(n); p->ax = t_jet(n); p->tx = t_jet(n); p->sx = t_jet(n);
    p->say = t_jet(n); p->cay = t_jet(n); p->ay = t_jet(n); p->ty = t_jet(n); p->sy = t_jet(n);
    p->saz = t_jet(n); p->caz = t_jet(n); p->az = t_jet(n); p->tz = t_jet(n); p->sz = t_jet(n);
    return p;
}

void ode (components *v, series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    mpfr_mul(p->ax[k], x[k], p->a, RND);
    mpfr_mul(p->ay[k], y[k], p->a, RND);
    mpfr_mul(p->az[k], z[k], p->a, RND);
    //  x' = sin(Ay) - Btan(x)
    mpfr_fms(v->x, p->b, *t_tan_sec2(p->tx, p->sx, x, k, TRIG).a, *t_sin_cos(p->say, p->cay, p->ay, k, TRIG).a, RND);
    mpfr_neg(v->x, v->x, RND);
    //  y' = sin(Az) - Btan(y)
    mpfr_fms(v->y, p->b, *t_tan_sec2(p->ty, p->sy, y, k, TRIG).a, *t_sin_cos(p->saz, p->caz, p->az, k, TRIG).a, RND);
    mpfr_neg(v->y, v->y, RND);
    //  z' = sin(Ax) - Btan(z)
    mpfr_fms(v->z, p->b, *t_tan_sec2(p->tz, p->sz, z, k, TRIG).a, *t_sin_cos(p->sax, p->cax, p->ax, k, TRIG).a, RND);
    mpfr_neg(v->z, v->z, RND);
}

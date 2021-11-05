/*
 * Rabinovich–Fabrikant System
 *
 * Example: ./tsm-rf-dbg 9 32 10 .01 50000 .05 -.05 .3 .28713 .1
 *          ./tsm-rf-dbg 9 32 16 .01 50000 .05 -.05 .3 .105 .1 | ./plot3d.py
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t alpha, gamma, D1, D4; series a, b, c; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 11);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->alpha, &p->gamma);
    mpfr_init_set_si(p->D1, 1, RND);
    mpfr_init_set_si(p->D4, 4, RND);
    p->a = t_jet(n); p->b = t_jet(n); p->c = t_jet(n);
    return p;
}

void ode (components *v, series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = y(z - 1 + x^2) + Gx
    mpfr_sub(p->a[k], *t_sqr(x, k), *t_const(p->D1, k), RND);
    mpfr_add(p->a[k], z[k], p->a[k], RND);
    mpfr_fma(v->x, p->gamma, x[k], *t_prod(y, p->a, k), RND);
    //  y' = x(3z + 1 - x^2) + Gy
    mpfr_fms(p->b[k], p->D4, z[k], p->a[k], RND);
    mpfr_fma(v->y, p->gamma, y[k], *t_prod(x, p->b, k), RND);
    //  z' = -2z(A + xy)
    mpfr_add(p->c[k], *t_prod(x, y, k), *t_const(p->alpha, k), RND);
    mpfr_mul_2si(v->z, *t_prod(z, p->c, k), 1, RND);
    mpfr_neg(v->z, v->z, RND);
}

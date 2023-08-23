/*
 * Rabinovich–Fabrikant System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { mpfr_t alpha, gamma, D1, D4; series a, b, c, _ALPHA, _1; };

parameters *get_p (int argc, char **argv, int n) {
    CHECK(argc == 11);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->alpha, &p->gamma);
    mpfr_init_set_si(p->D1, 1, RND);
    mpfr_init_set_si(p->D4, 4, RND);
    p->a = t_jet(n);
    p->b = t_jet(n);
    p->c = t_jet(n);
    p->_ALPHA = t_const(n, p->alpha);
    p->_1 = t_const(n, p->D1);
    return p;
}

void ode (triplet *v_k, series x, series y, series z, parameters *p, int k) {
    //  x' = y(z - 1 + x^2) + Gx
    mpfr_sub(p->a[k], *t_sqr(x, k), p->_1[k], RND);
    mpfr_add(p->a[k], z[k], p->a[k], RND);
    mpfr_fma(v_k->x, p->gamma, x[k], *t_mul(y, p->a, k), RND);
    //  y' = x(3z + 1 - x^2) + Gy
    mpfr_fms(p->b[k], p->D4, z[k], p->a[k], RND);
    mpfr_fma(v_k->y, p->gamma, y[k], *t_mul(x, p->b, k), RND);
    //  z' = -2z(A + xy)
    mpfr_add(p->c[k], *t_mul(x, y, k), p->_ALPHA[k], RND);
    mpfr_mul_2si(v_k->z, *t_mul(z, p->c, k), 1, RND);
    mpfr_neg(v_k->z, v_k->z, RND);
}

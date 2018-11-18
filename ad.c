/*
 * (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"
#include "ad.h"

void set_ad_status (mpfr_t *jet, ad_status s) {
    if (s == VARIABLE) {
        mpfr_set_ui(jet[1], 1, RND);
    } else {
        mpfr_set_ui(jet[1], 0, RND);
    }
}

void jet_output (mpfr_t *jet, long n, char* f_colour, char *fk_colour) {
    mpfr_printf("%s%9.6RNf ", f_colour, jet[0]);
    for (int i = 1; i < n; i++) {
        mpfr_printf("%s%9.6RNf ", fk_colour, jet[i]);
    }
    printf("%s\n", f_colour);
}

void jet_to_derivs (mpfr_t *jet, long n) {
    for (int i = 1, fac = 1; i < n; i++) {
        fac *= i;
        mpfr_mul_ui(jet[i], jet[i], fac, RND);
    }
}

void derivative_output (mpfr_t *jet, long n, char* f_colour, char *fk_colour) {
    jet_to_derivs(jet, n);
    jet_output(jet, n, f_colour, fk_colour);
}

void ad_bisect (model m, mpfr_t *xa, mpfr_t *xb, int max_it, mpfr_t f_tol, mpfr_t x_tol, mpfr_t *xc, mpfr_t *fa, mpfr_t *fc) {
    mpfr_t delta, _;
    int counter = 0;
    mpfr_inits(delta, _, NULL);
    m(fa, xa, 1);
    mpfr_sub(delta, xb[0], xa[0], RND);
    while(mpfr_cmp_abs(fc[0], f_tol) >= 0 || mpfr_cmp_abs(delta, x_tol) >= 0) {
        mpfr_add(xc[0], xa[0], xb[0], RND);
        mpfr_div_ui(xc[0], xc[0], 2, RND);
        m(fc, xc, 1);
        mpfr_mul(_, fa[0], fc[0], RND);
        if (mpfr_sgn(_) < 0) {
            mpfr_set(xb[0], xc[0], RND);
        } else {
            mpfr_set(xa[0], xc[0], RND);
        }
        mpfr_sub(delta, xb[0], xa[0], RND);
        if (++counter > max_it) break;
    }
    mpfr_fprintf(stderr, "%3d %20.12RNe %20.12RNe %20.12RNe\n", counter, xc[0], fc[0], delta);
    mpfr_clears(delta, _, NULL);
}

void ad_newton (model m, mpfr_t *f, mpfr_t *x, int max_it, mpfr_t f_tol, mpfr_t x_tol) {
    mpfr_t delta;
    int counter = 0;
    mpfr_init_set_ui(delta, 1, RND);
    while(mpfr_cmp_abs(f[0], f_tol) >= 0 || mpfr_cmp_abs(delta, x_tol) >= 0) {
        m(f, x, 2);
        mpfr_div(delta, f[0], f[1], RND);
        mpfr_sub(x[0], x[0], delta, RND);
        if (++counter > max_it) break;
    }
    mpfr_fprintf(stderr, "%3d %20.12RNe %20.12RNe %20.12RNe\n", counter, x[0], f[0], delta);
    mpfr_clear(delta);
}

void ad_householder (model m, mpfr_t *f, mpfr_t *x, long n, int max_it, mpfr_t f_tol, mpfr_t x_tol, mpfr_t *f_reciprocal, mpfr_t *w1) {
    mpfr_t delta;
    int counter = 0;
    mpfr_init_set_ui(delta, 1, RND);
    while(mpfr_cmp_abs(f[0], f_tol) >= 0 || mpfr_cmp_abs(delta, x_tol) >= 0) {
        m(f, x, n);
        ad_quotient(f_reciprocal, w1, f, n);
        jet_to_derivs(f_reciprocal, n);
        mpfr_div(delta, f_reciprocal[n - 2], f_reciprocal[n - 1], RND);
        mpfr_mul_ui(delta, delta, n - 1, RND);
        mpfr_add(x[0], x[0], delta, RND);
        if (++counter > max_it) break;
    }
    mpfr_fprintf(stderr, "%3d %20.12RNe %20.12RNe %20.12RNe\n", counter, x[0], f[0], delta);
    mpfr_clear(delta);
}

void ad_scale (mpfr_t *s, mpfr_t *u, mpfr_t a, int n) {
    assert(s != u);
    assert(sizeof *s == sizeof *u);
    assert(sizeof *a == sizeof (mpfr_t));
    for (int k = 0; k < n; k++) {
        mpfr_mul(s[k], u[k], a, RND);
    }
}

void ad_plus (mpfr_t *p, mpfr_t *u, mpfr_t *v, int n) {
    assert(p != u && p != v);
    assert(sizeof *p == sizeof *u && sizeof *p == sizeof *v);
    for (int k = 0; k < n; k++) {
        mpfr_add(p[k], u[k], v[k], RND);
    }
}

void ad_minus (mpfr_t *p, mpfr_t *u, mpfr_t *v, int n) {
    assert(p != u && p != v);
    assert(sizeof *p == sizeof *u && sizeof *p == sizeof *v);
    for (int k = 0; k < n; k++) {
        mpfr_sub(p[k], u[k], v[k], RND);
    }
}

void ad_square (mpfr_t *s, mpfr_t *u, int n) {
    assert(s != u);
    assert(sizeof *s == sizeof *u);
    for (int k = 0; k < n; k++) {
        t_square(&s[k], u, k);
    }
}

void ad_product (mpfr_t *p, mpfr_t *u, mpfr_t *v, int n) {
    assert(p != u && p != v);
    assert(sizeof *p == sizeof *u && sizeof *p == sizeof *v);
    for (int k = 0; k < n; k++) {
        t_product(&p[k], u, v, k);
    }
}

void ad_quotient (mpfr_t *q, mpfr_t *u, mpfr_t *v, int n) {
    for (int k = 0; k < n; k++) {
        t_quotient(q, u, v, k);
    }
}

void ad_exp (mpfr_t *e, mpfr_t *u, int n) {
    mpfr_t _;
    mpfr_init(_);
    for (int k = 0; k < n; k++) {
        t_exp(e, u, k, &_);
    }
    mpfr_clear(_);
}

void ad_sin_cos (mpfr_t *s, mpfr_t *c, mpfr_t *u, int n) {
    mpfr_t _;
    mpfr_init(_);
    for (int k = 0; k < n; k++) {
        t_sin_cos(s, c, u, k, &_, TRIG);
    }
    mpfr_clear(_);
}

void ad_sinh_cosh (mpfr_t *s, mpfr_t *c, mpfr_t *u, int n) {
    mpfr_t _;
    mpfr_init(_);
    for (int k = 0; k < n; k++) {
        t_sin_cos(s, c, u, k, &_, HYP);
    }
    mpfr_clear(_);
}

void ad_tan_sec2 (mpfr_t *t, mpfr_t *s2, mpfr_t *u, int n) {
    mpfr_t _;
    mpfr_init(_);
    for (int k = 0; k < n; k++) {
        t_tan_sec2(t, s2, u, k, &_, TRIG);
    }
    mpfr_clear(_);
}

void ad_tanh_sech2 (mpfr_t *t, mpfr_t *s2, mpfr_t *u, int n) {
    mpfr_t _;
    mpfr_init(_);
    for (int k = 0; k < n; k++) {
        t_tan_sec2(t, s2, u, k, &_, HYP);
    }
    mpfr_clear(_);
}

void ad_power (mpfr_t *p, mpfr_t *u, mpfr_t a, int n) {
    mpfr_t _1, _2;
    mpfr_inits(_1, _2, NULL);
    for (int k = 0; k < n; k++) {
        t_power(p, u, a, k, &_1, &_2);
    }
    mpfr_clears(_1, _2, NULL);
}

void ad_ln (mpfr_t *l, mpfr_t *u, int n) {
    mpfr_t _;
    mpfr_init(_);
    for (int k = 0; k < n; k++) {
        t_ln(l, u, k, &_);
    }
    mpfr_clear(_);
}

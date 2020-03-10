/*
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"
#include "ad.h"

void set_ad_status (mpfr_t *jet, ad_status s) {
    s == VARIABLE ? mpfr_set_ui(jet[1], 1, RND) : mpfr_set_ui(jet[1], 0, RND);
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

void ad_solver_output (int counter, mpfr_t x, mpfr_t f, mpfr_t delta) {
    mpfr_fprintf(stderr, "%3d %20.12RNe %20.12RNe %20.12RNe\n", counter, x, f, delta);
}

void ad_bisect (model function, mpfr_t *a, mpfr_t *b, int max_it, mpfr_t f_tol, mpfr_t x_tol, mpfr_t *c, mpfr_t *fa, mpfr_t *fc, mode degree) {
    mpfr_t delta, _;
    int counter = 0;
    mpfr_inits(delta, _, NULL);
    function(fa, a, degree + 1);
    jet_to_derivs(fa, degree + 1);
    mpfr_sub(delta, b[0], a[0], RND);
    while(mpfr_cmp_abs(fc[degree], f_tol) >= 0 || mpfr_cmp_abs(delta, x_tol) >= 0) {
        mpfr_add(c[0], a[0], b[0], RND);
        mpfr_div_ui(c[0], c[0], 2, RND);
        function(fc, c, degree + 1);
        jet_to_derivs(fc, degree + 1);
        mpfr_mul(_, fa[degree], fc[degree], RND);
        if (mpfr_sgn(_) < 0) {
            mpfr_set(b[0], c[0], RND);
        } else {
            mpfr_set(a[0], c[0], RND);
        }
        mpfr_sub(delta, b[0], a[0], RND);
        if (++counter > max_it) break;
    }
    ad_solver_output(counter, c[0], fc[degree], delta);
    mpfr_clears(delta, _, NULL);
}

void ad_newton (model function, mpfr_t *f, mpfr_t *x, int max_it, mpfr_t f_tol, mpfr_t x_tol, mode degree) {
    mpfr_t delta;
    int counter = 0;
    mpfr_init_set_ui(delta, 1, RND);
    while(mpfr_cmp_abs(f[degree], f_tol) >= 0 || mpfr_cmp_abs(delta, x_tol) >= 0) {
        function(f, x, degree + 2);
        jet_to_derivs(f, degree + 2);
        mpfr_div(delta, f[degree], f[degree + 1], RND);
        mpfr_sub(x[0], x[0], delta, RND);
        if (++counter > max_it) break;
    }
    ad_solver_output(counter, x[0], f[degree], delta);
    mpfr_clear(delta);
}

void ad_householder (model function, mpfr_t *f, mpfr_t *x, long n, int max_it, mpfr_t f_tol, mpfr_t x_tol, mpfr_t *f_recip, mpfr_t *w1, mode degree) {
    mpfr_t delta;
    int counter = 0;
    mpfr_init_set_ui(delta, 1, RND);
    while(mpfr_cmp_abs(f[degree], f_tol) >= 0 || mpfr_cmp_abs(delta, x_tol) >= 0) {
        function(f, x, n);
        ad_quot(f_recip, w1, f, n);
        jet_to_derivs(f_recip, n);
        mpfr_div(delta, f_recip[n - 2], f_recip[n - 1], RND);
        mpfr_mul_ui(delta, delta, n - 1, RND);
        mpfr_add(x[0], x[0], delta, RND);
        if (++counter > max_it) break;
    }
    ad_solver_output(counter, x[0], f[degree], delta);
    mpfr_clear(delta);
}

mpfr_t *ad_set (mpfr_t *b, mpfr_t *a, int n) {
    for (int k = 0; k < n; k++) {
        mpfr_set(b[k], a[k], RND);
    }
    return b;
}

mpfr_t *ad_scale (mpfr_t *s, mpfr_t *u, mpfr_t a, int n) {
    for (int k = 0; k < n; k++) {
        mpfr_mul(s[k], u[k], a, RND);
    }
    return s;
}

mpfr_t *ad_plus (mpfr_t *p, mpfr_t *u, mpfr_t *v, int n) {
    for (int k = 0; k < n; k++) {
        mpfr_add(p[k], u[k], v[k], RND);
    }
    return p;
}

mpfr_t *ad_minus (mpfr_t *p, mpfr_t *u, mpfr_t *v, int n) {
    for (int k = 0; k < n; k++) {
        mpfr_sub(p[k], u[k], v[k], RND);
    }
    return p;
}

mpfr_t *ad_abs (mpfr_t *a, mpfr_t *u, int n) {
    for (int k = 0; k < n; k++) {
        t_abs(&a[k], u, k);
    }
    return a;
}

mpfr_t *ad_prod (mpfr_t *p, mpfr_t *u, mpfr_t *v, int n) {
    for (int k = 0; k < n; k++) {
        t_prod(&p[k], u, v, k);
    }
    return p;
}

mpfr_t *ad_quot (mpfr_t *q, mpfr_t *u, mpfr_t *v, int n) {
    for (int k = 0; k < n; k++) {
        t_quot(q, u, v, k);
    }
    return q;
}

mpfr_t *ad_sqr (mpfr_t *s, mpfr_t *u, int n) {
    for (int k = 0; k < n; k++) {
        t_sqr(&s[k], u, k);
    }
    return s;
}

mpfr_t *ad_sqrt (mpfr_t *r, mpfr_t *u, int n) {
    for (int k = 0; k < n; k++) {
        t_sqrt(r, u, k);
    }
    return r;
}

mpfr_t *ad_exp (mpfr_t *e, mpfr_t *u, int n) {
    mpfr_t _;
    mpfr_init(_);
    for (int k = 0; k < n; k++) {
        t_exp(e, u, k, &_);
    }
    mpfr_clear(_);
    return e;
}

tuple ad_sin_cos (mpfr_t *s, mpfr_t *c, mpfr_t *u, int n, geometry g) {
    mpfr_t _;
    mpfr_init(_);
    for (int k = 0; k < n; k++) {
        t_sin_cos(s, c, u, k, &_, g);
    }
    mpfr_clear(_);
    return (tuple){s, c};
}

tuple ad_tan_sec2 (mpfr_t *t, mpfr_t *s2, mpfr_t *u, int n, geometry g) {
    mpfr_t _;
    mpfr_init(_);
    for (int k = 0; k < n; k++) {
        t_tan_sec2(t, s2, u, k, &_, g);
    }
    mpfr_clear(_);
    return (tuple){t, s2};
}

mpfr_t *ad_pwr (mpfr_t *p, mpfr_t *u, double a, int n) {
    mpfr_t _, __;
    mpfr_inits(_, __, NULL);
    for (int k = 0; k < n; k++) {
        t_pwr(p, u, a, k, &_, &__);
    }
    mpfr_clears(_, __, NULL);
    return p;
}

mpfr_t *ad_ln (mpfr_t *l, mpfr_t *u, int n) {
    mpfr_t _;
    mpfr_init(_);
    for (int k = 0; k < n; k++) {
        t_ln(l, u, k, &_);
    }
    mpfr_clear(_);
    return l;
}

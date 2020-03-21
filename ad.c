/*
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"
#include "ad.h"

static mpfr_t delta;

void ad_tempvars (void) {
    t_tempvars();
    mpfr_init(delta);
}

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

void ad_newton (model function, mpfr_t *f, mpfr_t *x, int max_it, mpfr_t f_tol, mpfr_t x_tol, mode degree) {
    int counter = 0;
    mpfr_set_ui(delta, 1, RND);
    while(mpfr_cmp_abs(f[degree], f_tol) >= 0 || mpfr_cmp_abs(delta, x_tol) >= 0) {
        function(f, x, degree + 2);
        jet_to_derivs(f, degree + 2);
        mpfr_div(delta, f[degree], f[degree + 1], RND);
        mpfr_sub(x[0], x[0], delta, RND);
        if (++counter > max_it) break;
    }
    mpfr_fprintf(stderr, "%3d %11.3RNe %11.3RNe %11.3RNe ", counter, x[0], &delta, f[degree]);
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

mpfr_t *ad_minus (mpfr_t *m, mpfr_t *u, mpfr_t *v, int n) {
    for (int k = 0; k < n; k++) {
        mpfr_sub(m[k], u[k], v[k], RND);
    }
    return m;
}

mpfr_t *ad_neg (mpfr_t *m, mpfr_t *u, int n) {
    for (int k = 0; k < n; k++) {
        mpfr_neg(m[k], u[k], RND);
    }
    return m;
}

mpfr_t *ad_abs (mpfr_t *a, mpfr_t *u, int n) {
    for (int k = 0; k < n; k++) {
        mpfr_set(a[k], *t_abs(u, k), RND);
    }
    return a;
}

mpfr_t *ad_prod (mpfr_t *p, mpfr_t *u, mpfr_t *v, int n) {
    for (int k = 0; k < n; k++) {
        mpfr_set(p[k], *t_prod(u, v, k), RND);
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
        mpfr_set(s[k], *t_sqr(u, k), RND);
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
    for (int k = 0; k < n; k++) {
        t_exp(e, u, k);
    }
    return e;
}

tuple ad_sin_cos (mpfr_t *s, mpfr_t *c, mpfr_t *u, int n, geometry g) {
    for (int k = 0; k < n; k++) {
        t_sin_cos(s, c, u, k, g);
    }
    return (tuple){s, c};
}

tuple ad_tan_sec2 (mpfr_t *t, mpfr_t *s2, mpfr_t *u, int n, geometry g) {
    for (int k = 0; k < n; k++) {
        t_tan_sec2(t, s2, u, k, g);
    }
    return (tuple){t, s2};
}

mpfr_t *ad_pwr (mpfr_t *p, mpfr_t *u, double a, int n) {
    for (int k = 0; k < n; k++) {
        t_pwr(p, u, a, k);
    }
    return p;
}

mpfr_t *ad_ln (mpfr_t *l, mpfr_t *u, int n) {
    for (int k = 0; k < n; k++) {
        t_ln(l, u, k);
    }
    return l;
}

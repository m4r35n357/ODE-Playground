/*
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"
#include "ad.h"

static int order;

static mpfr_t delta;

void ad_tempvars (int n) {
    t_tempvars(9);
    order = n;
    mpfr_init(delta);
}

series ad_series_c (int n, mpfr_t value) {
    series s = t_jet(n);
    mpfr_set(s[0], value, RND);
    return s;
}

series ad_series_v (int n, mpfr_t value) {
    assert(n > 1);
    series s = ad_series_c(n, value);
    mpfr_set_ui(s[1], 1, RND);
    return s;
}

void jet_output (series s, char* f_colour, char *fk_colour) {
    mpfr_printf("%s%9.6RNf ", f_colour, s[0]);
    for (int i = 1; i < order; i++) {
        mpfr_printf("%s%9.6RNf ", fk_colour, s[i]);
    }
    printf("%s\n", f_colour);
}

void jet_to_derivs (series s) {
    for (int i = 1, fac = 1; i < order; i++) {
        fac *= i;
        mpfr_mul_si(s[i], s[i], fac, RND);
    }
}

void derivative_output (series jet, char* f_colour, char *fk_colour) {
    jet_to_derivs(jet);
    jet_output(jet, f_colour, fk_colour);
}

series ad_set (series b, series a) {
    for (int k = 0; k < order; k++) {
        mpfr_set(b[k], a[k], RND);
    }
    return b;
}

series ad_scale (series s, series u, mpfr_t a) {
    for (int k = 0; k < order; k++) {
        mpfr_mul(s[k], u[k], a, RND);
    }
    return s;
}

series ad_plus (series p, series u, series v) {
    for (int k = 0; k < order; k++) {
        mpfr_add(p[k], u[k], v[k], RND);
    }
    return p;
}

series ad_minus (series m, series u, series v) {
    for (int k = 0; k < order; k++) {
        mpfr_sub(m[k], u[k], v[k], RND);
    }
    return m;
}

series ad_neg (series m, series u) {
    for (int k = 0; k < order; k++) {
        mpfr_neg(m[k], u[k], RND);
    }
    return m;
}

series ad_abs (series a, series u) {
    for (int k = 0; k < order; k++) {
        mpfr_set(a[k], *t_abs(u, k), RND);
    }
    return a;
}

series ad_prod (series p, series u, series v) {
    for (int k = 0; k < order; k++) {
        mpfr_set(p[k], *t_prod(u, v, k), RND);
    }
    return p;
}

series ad_quot (series q, series u, series v) {
    for (int k = 0; k < order; k++) {
        t_quot(q, u, v, k);
    }
    return q;
}

series ad_inv (series i, series v) {
    for (int k = 0; k < order; k++) {
        t_inv(i, v, k);
    }
    return i;
}

series ad_sqr (series s, series u) {
    for (int k = 0; k < order; k++) {
        mpfr_set(s[k], *t_sqr(u, k), RND);
    }
    return s;
}

series ad_sqrt (series r, series u) {
    for (int k = 0; k < order; k++) {
        t_sqrt(r, u, k);
    }
    return r;
}

series ad_exp (series e, series u) {
    for (int k = 0; k < order; k++) {
        t_exp(e, u, k);
    }
    return e;
}

pair ad_sin_cos (series s, series c, series u, geometry g) {
    for (int k = 0; k < order; k++) {
        t_sin_cos(s, c, u, k, g);
    }
    return (pair){s, c};
}

pair ad_tan_sec2 (series t, series s2, series u, geometry g) {
    for (int k = 0; k < order; k++) {
        t_tan_sec2(t, s2, u, k, g);
    }
    return (pair){t, s2};
}

series ad_pwr (series p, series u, mpfr_t a) {
    for (int k = 0; k < order; k++) {
        t_pwr(p, u, a, k);
    }
    return p;
}

series ad_ln (series l, series u) {
    for (int k = 0; k < order; k++) {
        t_ln(l, u, k);
    }
    return l;
}

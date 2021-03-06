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

series ad_series_c (int n, mpfr_t value) {
    series s = t_series(n);
    mpfr_set(s.jet[0], value, RND);
    return s;
}

series ad_series_v (int n, mpfr_t value) {
    assert(n > 1);
    series s = ad_series_c(n, value);
    mpfr_set_ui(s.jet[1], 1, RND);
    return s;
}

void jet_output (series s, long n, char* f_colour, char *fk_colour) {
    mpfr_printf("%s%9.6RNf ", f_colour, s.jet[0]);
    for (int i = 1; i < n; i++) {
        mpfr_printf("%s%9.6RNf ", fk_colour, s.jet[i]);
    }
    printf("%s\n", f_colour);
}

void jet_to_derivs (series s, long n) {
    for (int i = 1, fac = 1; i < n; i++) {
        fac *= i;
        mpfr_mul_ui(s.jet[i], s.jet[i], fac, RND);
    }
}

void derivative_output (series jet, long n, char* f_colour, char *fk_colour) {
    jet_to_derivs(jet, n);
    jet_output(jet, n, f_colour, fk_colour);
}

series ad_set (series b, series a) {
    assert(b.size == a.size);
    for (int k = 0; k < b.size; k++) {
        mpfr_set(b.jet[k], a.jet[k], RND);
    }
    return b;
}

series ad_scale (series s, series u, mpfr_t a) {
    for (int k = 0; k < u.size; k++) {
        mpfr_mul(s.jet[k], u.jet[k], a, RND);
    }
    return s;
}

series ad_plus (series p, series u, series v) {
    for (int k = 0; k < u.size; k++) {
        mpfr_add(p.jet[k], u.jet[k], v.jet[k], RND);
    }
    return p;
}

series ad_minus (series m, series u, series v) {
    for (int k = 0; k < u.size; k++) {
        mpfr_sub(m.jet[k], u.jet[k], v.jet[k], RND);
    }
    return m;
}

series ad_neg (series m, series u) {
    for (int k = 0; k < u.size; k++) {
        mpfr_neg(m.jet[k], u.jet[k], RND);
    }
    return m;
}

series ad_abs (series a, series u) {
    for (int k = 0; k < u.size; k++) {
        mpfr_set(a.jet[k], *t_abs(u, k), RND);
    }
    return a;
}

series ad_prod (series p, series u, series v) {
    for (int k = 0; k < u.size; k++) {
        mpfr_set(p.jet[k], *t_prod(u, v, k), RND);
    }
    return p;
}

series ad_quot (series q, series u, series v) {
    for (int k = 0; k < u.size; k++) {
        t_quot(q, u, v, k);
    }
    return q;
}

series ad_inv (series i, series v) {
    for (int k = 0; k < v.size; k++) {
        t_inv(i, v, k);
    }
    return i;
}

series ad_sqr (series s, series u) {
    for (int k = 0; k < u.size; k++) {
        mpfr_set(s.jet[k], *t_sqr(u, k), RND);
    }
    return s;
}

series ad_sqrt (series r, series u) {
    for (int k = 0; k < u.size; k++) {
        t_sqrt(r, u, k);
    }
    return r;
}

series ad_exp (series e, series u) {
    for (int k = 0; k < u.size; k++) {
        t_exp(e, u, k);
    }
    return e;
}

tuple ad_sin_cos (series s, series c, series u, geometry g) {
    for (int k = 0; k < u.size; k++) {
        t_sin_cos(s, c, u, k, g);
    }
    return (tuple){s.jet, c.jet, u.size};
}

tuple ad_tan_sec2 (series t, series s2, series u, geometry g) {
    for (int k = 0; k < u.size; k++) {
        t_tan_sec2(t, s2, u, k, g);
    }
    return (tuple){t.jet, s2.jet, u.size};
}

series ad_pwr (series p, series u, mpfr_t a) {
    for (int k = 0; k < u.size; k++) {
        t_pwr(p, u, a, k);
    }
    return p;
}

series ad_ln (series l, series u) {
    for (int k = 0; k < u.size; k++) {
        t_ln(l, u, k);
    }
    return l;
}

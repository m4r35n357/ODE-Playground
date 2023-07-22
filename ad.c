/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "taylor-ode.h"
#include "ad.h"

static int order;

void ad_init (int n) {
    order = n;
}

series ad_const (mpfr_t a) {
    return t_const(order, a);
}

series ad_scale (series s, series u, mpfr_t a) {
    for (int k = 0; k < order; k++) {
        mpfr_mul(s[k], u[k], a, RND);
    }
    return s;
}

series ad_add (series p, series u, series v) {
    for (int k = 0; k < order; k++) {
        mpfr_add(p[k], u[k], v[k], RND);
    }
    return p;
}

series ad_sub (series m, series u, series v) {
    for (int k = 0; k < order; k++) {
        mpfr_sub(m[k], u[k], v[k], RND);
    }
    return m;
}

series ad_abs (series a, series u) {
    for (int k = 0; k < order; k++) {
        mpfr_set(a[k], *t_abs(u, k), RND);
    }
    return a;
}

series ad_mul (series p, series u, series v) {
    for (int k = 0; k < order; k++) {
        mpfr_set(p[k], *t_mul(u, v, k), RND);
    }
    return p;
}

series ad_div (series q, series u, series v) {
    for (int k = 0; k < order; k++) {
        t_div(q, u, v, k);
    }
    return q;
}

series ad_inv (series i, series v) {
    for (int k = 0; k < order; k++) {
        t_div(i, NULL, v, k);
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

pair ad_sin_cos (series s, series c, series u, bool trig) {
    for (int k = 0; k < order; k++) {
        t_sin_cos(s, c, u, k, trig);
    }
    return (pair){ .a = s, .b = c };
}

pair ad_tan_sec2 (series t, series s2, series u, bool trig) {
    for (int k = 0; k < order; k++) {
        t_tan_sec2(t, s2, u, k, trig);
    }
    return (pair){ .a = t, .b = s2 };
}

series ad_pwr (series p, series u, mpfr_t a) {
    for (int k = 0; k < order; k++) {
        t_pwr(p, u, a, k);
    }
    return p;
}

series ad_ipwr (series p, series u, int a) {
    for (int k = 0; k < order; k++) {
        t_ipwr(p, u, a, k);
    }
    return p;
}

series ad_ln (series l, series u) {
    for (int k = 0; k < order; k++) {
        t_ln(l, u, k);
    }
    return l;
}

void ad_asin (series as, series du_df, series u, bool trig) {
    for (int k = 0; k < order; k++) {
        t_asin(as, du_df, u, k, trig);
    }
}

void ad_acos (series ac, series du_df, series u, bool trig) {
    for (int k = 0; k < order; k++) {
        t_acos(ac, du_df, u, k, trig);
    }
}

void ad_atan (series at, series du_df, series u, bool trig) {
    for (int k = 0; k < order; k++) {
        t_atan(at, du_df, u, k, trig);
    }
}

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

void ad_init (int n) {
    order = n;
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

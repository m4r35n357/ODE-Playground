/*
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"
#include "ad.h"

static int order;

void ad_init (int n) {
    order = n;
}

void jet_output (series s, char* f_colour, char *fk_colour) {
    printf("%s%9.6Lf ", f_colour, s[0]);
    for (int i = 1; i < order; i++) {
        printf("%s%9.6Lf ", fk_colour, s[i]);
    }
    printf("%s\n", f_colour);
}

void jet_to_derivs (series s) {
    for (int i = 1, fac = 1; i < order; i++) {
        fac *= i;
        s[i] *= fac;
    }
}

void derivative_output (series jet, char* f_colour, char *fk_colour) {
    jet_to_derivs(jet);
    jet_output(jet, f_colour, fk_colour);
}

series ad_set (series b, series a) {
    for (int k = 0; k < order; k++) {
        b[k] = a[k];
    }
    return b;
}

series ad_scale (series s, series u, real a) {
    for (int k = 0; k < order; k++) {
        s[k] = u[k] * a;
    }
    return s;
}

series ad_plus (series p, series u, series v) {
    for (int k = 0; k < order; k++) {
        p[k] = u[k] + v[k];
    }
    return p;
}

series ad_minus (series m, series u, series v) {
    for (int k = 0; k < order; k++) {
        m[k] = u[k] - v[k];
    }
    return m;
}

series ad_neg (series m, series u) {
    for (int k = 0; k < order; k++) {
        m[k] = - u[k];
    }
    return m;
}

series ad_abs (series a, series u) {
    for (int k = 0; k < order; k++) {
        a[k] = t_abs(u, k);
    }
    return a;
}

series ad_prod (series p, series u, series v) {
    for (int k = 0; k < order; k++) {
        p[k] = t_prod(u, v, k);
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
        s[k] = t_sqr(u, k);
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

void ad_sin_cos (series s, series c, series u, geometry g) {
    for (int k = 0; k < order; k++) {
        t_sin_cos(s, c, u, k, g);
    }
}

void ad_tan_sec2 (series t, series s2, series u, geometry g) {
    for (int k = 0; k < order; k++) {
        t_tan_sec2(t, s2, u, k, g);
    }
}

series ad_pwr (series p, series u, real a) {
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

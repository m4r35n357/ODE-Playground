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

series t_jet_v (int n, mpfr_t value) {
    assert(n > 1);
    series s = t_jet_c(n, value);
    mpfr_set_ui(s.a[1], 1, RND);
    return s;
}

void jet_output (series jet, long n, char* f_colour, char *fk_colour) {
    mpfr_printf("%s%9.6RNf ", f_colour, jet.a[0]);
    for (int i = 1; i < n; i++) {
        mpfr_printf("%s%9.6RNf ", fk_colour, jet.a[i]);
    }
    printf("%s\n", f_colour);
}

void jet_to_derivs (series jet, long n) {
    for (int i = 1, fac = 1; i < n; i++) {
        fac *= i;
        mpfr_mul_ui(jet.a[i], jet.a[i], fac, RND);
    }
}

void derivative_output (series jet, long n, char* f_colour, char *fk_colour) {
    jet_to_derivs(jet, n);
    jet_output(jet, n, f_colour, fk_colour);
}

void ad_newton (model function, series f, series x, int max_it, mpfr_t f_tol, mpfr_t x_tol, mode degree) {
    int counter = 0;
    mpfr_set_ui(delta, 1, RND);
    while(mpfr_cmp_abs(f.a[degree], f_tol) >= 0 || mpfr_cmp_abs(delta, x_tol) >= 0) {
        function(f, x);
        jet_to_derivs(f, NEWTON + degree);
        mpfr_div(delta, f.a[degree], f.a[degree + 1], RND);
        mpfr_sub(x.a[0], x.a[0], delta, RND);
        if (++counter > max_it) break;
    }
    mpfr_fprintf(stderr, "%3d %11.3RNe %11.3RNe %11.3RNe ", counter, x.a[0], &delta, f.a[degree]);
}

series ad_set (series b, series a) {
    assert(b.size == a.size);
    for (int k = 0; k < b.size; k++) {
        mpfr_set(b.a[k], a.a[k], RND);
    }
    return b;
}

series ad_scale (series s, series u, mpfr_t a) {
    for (int k = 0; k < u.size; k++) {
        mpfr_mul(s.a[k], u.a[k], a, RND);
    }
    return s;
}

series ad_plus (series p, series u, series v) {
    for (int k = 0; k < u.size; k++) {
        mpfr_add(p.a[k], u.a[k], v.a[k], RND);
    }
    return p;
}

series ad_minus (series m, series u, series v) {
    for (int k = 0; k < u.size; k++) {
        mpfr_sub(m.a[k], u.a[k], v.a[k], RND);
    }
    return m;
}

series ad_neg (series m, series u) {
    for (int k = 0; k < u.size; k++) {
        mpfr_neg(m.a[k], u.a[k], RND);
    }
    return m;
}

series ad_abs (series a, series u) {
    for (int k = 0; k < u.size; k++) {
        mpfr_set(a.a[k], *t_abs(u, k), RND);
    }
    return a;
}

series ad_prod (series p, series u, series v) {
    for (int k = 0; k < u.size; k++) {
        mpfr_set(p.a[k], *t_prod(u, v, k), RND);
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
        mpfr_set(s.a[k], *t_sqr(u, k), RND);
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
    return (tuple){s.a, c.a, u.size};
}

tuple ad_tan_sec2 (series t, series s2, series u, geometry g) {
    for (int k = 0; k < u.size; k++) {
        t_tan_sec2(t, s2, u, k, g);
    }
    return (tuple){t.a, s2.a, u.size};
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

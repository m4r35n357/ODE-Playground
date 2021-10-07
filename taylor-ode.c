/*
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

const int BASE = 10;

const mpfr_rnd_t RND = MPFR_RNDN;

static char fs[42];

static mpfr_t D1, D_1, D2, D_2, _horner, _, _fk, _const, _abs, _prod, _sqr;

void t_init (int dp) {
    if (dp == 0) {
        sprintf(fs, "%%.RNe %%.RNe %%.RNe %%.9RNe\n");
    } else {
        sprintf(fs, "%%+.%uRNe %%+.%uRNe %%+.%uRNe %%+.9RNe\n", dp, dp, dp);
    }
    mpfr_inits(_horner, _, _fk, _const, _abs, _prod, _sqr, NULL);
    mpfr_init_set_ui(D1, 1, RND);
    mpfr_init_set_si(D_1, -1, RND);
    mpfr_init_set_ui(D2, 2, RND);
    mpfr_init_set_si(D_2, -2, RND);
}

void t_output (mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t h, int step) {
    mpfr_mul_si(_, h, step, RND);
    mpfr_printf(fs, x, y, z, _);
}

void t_params (char **argv, int argc, ...) {
    va_list model;
    va_start(model, argc);
    for (int i = 9; i < argc; i++) {
        mpfr_init_set_str(*va_arg(model, mpfr_t *), argv[i], BASE, RND);
    }
    va_end(model);
}

series t_jet (int n) {
    mpfr_t *s = calloc((size_t)n, sizeof (mpfr_t));
    if (s == NULL) exit(2);
    for (int i = 0; i < n; i++) {
        mpfr_init_set_ui(s[i], 0, RND);
    }
    return s;
}

mpfr_t *t_horner (series s, int n, mpfr_t h) {
    mpfr_set_zero(_horner, 1);
    for (int i = n; i >= 0; i--) {
        mpfr_fma(_horner, _horner, h, s[i], RND);
    }
    if (mpfr_number_p(_horner) == 0) exit(1);
    return &_horner;
}

mpfr_t *t_const (mpfr_t value, int k){
    if (k == 0) {
        mpfr_set(_const, value, RND);
    } else {
        mpfr_set_zero(_const, 1);
    }
    return &_const;
}

mpfr_t *t_abs (series u, int k) {
    if (mpfr_sgn(u[0]) < 0) {
        mpfr_neg(_abs, u[k], RND);
    } else {
        mpfr_set(_abs, u[k], RND);
    }
    return &_abs;
}

static mpfr_t *cauchy (series a, series b, int k, int j_lower, int j_upper) {
    mpfr_set_zero(_, 1);
    for (int j = j_lower; j <= j_upper; j++) {
        mpfr_fma(_, a[j], b[k - j], _, RND);
    }
    return &_;
}

mpfr_t *t_prod (series u, series v, int k) {
    mpfr_set(_prod, *cauchy(u, v, k, 0, k), RND);
    return &_prod;
}

mpfr_t *t_sqr (series u, int k) {
    mpfr_set(_sqr, *cauchy(u, u, k, 0, k), RND);
    return &_sqr;
}

mpfr_t *t_quot (series q, series u, series v, int k) {
    assert(mpfr_zero_p(v[0]) == 0);
    assert(q != u && q != v);
    if (k == 0) {
        mpfr_set(q[0], u[0], RND);
    } else {
        mpfr_sub(q[k], u[k], *cauchy(q, v, k, 0, k - 1), RND);
    }
    mpfr_div(q[k], q[k], v[0], RND);
    return &q[k];
}

mpfr_t *t_inv (series i, series v, int k) {
    assert(mpfr_zero_p(v[0]) == 0);
    assert(i != v);
    if (k == 0) {
        mpfr_set(i[0], D1, RND);
    } else {
        mpfr_neg(i[k], *cauchy(i, v, k, 0, k - 1), RND);
    }
    mpfr_div(i[k], i[k], v[0], RND);
    return &i[k];
}

mpfr_t *t_sqrt (series r, series u, int k) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(r != u);
    if (k == 0) {
        mpfr_sqrt(r[0], u[0], RND);
    } else {
        mpfr_sub(r[k], u[k], *cauchy(r, r, k, 1, k - 1), RND);
        mpfr_div_2ui(r[k], r[k], 1, RND);
        mpfr_div(r[k], r[k], r[0], RND);
    }
    return &r[k];
}

static mpfr_t *f_k (series df_du, series u, int k, int j_lower, int j_upper) {
    mpfr_set_zero(_, 1);
    for (int j = j_lower; j <= j_upper; j++) {
        mpfr_mul_si(_fk, u[k - j], k - j, RND);
        mpfr_fma(_, df_du[j], _fk, _, RND);
    }
    mpfr_div_si(_, _, k, RND);
    return &_;
}

mpfr_t *t_exp (series e, series u, int k) {
    assert(e != u);
    if (k == 0) {
        mpfr_exp(e[0], u[0], RND);
    } else {
        mpfr_set(e[k], *f_k(e, u, k, 0, k - 1), RND);
    }
    return &e[k];
}

pair t_sin_cos (series s, series c, series u, int k, geometry g) {
    assert(s != c && s != u && c != u);
    if (k == 0) {
        g == TRIG ? mpfr_sin_cos(s[0], c[0], u[0], RND) : mpfr_sinh_cosh(s[0], c[0], u[0], RND);
    } else {
        mpfr_set(s[k], *f_k(c, u, k, 0, k - 1), RND);
        mpfr_mul(c[k], *f_k(s, u, k, 0, k - 1), g == TRIG ? D_1 : D1, RND);
    }
    return (pair){&s[k], &c[k]};
}

pair t_tan_sec2 (series t, series s, series u, int k, geometry g) {
    assert(t != s && t != u && s != u);
    if (k == 0) {
        g == TRIG ? mpfr_tan(t[0], u[0], RND) : mpfr_tanh(t[0], u[0], RND);
        mpfr_sqr(s[0], t[0], RND);
        g == TRIG ? mpfr_add(s[0], D1, s[0], RND) : mpfr_sub(s[0], D1, s[0], RND);
    } else {
        mpfr_set(t[k], *f_k(s, u, k, 0, k - 1), RND);
        mpfr_mul(s[k], *f_k(t, t, k, 0, k - 1), g == TRIG ? D2 : D_2, RND);
    }
    return (pair){&t[k], &s[k]};
}

mpfr_t *t_pwr (series p, series u, mpfr_t powr, int k) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(p != u);
    if (k == 0) {
        mpfr_pow(p[0], u[0], powr, RND);
    } else {
        mpfr_mul(p[k], powr, *f_k(p, u, k, 0, k - 1), RND);
        mpfr_sub(p[k], p[k], *f_k(u, p, k, 1, k - 1), RND);
        mpfr_div(p[k], p[k], u[0], RND);
    }
    return &p[k];
}

mpfr_t *t_ln (series l, series u, int k) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(l != u);
    if (k == 0) {
        mpfr_log(l[0], u[0], RND);
    } else {
        mpfr_sub(l[k], u[k], *f_k(u, l, k, 1, k - 1), RND);
        mpfr_div(l[k], l[k], u[0], RND);
    }
    return &l[k];
}

void tsm (int argc, char **argv, int n, mpfr_t h, int steps, mpfr_t x0, mpfr_t y0, mpfr_t z0) {
    series x = t_jet(n + 1); mpfr_set(x[0], x0, RND);
    series y = t_jet(n + 1); mpfr_set(y[0], y0, RND);
    series z = t_jet(n + 1); mpfr_set(z[0], z0, RND);
    void *p = get_p(argc, argv, n);
    components *c_dot = malloc(sizeof (components));
    mpfr_inits(c_dot->x, c_dot->y, c_dot->z, NULL);
    for (int step = 0; step < steps; step++) {
        for (int k = 0; k < n; k++) {
            ode(x, y, z, c_dot, p, k);
            mpfr_div_si(x[k + 1], c_dot->x, k + 1, RND);
            mpfr_div_si(y[k + 1], c_dot->y, k + 1, RND);
            mpfr_div_si(z[k + 1], c_dot->z, k + 1, RND);
        }
        t_output(x[0], y[0], z[0], h, step);
        mpfr_set(x[0], *t_horner(x, n, h), RND);
        mpfr_set(y[0], *t_horner(y, n, h), RND);
        mpfr_set(z[0], *t_horner(z, n, h), RND);
    }
    t_output(x[0], y[0], z[0], h, steps);
}

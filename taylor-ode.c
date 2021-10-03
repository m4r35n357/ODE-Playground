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

static mpfr_t D0, D1, D_1, D2, D_2, _, __, ___;

void t_tempvars (long dp) {
    dp == 0 ? sprintf(fs, "%%.RNe %%.RNe %%.RNe %%.9RNe\n") : sprintf(fs, "%%+.%luRNe %%+.%luRNe %%+.%luRNe %%+.9RNe\n", dp, dp, dp);
    mpfr_inits(_, __, ___, NULL);
    mpfr_init_set_ui(D0, 0, RND);
    mpfr_init_set_ui(D1, 1, RND);
    mpfr_init_set_si(D_1, -1, RND);
    mpfr_init_set_ui(D2, 2, RND);
    mpfr_init_set_si(D_2, -2, RND);
}

void t_output (mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t h, long step) {
    mpfr_mul_si(_, h, step, RND);
    mpfr_printf(fs, x, y, z, _);
}

void t_params (char **argv, int argc, ...) {
    va_list model_params;
    va_start(model_params, argc);
    for (int i = 9; i < argc; i++) {
        mpfr_init_set_str(*va_arg(model_params, mpfr_t *), argv[i], BASE, RND);
    }
    va_end(model_params);
}

series t_jet (long n) {
    mpfr_t *s = calloc((size_t)n, sizeof (mpfr_t));
    if (s == NULL) exit(2);
    for (long i = 0; i < n; i++) {
        mpfr_init_set_ui(s[i], 0, RND);
    }
    return s;
}

mpfr_t *t_horner (series s, long n, mpfr_t h) {
    mpfr_set_zero(_, 1);
    for (long i = n; i >= 0; i--) {
        mpfr_fma(_, _, h, s[i], RND);
    }
    if (mpfr_number_p(_) == 0) exit(1);
    return &_;
}

mpfr_t *t_const (mpfr_t *value, int k){
    mpfr_set_zero(_, 1);
    return k == 0 ? value : &_;
}

mpfr_t *t_abs (series u, int k) {
    mpfr_sgn(u[0]) < 0 ? mpfr_neg(_, u[k], RND) : mpfr_set(_, u[k], RND);
    return &_;
}

static mpfr_t *cauchy (mpfr_t *ck, series a, series b, int k, int j_lower, int j_upper) {
    mpfr_set_zero(*ck, 1);
    for (int j = j_lower; j <= j_upper; j++) {
        mpfr_fma(*ck, a[j], b[k - j], *ck, RND);
    }
    return ck;
}

mpfr_t *t_prod (series u, series v, int k) {
    return cauchy(&__, u, v, k, 0, k);
}

mpfr_t *t_sqr (series u, int k) {
    return cauchy(&__, u, u, k, 0, k);
}

mpfr_t *t_quot (series q, series u, series v, int k) {
    assert(mpfr_zero_p(v[0]) == 0);
    assert(q != u && q != v);
    k == 0 ? mpfr_set(q[k], u[0], RND) : mpfr_sub(q[k], u[k], *cauchy(&__, q, v, k, 0, k - 1), RND);
    mpfr_div(q[k], q[k], v[0], RND);
    return &q[k];
}

mpfr_t *t_inv (series i, series v, int k) {
    assert(mpfr_zero_p(v[0]) == 0);
    assert(i != v);
    k == 0 ? mpfr_set(i[k], D1, RND) : mpfr_neg(i[k], *cauchy(&__, i, v, k, 0, k - 1), RND);
    mpfr_div(i[k], i[k], v[0], RND);
    return &i[k];
}

mpfr_t *t_sqrt (series r, series u, int k) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(r != u);
    if (k == 0) {
        mpfr_sqrt(r[k], u[0], RND);
    } else {
        mpfr_sub(r[k], u[k], *cauchy(&__, r, r, k, 1, k - 1), RND);
        mpfr_div_2ui(r[k], r[k], 1, RND);
        mpfr_div(r[k], r[k], r[0], RND);
    }
    return &r[k];
}

static mpfr_t *f_k (mpfr_t *fk, series df_du, series u, int k, int j_lower, int j_upper, mpfr_t factor) {
    mpfr_set_zero(*fk, 1);
    for (int j = j_lower; j <= j_upper; j++) {
        mpfr_mul_si(_, u[k - j], k - j, RND);
        mpfr_fma(*fk, df_du[j], _, *fk, RND);
    }
    mpfr_div_si(*fk, *fk, k, RND);
    mpfr_mul(*fk, *fk, factor, RND);
    return fk;
}

mpfr_t *t_exp (series e, series u, int k) {
    assert(e != u);
    if (k == 0) {
        mpfr_exp(e[k], u[0], RND);
    } else {
        f_k(&e[k], e, u, k, 0, k - 1, D1);
    }
    return &e[k];
}

pair t_sin_cos (series s, series c, series u, int k, geometry g) {
    assert(s != c && s != u && c != u);
    if (k == 0) {
        g == TRIG ? mpfr_sin_cos(s[k], c[k], u[0], RND) : mpfr_sinh_cosh(s[k], c[k], u[0], RND);
    } else {
        f_k(&s[k], c, u, k, 0, k - 1, D1);
        f_k(&c[k], s, u, k, 0, k - 1, g == TRIG ? D_1 : D1);
    }
    return (pair){&s[k], &c[k]};
}

pair t_tan_sec2 (series t, series s, series u, int k, geometry g) {
    assert(t != s && t != u && s != u);
    if (k == 0) {
        g == TRIG ? mpfr_tan(t[k], u[0], RND) : mpfr_tanh(t[k], u[0], RND);
        g == TRIG ? mpfr_sec(s[k], u[0], RND) : mpfr_sech(s[k], u[0], RND);
        mpfr_sqr(s[k], s[k], RND);
    } else {
        f_k(&t[k], s, u, k, 0, k - 1, D1);
        f_k(&s[k], t, t, k, 0, k - 1, g == TRIG ? D2 : D_2);
    }
    return (pair){&t[k], &s[k]};
}

mpfr_t *t_pwr (series p, series u, mpfr_t a, int k) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(p != u);
    if (k == 0) {
        mpfr_pow(p[k], u[0], a, RND);
    } else {
        mpfr_sub(p[k], *f_k(&__, p, u, k, 0, k - 1, a), *f_k(&___, u, p, k, 1, k - 1, D1), RND);
        mpfr_div(p[k], p[k], u[0], RND);
    }
    return &p[k];
}

mpfr_t *t_ln (series l, series u, int k) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(l != u);
    if (k == 0) {
        mpfr_log(l[k], u[0], RND);
    } else {
        mpfr_sub(l[k], u[k], *f_k(&__, u, l, k, 1, k - 1, D1), RND);
        mpfr_div(l[k], l[k], u[0], RND);
    }
    return &l[k];
}

void tsm (int argc, char **argv, long n, mpfr_t h, long steps, mpfr_t x0, mpfr_t y0, mpfr_t z0) {
    series x = t_jet(n + 1); mpfr_set(x[0], x0, RND);
    series y = t_jet(n + 1); mpfr_set(y[0], y0, RND);
    series z = t_jet(n + 1); mpfr_set(z[0], z0, RND);
    void *p = get_p(argc, argv, n);
    components *c_dot = malloc(sizeof (components));
    mpfr_inits(c_dot->x, c_dot->y, c_dot->z, NULL);
    for (long step = 0; step < steps; step++) {
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

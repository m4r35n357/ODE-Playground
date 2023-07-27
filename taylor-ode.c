/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>
#include <mpfr.h>
#include "taylor-ode.h"

const int BASE = 10;

const mpfr_rnd_t RND = MPFR_RNDN;

static mpfr_t _, _a, _m, _s;

static char template[60];

void t_params (char **argv, int argc, ...) {
    PRINT_ARGS(argc, argv);
    va_list model;
    va_start(model, argc);
    for (int i = 9; i < argc; i++) mpfr_init_set_str(*va_arg(model, mpfr_t *), argv[i], BASE, RND);
    va_end(model);
}

series t_jet (int n) {
    mpfr_t *s = calloc((size_t)n, sizeof (mpfr_t)); CHECK(s);
    for (int i = 0; i <= n; i++) mpfr_init(s[i]);
    return s;
}

mpfr_t *t_const (int n, mpfr_t a) {
    series c = t_jet(n);
    mpfr_set(c[0], a, RND);
    for (int k = 1; k < n; k++) mpfr_set_zero(c[k], 1);
    return c;
}

void t_out (mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t h, int step, clock_t since) {
    mpfr_mul_si(_, h, step, RND);
    mpfr_printf(template, x, y, z, _, (double)(clock() - since) / CLOCKS_PER_SEC);
}

mpfr_t *t_horner (series s, int n, mpfr_t h) {
    mpfr_set_zero(_, 1);
    for (int i = n; i >= 0; i--) mpfr_fma(_, _, h, s[i], RND);
    CHECK(mpfr_number_p(_) != 0);
    return &_;
}

void tsm (int places, int n, mpfr_t h, int steps, mpfr_t x0, mpfr_t y0, mpfr_t z0, parameters *p, clock_t t0) {
    sprintf(template, "%%+.%uRNe %%+.%uRNe %%+.%uRNe %%+.9RNe %%.3f\n", places, places, places);
    mpfr_inits(_, _a, _m, _s, NULL);
    series x = t_const(n + 1, x0), y = t_const(n + 1, y0), z = t_const(n + 1, z0);
    triplet *v_k = malloc(sizeof (triplet)); CHECK(v_k);
    mpfr_inits(v_k->x, v_k->y, v_k->z, NULL);
    for (int step = 0; step < steps; step++) {
        for (int k = 0; k < n; k++) {
            ode(v_k, x, y, z, p, k);
            mpfr_div_si(x[k + 1], v_k->x, k + 1, RND);
            mpfr_div_si(y[k + 1], v_k->y, k + 1, RND);
            mpfr_div_si(z[k + 1], v_k->z, k + 1, RND);
        }
        t_out(x[0], y[0], z[0], h, step, t0);
        mpfr_swap(x[0], *t_horner(x, n, h));
        mpfr_swap(y[0], *t_horner(y, n, h));
        mpfr_swap(z[0], *t_horner(z, n, h));
    }
    t_out(x[0], y[0], z[0], h, steps, t0);
}

static int _half_ (int k) {
    return 1 + (k - (k % 2 ? 1 : 2)) / 2;
}

static mpfr_t *_prod_ (mpfr_t *_p, series u, series v, int k, int k0, int k1) {
    mpfr_set_zero(*_p, 1);
    for (int j = k0; j < k1; j++) {
        mpfr_fma(*_p, u[j], v[k - j], *_p, RND);
    }
    return _p;
}

static mpfr_t *_exp_ (mpfr_t *_e, series df_du, series u, int k, mpfr_t *_du_dt, int factor) {
    mpfr_set_zero(*_e, 1);
    for (int j = 0; j < k; j++) {
        mpfr_mul_si(*_du_dt, u[k - j], k - j, RND);
        mpfr_fma(*_e, df_du[j], *_du_dt, *_e, RND);
    }
    mpfr_div_si(*_e, *_e, k, RND);
    if (factor != 1) mpfr_mul_si(*_e, *_e, factor, RND);
    return _e;
}

static mpfr_t *_log_ (mpfr_t *_l, series f, series du_df, series u, int k, mpfr_t *_df_dt, bool neg) {
    mpfr_set_zero(*_l, 1);
    for (int j = 1; j < k; j++) {
        mpfr_mul_si(*_df_dt, f[k - j], k - j, RND);
        mpfr_fma(*_l, du_df[j], *_df_dt, *_l, RND);
    }
    mpfr_div_si(*_l, *_l, k, RND);
    neg ? mpfr_add(*_l, u[k], *_l, RND) : mpfr_sub(*_l, u[k], *_l, RND);
    mpfr_div(*_l, *_l, du_df[0], RND);
    return _l;
}

mpfr_t *t_abs (series u, int k) {
    CHECK(mpfr_zero_p(u[0]) == 0);
    mpfr_sgn(u[0]) < 0 ? mpfr_neg(_a, u[k], RND) : mpfr_set(_a, u[k], RND);
    return &_a;
}

mpfr_t *t_mul (series u, series v, int k) {
    return _prod_(&_m, u, v, k, 0, k + 1);
}

mpfr_t *t_sqr (series u, int k) {
    mpfr_mul_2si(_s, *_prod_(&_s, u, u, k, 0, _half_(k)), 1, RND);
    if (!(k % 2)) mpfr_fma(_s, u[k / 2], u[k / 2], _s, RND);
    return &_s;
}

mpfr_t *t_div (series q, series u, series v, int k) {
    CHECK(mpfr_zero_p(v[0]) == 0); CHECK(q != u && q != v);
    if (!k) u ? mpfr_set(q[k], u[k], RND) : mpfr_set_si(q[k], 1, RND);
    else {
        _prod_(&q[k], q, v, k, 0, k);
        u ? mpfr_sub(q[k], u[k], q[k], RND) : mpfr_neg(q[k], q[k], RND);
    }
    mpfr_div(q[k], q[k], v[0], RND);
    return &q[k];
}

mpfr_t *t_sqrt (series r, series u, int k) {
    CHECK(mpfr_sgn(u[0]) > 0); CHECK(r != u);
    if (!k) mpfr_sqrt(r[k], u[k], RND);
    else {
        mpfr_mul_2si(r[k], *_prod_(&r[k], r, r, k, 1, _half_(k)), 1, RND);
        if (!(k % 2)) mpfr_fma(r[k], r[k / 2], r[k / 2], r[k], RND);
        mpfr_sub(r[k], u[k], r[k], RND);
        mpfr_div_2si(r[k], r[k], 1, RND);
        mpfr_div(r[k], r[k], r[0], RND);
    }
    return &r[k];
}

mpfr_t *t_exp (series e, series u, int k) {
    CHECK(e != u);
    if (!k) mpfr_exp(e[k], u[k], RND); else return _exp_(&e[k], e, u, k, &_, 1);
    return &e[k];
}

pair t_sin_cos (series s, series c, series u, int k, bool trig) {
    CHECK(s != c && s != u && c != u);
    if (!k) {
        trig ? mpfr_sin_cos(s[k], c[k], u[k], RND) : mpfr_sinh_cosh(s[k], c[k], u[k], RND);
    } else {
        _exp_(&s[k], c, u, k, &_, 1);
        _exp_(&c[k], s, u, k, &_, trig ? -1 : 1);
    };
    return (pair){.a = &s[k], .b = &c[k]};
}

pair t_tan_sec2 (series t, series s, series u, int k, bool trig) {
    CHECK(t != s && t != u && s != u);
    if (!k) {
        trig ? mpfr_tan(t[k], u[k], RND) : mpfr_tanh(t[k], u[k], RND);
        mpfr_sqr(s[k], t[k], RND);
        trig ? mpfr_add_si(s[k], s[k], 1, RND) : mpfr_si_sub(s[k], 1, s[k], RND);
    } else {
        _exp_(&t[k], s, u, k, &_, 1);
        _exp_(&s[k], t, t, k, &_, trig ? 2 : -2);
    };
    return (pair){.a = &t[k], .b = &s[k]};
}

mpfr_t *t_pwr (series p, series u, mpfr_t a, int k) {
    CHECK(mpfr_sgn(u[0]) > 0); CHECK(p != u);
    if (!k) mpfr_pow(p[k], u[k], a, RND);
    else {
        mpfr_set_zero(p[k], 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_si(_, a, k - j, RND);
            mpfr_sub_si(_, _, j, RND);
            mpfr_mul(_, _, u[k - j], RND);
            mpfr_fma(p[k], p[j], _, p[k], RND);
        }
        mpfr_div_si(p[k], p[k], k, RND);
        mpfr_div(p[k], p[k], u[0], RND);
    }
    return &p[k];
}

mpfr_t *t_ln (series l, series u, int k) {
    CHECK(mpfr_sgn(u[0]) > 0); CHECK(l != u);
    if (!k) mpfr_log(l[k], u[k], RND); else _log_(&l[k], l, u, u, k, &_, false);
    return &l[k];
}

pair t_asin (series a, series g, series u, int k, bool trig) {
    CHECK(trig ? mpfr_cmpabs_ui(u[0], 1) <= 0 : true); CHECK(a != g && a != u && g != u);
    if (!k) {
        trig ? mpfr_asin(a[k], u[k], RND) : mpfr_asinh(a[k], u[k], RND);
        mpfr_sqr(g[k], u[k], RND);
        trig ? mpfr_si_sub(g[0], 1, g[k], RND) : mpfr_add_si(g[k], g[k], 1, RND);
        mpfr_sqrt(g[k], g[k], RND);
    } else {
        _log_(&a[k], a, g, u, k, &_, false);
        _exp_(&g[k], u, a, k, &_, trig ? -1 : 1);
    };
    return (pair){.a = &a[k], .b = &g[k]};
}

pair t_acos (series a, series g, series u, int k, bool trig) {
    CHECK(trig ? mpfr_cmpabs_ui(u[0], 1) <= 0 : mpfr_cmp_si(u[0], 1) >= 0); CHECK(a != g && a != u && g != u);
    if (!k) {
        trig ? mpfr_acos(a[k], u[k], RND) : mpfr_acosh(a[k], u[k], RND);
        mpfr_sqr(g[k], u[k], RND);
        trig ? mpfr_si_sub(g[k], 1, g[k], RND) : mpfr_sub_si(g[k], g[k], 1, RND);
        mpfr_sqrt(g[k], g[k], RND);
        if (trig) mpfr_neg(g[k], g[k], RND);
    } else {
        _log_(&a[k], a, g, u, k, &_, trig);
        _exp_(&g[k], u, a, k, &_, 1);
    };
    return (pair){.a = &a[k], .b = &g[k]};
}

pair t_atan (series a, series g, series u, int k, bool trig) {
    CHECK(trig ? true : mpfr_cmpabs_ui(u[0], 1) <= 0); CHECK(a != g && a != u && g != u);
    if (!k) {
        trig ? mpfr_atan(a[k], u[k], RND) : mpfr_atanh(a[k], u[k], RND);
        mpfr_sqr(g[k], u[k], RND);
        trig ? mpfr_add_si(g[k], g[k], 1, RND) : mpfr_si_sub(g[k], 1, g[k], RND);
    } else {
        _log_(&a[k], a, g, u, k, &_, false);
        _exp_(&g[k], u, u, k, &_, trig ? 2 : -2);
    };
    return (pair){.a = &a[k], .b = &g[k]};
}

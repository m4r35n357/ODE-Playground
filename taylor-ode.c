/*
 * Arbitrary-Order Taylor Series Integrator, together with a well-tested set of recurrence relations
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>
#include <mpfr.h>
#include "taylor-ode.h"

static real __, _a, _m, _h, _p, _fk, D0, D1, D_1, D2, D_2;

static char format[60];

void tsm_get_p (char **argv, int argc, ...) {
    va_list _;
    va_start(_, argc);
    for (int i = 9; i < argc; i++) mpfr_init_set_str(*va_arg(_, real *), argv[i], BASE, RND);
    va_end(_);
}

series tsm_jet (int k) {
    CHECK(k > 0);
    series _ = malloc((size_t)k * sizeof (real)); CHECK(_);
    for (int i = 0; i < k; i++) mpfr_init(_[i]);
    return _;
}

series tsm_const (int k, real a) {
    series _ = tsm_jet(k);
    for (int i = 0; i < k; i++) mpfr_set(_[i], i ? D0 : a, RND);
    return _;
}

void tsm_init(int dp) {
    sprintf(format, "%%+.%uRNe %%+.%uRNe %%+.%uRNe %%+.9RNe %%.3f\n", dp, dp, dp);
    mpfr_init_set_si(D0, 0, RND);
    mpfr_init_set_si(D1, 1, RND);
    mpfr_init_set_si(D2, 2, RND);
    mpfr_init_set_si(D_1, -1, RND);
    mpfr_init_set_si(D_2, -2, RND);
}

real *horner (series s, int n, real h) {
    mpfr_set_zero(__, 1);
    for (int i = n; i >= 0; i--) mpfr_fma(__, __, h, s[i], RND);
    CHECK(mpfr_number_p(__) != 0);
    return &__;
}

void _out_ (real x, real y, real z, real h, int step, clock_t since) {
    mpfr_mul_si(__, h, step, RND);
    mpfr_printf(format, x, y, z, __, (double)(clock() - since) / CLOCKS_PER_SEC);
}

void tsm_stdout (int n, real h, int steps, series3 *j, parameters *p, clock_t t0) {
    triplet *v_k = malloc(sizeof (triplet)); CHECK(v_k);
    mpfr_inits(v_k->x, v_k->y, v_k->z, __, _a, _m, _h, _p, _fk, NULL);
    for (int step = 0; step < steps; step++) {
        for (int k = 0; k < n; k++) {
            ode(v_k, j->x, j->y, j->z, p, k);
            mpfr_div_si(j->x[k + 1], v_k->x, k + 1, RND);
            mpfr_div_si(j->y[k + 1], v_k->y, k + 1, RND);
            mpfr_div_si(j->z[k + 1], v_k->z, k + 1, RND);
        }
        _out_(j->x[0], j->y[0], j->z[0], h, step, t0);
        mpfr_swap(j->x[0], *horner(j->x, n, h));
        mpfr_swap(j->y[0], *horner(j->y, n, h));
        mpfr_swap(j->z[0], *horner(j->z, n, h));
    }
    _out_(j->x[0], j->y[0], j->z[0], h, steps, t0);
}

static real *_cauchy_ (real *_, series b, series a, int k, int k0, int k1) {
    mpfr_set_zero(*_, 1);
    for (int j = k0; j <= k1; j++) mpfr_fma(*_, b[j], a[k - j], *_, RND);
    return _;
}

static real *_half_ (series a, int k, int k0) {
    mpfr_mul_2si(_h, *_cauchy_(&_h, a, a, k, k0, (k - (k % 2 ? 1 : 2)) / 2), 1, RND);
    if (!(k % 2)) mpfr_fma(_h, a[k / 2], a[k / 2], _h, RND);
    return &_h;
}

static real *_chain_ (real *_, series b, series a, int k, int k0, real scale) {
    mpfr_set_zero(*_, 1);
    for (int j = k0; j < k; j++) {
        mpfr_mul_si(__, a[k - j], k - j, RND);
        mpfr_fma(*_, b[j], __, *_, RND);
    }
    mpfr_div_si(*_, *_, k, RND);
    mpfr_mul(*_, *_, scale, RND);
    return _;
}

static real *_fk_ (real *_, series df_du, series u, int k, real scale) {
    return _chain_(_, df_du, u, k, 0, scale);
}

static real *_uk_ (real *_, series df_du, series u, int k, real *f_k, real sign) {
    mpfr_sub(*_, *f_k, *_chain_(&_fk, df_du, u, k, 1, sign), RND);
    mpfr_div(*_, *_, df_du[0], RND);
    return _;
}

real *t_abs (series u, int k) {
    CHECK(mpfr_zero_p(u[0]) == 0);
    mpfr_sgn(u[0]) < 0 ? mpfr_neg(_a, u[k], RND) : mpfr_set(_a, u[k], RND);
    return &_a;
}

real *t_mul (series u, series w, int k) {
    return _cauchy_(&_m, u, w, k, 0, k);
}

real *t_div (series q, series u, series w, int k) {
    CHECK(mpfr_zero_p(w[0]) == 0); CHECK(q != u && q != w);
    mpfr_sub(q[k], u[k], !k ? D0 : *_cauchy_(&q[k], q, w, k, 0, k - 1), RND);
    mpfr_div(q[k], q[k], w[0], RND);
    return &q[k];
}

real *t_sqr (series u, int k) {
    return _half_(u, k, 0);
}

real *t_sqrt (series r, series u, int k) {
    CHECK(mpfr_sgn(u[0]) > 0); CHECK(r != u);
    if (!k) mpfr_sqrt(r[k], u[k], RND);
    else {
        mpfr_sub(r[k], u[k], *_half_(r, k, 1), RND);
        mpfr_div_2si(r[k], r[k], 1, RND);
        mpfr_div(r[k], r[k], r[0], RND);
    }
    return &r[k];
}

real *t_exp (series e, series u, int k) {
    CHECK(e != u);
    if (!k) mpfr_exp(e[k], u[k], RND); else _fk_(&e[k], e, u, k, D1);
    return &e[k];
}

pair t_sin_cos (series s, series c, series u, int k, bool trig) {
    CHECK(s != c && s != u && c != u);
    if (!k) trig ? mpfr_sin_cos(s[k], c[k], u[k], RND) : mpfr_sinh_cosh(s[k], c[k], u[k], RND);
    else {
        _fk_(&s[k], c, u, k, D1);
        _fk_(&c[k], s, u, k, trig ? D_1 : D1);
    };
    return (pair){.a = &s[k], .b = &c[k]};
}

pair t_tan_sec2 (series t, series s, series u, int k, bool trig) {
    CHECK(t != s && t != u && s != u);
    if (!k) {
        trig ? mpfr_tan(t[k], u[k], RND) : mpfr_tanh(t[k], u[k], RND);
        trig ? mpfr_sec(s[k], u[k], RND) : mpfr_sech(s[k], u[k], RND);
        mpfr_sqr(s[k], s[k], RND);
    } else {
        _fk_(&t[k], s, u, k, D1);
        _fk_(&s[k], t, t, k, trig ? D2 : D_2);
    };
    return (pair){.a = &t[k], .b = &s[k]};
}

real *t_ln (series u, series e, int k) {
    CHECK(mpfr_sgn(e[0]) > 0); CHECK(u != e);
    if (!k) mpfr_log(u[k], e[k], RND); else _uk_(&u[k], e, u, k, &e[k], D1);
    return &u[k];
}

pair t_asin (series u, series c, series s, int k, bool trig) {
    CHECK(trig ? mpfr_cmpabs_ui(s[0], 1) <= 0 : true); CHECK(u != c && u != s && c != s);
    if (!k) {
        trig ? mpfr_asin(u[k], s[k], RND) : mpfr_asinh(u[k], s[k], RND);
        trig ? mpfr_cos(c[k], u[k], RND) : mpfr_cosh(c[k], u[k], RND);
    } else {
        _uk_(&u[k], c, u, k, &s[k], D1);
        _fk_(&c[k], s, u, k, trig ? D_1 : D1);
    };
    return (pair){.a = &u[k], .b = &c[k]};
}

pair t_acos (series u, series s, series c, int k, bool trig) {
    CHECK(trig ? mpfr_cmpabs_ui(c[0], 1) <= 0 : mpfr_cmp_si(c[0], 1) >= 0); CHECK(u != s && u != c && s != c);
    if (!k) {
        trig ? mpfr_acos(u[k], c[k], RND) : mpfr_acosh(u[k], c[k], RND);
        trig ? mpfr_sin(s[k], u[k], RND) : mpfr_sinh(s[k], u[k], RND);
        if (trig) mpfr_neg(s[k], s[k], RND);
    } else {
        _uk_(&u[k], s, u, k, &c[k], trig ? D_1 : D1);
        _fk_(&s[k], c, u, k, D1);
    };
    return (pair){.a = &u[k], .b = &s[k]};
}

pair t_atan (series u, series s, series t, int k, bool trig) {
    CHECK(trig ? true : mpfr_cmpabs_ui(t[0], 1) <= 0); CHECK(u != s && u != t && s != t);
    if (!k) {
        trig ? mpfr_atan(u[k], t[k], RND) : mpfr_atanh(u[k], t[k], RND);
        trig ? mpfr_sec(s[k], u[k], RND) : mpfr_sech(s[k], u[k], RND);
        mpfr_sqr(s[k], s[k], RND);
    } else {
        _uk_(&u[k], s, u, k, &t[k], D1);
        _fk_(&s[k], t, t, k, trig ? D2 : D_2);
    };
    return (pair){.a = &u[k], .b = &s[k]};
}

real *t_pwr (series p, series u, real a, int k) {
    CHECK(mpfr_sgn(u[0]) > 0); CHECK(p != u);
    if (!k) mpfr_pow(p[k], u[k], a, RND); else _uk_(&p[k], u, p, k, _fk_(&_p, p, u, k, a), D1);
    return &p[k];
}

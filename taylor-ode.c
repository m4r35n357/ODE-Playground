/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>
#include <mpfr.h>
#include "taylor-ode.h"

static mpfr_t __, _a, _m, _s, _p, D1, D_1, D2, D_2;

static char format[60];

void tsm_init(int dp) {
    sprintf(format, "%%+.%uRNe %%+.%uRNe %%+.%uRNe %%+.9RNe %%.3f\n", dp, dp, dp);
    mpfr_init_set_si(D1, 1, RND);
    mpfr_init_set_si(D2, 2, RND);
    mpfr_init_set_si(D_1, -1, RND);
    mpfr_init_set_si(D_2, -2, RND);
}

void t_out (mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t h, int step, clock_t since) {
    mpfr_mul_si(__, h, step, RND);
    mpfr_printf(format, x, y, z, __, (double)(clock() - since) / CLOCKS_PER_SEC);
}

void t_params (char **argv, int argc, ...) {
    va_list _;
    va_start(_, argc);
    for (int i = 9; i < argc; i++) {
        mpfr_init_set_str(*va_arg(_, mpfr_t *), argv[i], BASE, RND);
    }
    va_end(_);
}

series t_jet (int k) {
    CHECK(k > 0);
    series _ = malloc((size_t)k * sizeof (mpfr_t)); CHECK(_);
    for (int i = 0; i < k; i++) {
        mpfr_init(_[i]);
        mpfr_set_zero(_[i], 1);
    }
    return _;
}

series t_const (int k, mpfr_t a) {
    series _ = t_jet(k);
    mpfr_set(_[0], a, RND);
    return _;
}

mpfr_t *t_horner (series s, int n, mpfr_t h) {
    mpfr_set_zero(__, 1);
    for (int i = n; i >= 0; i--) {
        mpfr_fma(__, __, h, s[i], RND);
    }
    CHECK(mpfr_number_p(__) != 0);
    return &__;
}

void tsm (int n, mpfr_t h, int steps, series3 *j, parameters *p, clock_t t0) {
    triplet *v_k = malloc(sizeof (triplet)); CHECK(v_k);
    mpfr_inits(v_k->x, v_k->y, v_k->z, __, _a, _m, _s, _p, NULL);
    for (int step = 0; step < steps; step++) {
        for (int k = 0; k < n; k++) {
            ode(v_k, j->x, j->y, j->z, p, k);
            mpfr_div_si(j->x[k + 1], v_k->x, k + 1, RND);
            mpfr_div_si(j->y[k + 1], v_k->y, k + 1, RND);
            mpfr_div_si(j->z[k + 1], v_k->z, k + 1, RND);
        }
        t_out(j->x[0], j->y[0], j->z[0], h, step, t0);
        mpfr_swap(j->x[0], *t_horner(j->x, n, h));
        mpfr_swap(j->y[0], *t_horner(j->y, n, h));
        mpfr_swap(j->z[0], *t_horner(j->z, n, h));
    }
    t_out(j->x[0], j->y[0], j->z[0], h, steps, t0);
}

static int _half_ (int k) {
    return 1 + (k - (k % 2 ? 1 : 2)) / 2;
}

static mpfr_t *_prod_ (mpfr_t *_, series b, series a, int k, int k0, int k1) {
    mpfr_set_zero(*_, 1);
    for (int j = k0; j < k1; j++) {
        mpfr_fma(*_, b[j], a[k - j], *_, RND);
    }
    return _;
}

static mpfr_t *_fwd_ (mpfr_t *_, series b, series a, int k, mpfr_t *da_dt, mpfr_t scale) {
    mpfr_set_zero(*_, 1);
    for (int j = 0; j < k; j++) {
        mpfr_mul_si(*da_dt, a[k - j], k - j, RND);
        mpfr_fma(*_, b[j], *da_dt, *_, RND);
    }
    mpfr_div_si(*_, *_, k, RND);
    mpfr_mul(*_, *_, scale, RND);
    return _;
}

static mpfr_t *_rev_ (mpfr_t *_, series a, series b, mpfr_t *c_k, int k, mpfr_t *da_dt, bool neg) {
    mpfr_set_zero(*_, 1);
    for (int j = 1; j < k; j++) {
        mpfr_mul_si(*da_dt, a[k - j], k - j, RND);
        mpfr_fma(*_, b[j], *da_dt, *_, RND);
    }
    mpfr_div_si(*_, *_, k, RND);
    neg ? mpfr_add(*_, *c_k, *_, RND) : mpfr_sub(*_, *c_k, *_, RND);
    mpfr_div(*_, *_, b[0], RND);
    return _;
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
    if (!k) {
        u ? mpfr_set(q[k], u[k], RND) : mpfr_set_si(q[k], 1, RND);
    } else {
        _prod_(&q[k], q, v, k, 0, k);
        u ? mpfr_sub(q[k], u[k], q[k], RND) : mpfr_neg(q[k], q[k], RND);
    }
    mpfr_div(q[k], q[k], v[0], RND);
    return &q[k];
}

mpfr_t *t_sqrt (series r, series u, int k) {
    CHECK(mpfr_sgn(u[0]) > 0); CHECK(r != u);
    if (!k) {
        mpfr_sqrt(r[k], u[k], RND);
    } else {
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
    if (!k) {
        mpfr_exp(e[k], u[k], RND);
    } else {
        _fwd_(&e[k], e, u, k, &__, D1);
    };
    return &e[k];
}

pair t_sin_cos (series s, series c, series u, int k, bool trig) {
    CHECK(s != c && s != u && c != u);
    if (!k) {
        trig ? mpfr_sin_cos(s[k], c[k], u[k], RND) : mpfr_sinh_cosh(s[k], c[k], u[k], RND);
    } else {
        _fwd_(&s[k], c, u, k, &__, D1);
        _fwd_(&c[k], s, u, k, &__, trig ? D_1 : D1);
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
        _fwd_(&t[k], s, u, k, &__, D1);
        _fwd_(&s[k], t, t, k, &__, trig ? D2 : D_2);
    };
    return (pair){.a = &t[k], .b = &s[k]};
}

mpfr_t *t_ln (series l, series u, int k) {
    CHECK(mpfr_sgn(u[0]) > 0); CHECK(l != u);
    if (!k) {
        mpfr_log(l[k], u[k], RND);
    } else {
        _rev_(&l[k], l, u, &u[k], k, &__, false);
    };
    return &l[k];
}

pair t_asin (series u, series g, series s, int k, bool trig) {
    CHECK(trig ? mpfr_cmpabs_ui(s[0], 1) <= 0 : true); CHECK(u != g && u != s && g != s);
    if (!k) {
        trig ? mpfr_asin(u[k], s[k], RND) : mpfr_asinh(u[k], s[k], RND);
        trig ? mpfr_cos(g[k], u[k], RND) : mpfr_cosh(g[k], u[k], RND);
    } else {
        _rev_(&u[k], u, g, &s[k], k, &__, false);
        _fwd_(&g[k], s, u, k, &__, trig ? D_1 : D1);
    };
    return (pair){.a = &u[k], .b = &g[k]};
}

pair t_acos (series u, series g, series c, int k, bool trig) {
    CHECK(trig ? mpfr_cmpabs_ui(c[0], 1) <= 0 : mpfr_cmp_si(c[0], 1) >= 0); CHECK(u != g && u != c && g != c);
    if (!k) {
        trig ? mpfr_acos(u[k], c[k], RND) : mpfr_acosh(u[k], c[k], RND);
        trig ? mpfr_sin(g[k], u[k], RND) : mpfr_sinh(g[k], u[k], RND);
        if (trig) mpfr_neg(g[k], g[k], RND);
    } else {
        _rev_(&u[k], u, g, &c[k], k, &__, trig);
        _fwd_(&g[k], c, u, k, &__, D1);
    };
    return (pair){.a = &u[k], .b = &g[k]};
}

pair t_atan (series u, series g, series t, int k, bool trig) {
    CHECK(trig ? true : mpfr_cmpabs_ui(t[0], 1) <= 0); CHECK(u != g && u != t && g != t);
    if (!k) {
        trig ? mpfr_atan(u[k], t[k], RND) : mpfr_atanh(u[k], t[k], RND);
        trig ? mpfr_sec(g[k], u[k], RND) : mpfr_sech(g[k], u[k], RND);
        mpfr_sqr(g[k], g[k], RND);
    } else {
        _rev_(&u[k], u, g, &t[k], k, &__, false);
        _fwd_(&g[k], t, t, k, &__, trig ? D2 : D_2);
    };
    return (pair){.a = &u[k], .b = &g[k]};
}

mpfr_t *t_pwr (series p, series u, mpfr_t a, int k) {
    CHECK(mpfr_sgn(u[0]) > 0); CHECK(p != u);
    if (!k) {
        mpfr_pow(p[k], u[k], a, RND);
    } else {
        _rev_(&p[k], p, u, _fwd_(&_p, p, u, k, &__, a), k, &__, false);
    }
    return &p[k];
}

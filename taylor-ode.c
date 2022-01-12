/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

const int BASE = 10;

const mpfr_rnd_t RND = MPFR_RNDN;

static mpfr_t D1, D2, D_1, D_2, _, __, _a, _b, _const, _abs, _mul, _sqr;

static char f[60];

void t_init (int p) {
    p == 0 ? sprintf(f, "%%.RNe %%.RNe %%.RNe %%.9RNe\n") : sprintf(f, "%%+.%uRNe %%+.%uRNe %%+.%uRNe %%+.9RNe\n", p, p, p);
    mpfr_init_set_si(D1, 1, RND);
    mpfr_init_set_si(D2, 2, RND);
    mpfr_init_set_si(D_1, -1, RND);
    mpfr_init_set_si(D_2, -2, RND);
    mpfr_inits( _, __, _a, _b, _const, _abs, _mul, _sqr, NULL);
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
    mpfr_t *s = calloc((size_t)n + 1, sizeof (mpfr_t));
    if (s == NULL) {
        fprintf(stderr, "Allocation failure!\n");
        exit(1);
    }
    for (int i = 0; i <= n; i++) {
        mpfr_init(s[i]);
    }
    return s;
}

void t_output (mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t h, int step) {
    mpfr_mul_si(_, h, step, RND);
    mpfr_printf(f, x, y, z, _);
}

mpfr_t *t_horner (series s, int n, mpfr_t h) {
    mpfr_set_zero(_, 1);
    for (int i = n; i >= 0; i--) {
        mpfr_fma(_, _, h, s[i], RND);
    }
    if (mpfr_number_p(_) == 0) {
        fprintf(stderr, "Value error!\n");
        exit(2);
    }
    return &_;
}

void tsm (int argc, char **argv, int n, mpfr_t h, int steps, mpfr_t x0, mpfr_t y0, mpfr_t z0) {
    series x = t_jet(n); mpfr_set(x[0], x0, RND);
    series y = t_jet(n); mpfr_set(y[0], y0, RND);
    series z = t_jet(n); mpfr_set(z[0], z0, RND);
    void *p = get_p(argc, argv, n);
    components *vk = malloc(sizeof (components));
    mpfr_inits(vk->x, vk->y, vk->z, NULL);
    for (int step = 0; step < steps; step++) {
        for (int k = 0; k < n; k++) {
            ode(vk, x, y, z, p, k);
            mpfr_div_si(x[k + 1], vk->x, k + 1, RND);
            mpfr_div_si(y[k + 1], vk->y, k + 1, RND);
            mpfr_div_si(z[k + 1], vk->z, k + 1, RND);
        }
        t_output(x[0], y[0], z[0], h, step);
        mpfr_swap(x[0], *t_horner(x, n, h));
        mpfr_swap(y[0], *t_horner(y, n, h));
        mpfr_swap(z[0], *t_horner(z, n, h));
    }
    t_output(x[0], y[0], z[0], h, steps);
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

mpfr_t *t_mul (series u, series v, int k) {
    mpfr_set_zero(_mul, 1);
    for (int j = 0; j <= k; j++) {
        mpfr_fma(_mul, u[j], v[k - j], _mul, RND);
    }
    return &_mul;
}

mpfr_t *t_sqr (series u, int k) {
    mpfr_set_zero(_, 1);
    for (int j = 0; j <= (k - (k % 2 ? 1 : 2)) / 2; j++) {
        mpfr_fma(_, u[j], u[k - j], _, RND);
    }
    mpfr_mul_2si(_sqr, _, 1, RND);
    if (!(k % 2)) mpfr_fma(_sqr, u[k / 2], u[k / 2], _sqr, RND);
    return &_sqr;
}

mpfr_t *t_div (series q, series u, series v, int k) {
    assert(mpfr_zero_p(v[0]) == 0);
    assert(q != u && q != v);
    if (k == 0) {
        mpfr_set(_, u == NULL ? D1 : u[0], RND);
    } else {
        mpfr_set_zero(_, 1);
        for (int j = 0; j < k; j++) {
            mpfr_fma(_, q[j], v[k - j], _, RND);
        }
        u == NULL ? mpfr_neg(_, _, RND) : mpfr_sub(_, u[k], _, RND);
    }
    mpfr_div(q[k], _, v[0], RND);
    return &q[k];
}

mpfr_t *t_sqrt (series r, series u, int k) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(r != u);
    if (k == 0) {
        mpfr_sqrt(r[0], u[0], RND);
    } else {
        mpfr_set_zero(_, 1);
        for (int j = 1; j <= (k - (k % 2 ? 1 : 2)) / 2; j++) {
            mpfr_fma(_, r[j], r[k - j], _, RND);
        }
        mpfr_mul_2si(_, _, 1, RND);
        if (!(k % 2)) mpfr_fma(_, r[k / 2], r[k / 2], _, RND);
        mpfr_sub(_, u[k], _, RND);
        mpfr_div_2si(_, _, 1, RND);
        mpfr_div(r[k], _, r[0], RND);
    }
    return &r[k];
}

mpfr_t *t_exp (series e, series u, int k) {
    assert(e != u);
    if (k == 0) {
        mpfr_exp(e[0], u[0], RND);
    } else {
        mpfr_set_zero(_, 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_si(__, u[k - j], k - j, RND);
            mpfr_fma(_, e[j], __, _, RND);
        }
        mpfr_div_si(e[k], _, k, RND);
    }
    return &e[k];
}

pair t_sin_cos (series s, series c, series u, int k, geometry g) {
    assert(s != c && s != u && c != u);
    if (k == 0) {
        g == TRIG ? mpfr_sin_cos(s[0], c[0], u[0], RND) : mpfr_sinh_cosh(s[0], c[0], u[0], RND);
    } else {
        mpfr_set_zero(_a, 1);
        mpfr_set_zero(_b, 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_si(__, u[k - j], k - j, RND);
            mpfr_fma(_b, c[j], __, _b, RND);
            mpfr_fma(_a, s[j], __, _a, RND);
        }
        mpfr_div_si(s[k], _b, k, RND);
        mpfr_div_si(c[k], _a, k, RND);
        if (g == TRIG) mpfr_neg(c[k], c[k], RND);
    }
    return (pair) {&s[k], &c[k]};
}

pair t_tan_sec2 (series t, series s, series u, int k, geometry g) {
    assert(t != s && t != u && s != u);
    if (k == 0) {
        g == TRIG ? mpfr_tan(t[0], u[0], RND) : mpfr_tanh(t[0], u[0], RND);
        mpfr_sqr(s[0], t[0], RND);
        g == TRIG ? mpfr_add(s[0], D1, s[0], RND) : mpfr_sub(s[0], D1, s[0], RND);
    } else {
        mpfr_set_zero(_, 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_si(__, u[k - j], k - j, RND);
            mpfr_fma(_, s[j], __, _, RND);
        }
        mpfr_div_si(t[k], _, k, RND);
        mpfr_set_zero(_, 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_si(__, t[k - j], k - j, RND);
            mpfr_fma(_, t[j], __, _, RND);
        }
        mpfr_div_si(s[k], _, k, RND);
        mpfr_mul(s[k], s[k], g == TRIG ? D2 : D_2, RND);
    }
    return (pair) {&t[k], &s[k]};
}

mpfr_t *t_pwr (series p, series u, mpfr_t a, int k) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(p != u);
    if (k == 0) {
        mpfr_pow(p[0], u[0], a, RND);
    } else {
        mpfr_set_zero(_, 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_si(__, a, k - j, RND);
            mpfr_sub_si(__, __, j, RND);
            mpfr_mul(__, __, u[k - j], RND);
            mpfr_fma(_, p[j], __, _, RND);
        }
        mpfr_div_si(_, _, k, RND);
        mpfr_div(p[k], _, u[0], RND);
    }
    return &p[k];
}

mpfr_t *t_ln (series l, series u, int k) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(l != u);
    if (k == 0) {
        mpfr_log(l[0], u[0], RND);
    } else {
        mpfr_set_zero(_, 1);
        for (int j = 1; j < k; j++) {
            mpfr_mul_si(__, l[j], j, RND);
            mpfr_fma(_, u[k - j], __, _, RND);
        }
        mpfr_div_si(_, _, k, RND);
        mpfr_sub(_, u[k], _, RND);
        mpfr_div(l[k], _, u[0], RND);
    }
    return &l[k];
}

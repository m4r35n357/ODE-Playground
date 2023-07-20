/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <time.h>
#include <mpfr.h>
#include "taylor-ode.h"

const int BASE = 10;

const mpfr_rnd_t RND = MPFR_RNDN;

static mpfr_t _, __, _a, _b, _const, _abs, _mul, _sqr;

static char f[60];

void t_init (int places) {
    sprintf(f, "%%+.%uRNe %%+.%uRNe %%+.%uRNe %%+.9RNe %%.3f\n", places, places, places);
    mpfr_inits( _, __, _a, _b, _const, _abs, _mul, _sqr, NULL);
}

void t_params (char **argv, int argc, ...) {
    fprintf(stderr, "[ ");
	for (int i = 0; i < argc; i++) {
        fprintf(stderr, "%s ", argv[i]);
    }
    fprintf(stderr, "]\n");
    va_list model;
    va_start(model, argc);
    for (int i = 9; i < argc; i++) {
        mpfr_init_set_str(*va_arg(model, mpfr_t *), argv[i], BASE, RND);
    }
    va_end(model);
}

series t_jet (int n) {
    mpfr_t *s = calloc((size_t)n, sizeof (mpfr_t));
    if (!s) {
        fprintf(stderr, "Allocation failure!\n");
        exit(1);
    }
    for (int i = 0; i <= n; i++) {
        mpfr_init(s[i]);
    }
    return s;
}

void t_out (mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t h, int step, clock_t since) {
    mpfr_mul_si(_, h, step, RND);
    mpfr_printf(f, x, y, z, _, (double)(clock() - since) / CLOCKS_PER_SEC);
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

void tsm (int n, mpfr_t h, int steps, mpfr_t x0, mpfr_t y0, mpfr_t z0, void *p, clock_t t0) {
    series x = t_jet(n + 1); mpfr_set(x[0], x0, RND);
    series y = t_jet(n + 1); mpfr_set(y[0], y0, RND);
    series z = t_jet(n + 1); mpfr_set(z[0], z0, RND);
    components *vk = malloc(sizeof (components));
    mpfr_inits(vk->x, vk->y, vk->z, NULL);
    for (int step = 0; step < steps; step++) {
        for (int k = 0; k < n; k++) {
            ode(vk, x, y, z, p, k);
            mpfr_div_si(x[k + 1], vk->x, k + 1, RND);
            mpfr_div_si(y[k + 1], vk->y, k + 1, RND);
            mpfr_div_si(z[k + 1], vk->z, k + 1, RND);
        }
        t_out(x[0], y[0], z[0], h, step, t0);
        mpfr_swap(x[0], *t_horner(x, n, h));
        mpfr_swap(y[0], *t_horner(y, n, h));
        mpfr_swap(z[0], *t_horner(z, n, h));
    }
    t_out(x[0], y[0], z[0], h, steps, t0);
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
    assert(mpfr_zero_p(u[0]) == 0);
    if (mpfr_sgn(u[0]) < 0) {
        mpfr_neg(_abs, u[k], RND);
    } else {
        mpfr_set(_abs, u[k], RND);
    }
    return &_abs;
}

static mpfr_t *_prod (series u, series v, int k) {
    mpfr_set_zero(_mul, 1);
    for (int j = 0; j <= k; j++) {
        mpfr_fma(_mul, u[j], v[k - j], _mul, RND);
    }
    return &_mul;
}

mpfr_t *t_mul (series u, series v, int k) {
    return _prod(u, v, k);
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
        u ? mpfr_set(_, u[0], RND) : mpfr_set_si(_, 1, RND);
    } else {
        mpfr_set_zero(_, 1);
        for (int j = 0; j < k; j++) {
            mpfr_fma(_, q[j], v[k - j], _, RND);
        }
        u ? mpfr_sub(_, u[k], _, RND) : mpfr_neg(_, _, RND);
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

static mpfr_t *_exp (series e, series u, int k) {
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

pair t_sin_cos (series s, series c, series u, int k, geometry G) {
    assert(s != c && s != u && c != u);
    if (k == 0) {
        G == TRIG ? mpfr_sin_cos(s[0], c[0], u[0], RND) : mpfr_sinh_cosh(s[0], c[0], u[0], RND);
    } else {
        mpfr_set(s[k], *_exp(c, u, k), RND);
        mpfr_set(c[k], *_exp(s, u, k), RND);
        if (G == TRIG) mpfr_neg(c[k], c[k], RND);
    }
    return (pair){&s[k], &c[k]};
}

pair t_tan_sec2 (series t, series s, series u, int k, geometry G) {
    assert(t != s && t != u && s != u);
    if (k == 0) {
        G == TRIG ? mpfr_tan(t[0], u[0], RND) : mpfr_tanh(t[0], u[0], RND);
        mpfr_sqr(s[0], t[0], RND);
        G == TRIG ? mpfr_add_si(s[0], s[0], 1, RND) : mpfr_si_sub(s[0], 1, s[0], RND);
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
        mpfr_mul_si(s[k], s[k], G == TRIG ? 2 : -2, RND);
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

mpfr_t *t_ipwr (series p, series u, int a, int k) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(p != u);
    if (k == 0) {
        mpfr_pow_si(p[0], u[0], a, RND);
    } else {
        mpfr_set_zero(_, 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_si(__, u[k - j], a * (k - j) - j, RND);
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

pair t_asin (series a, series g, series u, int k, geometry G) {
    assert(G == TRIG ? mpfr_cmpabs_ui(u[0], 1) <= 0 : 1);
    assert(a != g && a != u && g != u);
    if (k == 0) {
        G == TRIG ? mpfr_asin(a[0], u[0], RND) : mpfr_asinh(a[0], u[0], RND);
        mpfr_sqr(g[0], u[0], RND);
        G == TRIG ? mpfr_si_sub(g[0], 1, g[0], RND) : mpfr_add_si(g[0], g[0], 1, RND);
        mpfr_sqrt(g[0], g[0], RND);
    } else {
        mpfr_set_zero(_, 1);
        for (int j = 1; j < k; j++) {
            mpfr_mul_si(__, a[j], j, RND);
            mpfr_fma(_, g[k - j], __, _, RND);
        }
        mpfr_div_si(_, _, k, RND);
        mpfr_sub(_, u[k], _, RND);
        mpfr_div(a[k], _, g[0], RND);
        mpfr_set_zero(_, 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_si(__, a[k - j], k - j, RND);
            mpfr_fma(_, u[j], __, _, RND);
        }
        mpfr_div_si(g[k], _, G == TRIG ? - k : k, RND);
    }
    return (pair) {&a[k], &g[k]};
}

pair t_acos (series a, series g, series u, int k, geometry G) {
    assert(G == TRIG ? mpfr_cmpabs_ui(u[0], 1) <= 0 : mpfr_cmp_si(u[0], 1) >= 0);
    assert(a != g && a != u && g != u);
    if (k == 0) {
        G == TRIG ? mpfr_acos(a[0], u[0], RND) : mpfr_acosh(a[0], u[0], RND);
        mpfr_sqr(g[0], u[0], RND);
        G == TRIG ? mpfr_si_sub(g[0], 1, g[0], RND) : mpfr_sub_si(g[0], g[0], 1, RND);
        mpfr_sqrt(g[0], g[0], RND);
        if (G == TRIG) mpfr_neg(g[0], g[0], RND);
    } else {
        mpfr_set_zero(_, 1);
        for (int j = 1; j < k; j++) {
            mpfr_mul_si(__, a[j], j, RND);
            mpfr_fma(_, g[k - j], __, _, RND);
        }
        mpfr_div_si(_, _, G == TRIG ? - k : k, RND);
        mpfr_sub(_, u[k], _, RND);
        mpfr_div(a[k], _, g[0], RND);
        mpfr_set_zero(_, 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_si(__, a[k - j], k - j, RND);
            mpfr_fma(_, u[j], __, _, RND);
        }
        mpfr_div_si(g[k], _, k, RND);
    }
    return (pair) {&a[k], &g[k]};
}

pair t_atan (series a, series g, series u, int k, geometry G) {
    assert(G == TRIG ? 1 : mpfr_cmpabs_ui(u[0], 1) <= 0);
    assert(a != g && a != u && g != u);
    if (k == 0) {
        G == TRIG ? mpfr_atan(a[0], u[0], RND) : mpfr_atanh(a[0], u[0], RND);
        mpfr_sqr(g[0], u[0], RND);
        G == TRIG ? mpfr_add_si(g[0], g[0], 1, RND) : mpfr_si_sub(g[0], 1, g[0], RND);
    } else {
        mpfr_set_zero(_, 1);
        for (int j = 1; j < k; j++) {
            mpfr_mul_si(__, a[j], j, RND);
            mpfr_fma(_, g[k - j], __, _, RND);
        }
        mpfr_div_si(_, _, k, RND);
        mpfr_sub(_, u[k], _, RND);
        mpfr_div(a[k], _, g[0], RND);
        mpfr_set_zero(_, 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_si(__, u[k - j], k - j, RND);
            mpfr_fma(_, u[j], __, _, RND);
        }
        mpfr_div_si(_, _, G == TRIG ? k : - k, RND);
        mpfr_mul_2si(g[k], _, 1, RND);
    }
    return (pair) {&a[k], &g[k]};
}

/*
 * Arbitrary-Order Taylor Series Integrator, together with a well-tested set of recurrence relations
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>
#include "taylor-ode.h"

static real __, _a, _m, _s, D0, D1, PI_2;

static char format[60];

void tsm_init (int dp) {
    sprintf(format, "%%+.%uRNe %%+.%uRNe %%+.%uRNe %%+.9RNe %%.3f\n", dp, dp, dp);
    mpfr_inits(__, _a, _m, _s, PI_2, NULL);
    mpfr_init_set_si(D0, 0, RND);
    mpfr_init_set_si(D1, 1, RND);
    mpfr_const_pi(PI_2, RND);
    mpfr_div_2si(PI_2, PI_2, 1, RND);
}

void tsm_get_p (char **argv, int argc, ...) {
    va_list _;
    va_start(_, argc);
    for (int i = 9; i < argc; i++) mpfr_init_set_str(*va_arg(_, real *), argv[i], BASE, RND);
    va_end(_);
}

series tsm_jet (int n) {
    CHECK(n > 0);
    series _ = malloc((size_t)n * sizeof (real)); CHECK(_);
    for (int i = 0; i < n; i++) {
        mpfr_init(_[i]);
        mpfr_set_zero(_[i], 1);
    }
    return _;
}

real *horner (series u, int o, real h) {
    mpfr_set_zero(__, 1);
    for (int i = o; i >= 0; i--) mpfr_fma(__, __, h, u[i], RND);
    CHECK(mpfr_number_p(__) != 0);
    return &__;
}

void _out_ (real x, real y, real z, real h, int step, clock_t since) {
    mpfr_mul_si(__, h, step, RND);
    mpfr_printf(format, x, y, z, __, (double)(clock() - since) / CLOCKS_PER_SEC);
}

void tsm (int o, real h, int steps, triplet *v, series x, series y, series z, const model *p, clock_t t0) {
    for (int step = 0; step < steps; step++) {
        for (int k = 0; k < o; k++) {
            ode(v, x, y, z, p, k);
            mpfr_div_si(x[k + 1], v->x, k + 1, RND);
            mpfr_div_si(y[k + 1], v->y, k + 1, RND);
            mpfr_div_si(z[k + 1], v->z, k + 1, RND);
        }
        _out_(x[0], y[0], z[0], h, step, t0);
        mpfr_swap(x[0], *horner(x, o, h));
        mpfr_swap(y[0], *horner(y, o, h));
        mpfr_swap(z[0], *horner(z, o, h));
    }
    _out_(x[0], y[0], z[0], h, steps, t0);
}

real *t_abs (const series u, int k) {
    if (!k) CHECK(mpfr_zero_p(u[0]) == 0);
    mpfr_sgn(u[0]) < 0 ? mpfr_neg(_a, u[k], RND) : mpfr_set(_a, u[k], RND);
    return &_a;
}

static real *_cauchy_ (real *_, const series b, const series a, int k, int k0, int k1) {
    mpfr_set_zero(*_, 1);
    for (int j = k0; j <= k1; j++) mpfr_fma(*_, b[j], a[k - j], *_, RND);
    return _;
}

real *t_mul (const series u, const series v, int k) {
    return _cauchy_(&_m, u, v, k, 0, k);
}

real *t_div (series q, const series u, const series v, int k) {
    if (k) {
        mpfr_sub(q[k], u ? u[k] : D0, *_cauchy_(q + k, q, v, k, 0, k - 1), RND);
        mpfr_div(q[k], q[k], v[0], RND);
    } else {
        CHECK(q != u && q != v);
        CHECK(mpfr_zero_p(v[0]) == 0);
        mpfr_div(q[k], u ? u[k] : D1, v[0], RND);
    }
    return q + k;
}

static real *_half_ (real *_, const series a, int k, int k0, bool even) {
    mpfr_mul_2si(*_, *_cauchy_(_, a, a, k, k0, (k - (even ? 1 : 2)) / 2), 1, RND);
    if (!even) mpfr_fma(*_, a[k / 2], a[k / 2], *_, RND);
    return _;
}

real *t_sqr (const series u, int k) {
    return _half_(&_s, u, k, 0, k % 2);
}

real *t_sqrt (series r, const series u, int k) {
    if (k) {
        mpfr_sub(r[k], u[k], *_half_(r + k, r, k, 1, k % 2), RND);
        mpfr_div_2si(r[k], r[k], 1, RND);
        mpfr_div(r[k], r[k], r[0], RND);
    } else {
        CHECK(r != u);
        CHECK(mpfr_sgn(u[0]) > 0);
        mpfr_sqrt(r[k], u[k], RND);
    }
    return r + k;
}

real *t_pwr (series p, const series u, real a, int k) {
    if (k) {
        mpfr_set_zero(p[k], 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_si(__, a, k - j, RND);
            mpfr_sub_si(__, __, j, RND);
            mpfr_mul(__, __, u[k - j], RND);
            mpfr_fma(p[k], p[j], __, p[k], RND);
        }
        mpfr_div_si(p[k], p[k], k, RND);
        mpfr_div(p[k], p[k], u[0], RND);
    } else {
        CHECK(p != u);
        CHECK(mpfr_sgn(u[0]) > 0);
        mpfr_pow(p[k], u[k], a, RND);
    }
    return p + k;
}

static real *_chain_ (real *_, const series dfdu, const series u, int k, const series fk, int scale) {
    mpfr_set_zero(*_, 1);
    for (int j = fk ? 1 : 0; j < k; j++) {
        mpfr_mul_si(__, u[k - j], k - j, RND);
        mpfr_fma(*_, dfdu[j], __, *_, RND);
    }
    mpfr_div_si(*_, *_, k, RND);
    mpfr_mul_si(*_, *_, scale, RND);
    if (fk) {
        mpfr_sub(*_, *fk, *_, RND);
        mpfr_div(*_, *_, dfdu[0], RND);
    }
    return _;  // f[k] if fk NULL (forward), u[k] if non-NULL (reverse)
}

real *t_exp (series e, const series u, int k) {
    if (k) {
        _chain_(e + k, e, u, k, NULL, 1);
    } else {
        CHECK(e != u);
        mpfr_exp(e[k], u[k], RND);
    }
    return e + k;
}

real *t_ln (series u, const series e, int k) {
    if (k) {
        _chain_(u + k, e, u, k, e + k, 1);
    } else {
        CHECK(u != e);
        CHECK(mpfr_sgn(e[0]) > 0);
        mpfr_log(u[k], e[k], RND);
    }
    return u + k;
}

pair t_sin_cos (series s, series c, const series u, int k, bool trig) {
    if (k) {
        mpfr_set_zero(s[k], 1);
        mpfr_set_zero(c[k], 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_si(__, u[k - j], k - j, RND);
            mpfr_fma(s[k], c[j], __, s[k], RND);
            mpfr_fma(c[k], s[j], __, c[k], RND);
        }
        mpfr_div_si(s[k], s[k], k, RND);
        mpfr_div_si(c[k], c[k], trig ? -k : k, RND);
    } else {
        CHECK(s != c && s != u && c != u);
        trig ? mpfr_sin_cos(s[k], c[k], u[k], RND) : mpfr_sinh_cosh(s[k], c[k], u[k], RND);
    }
    return (pair){ s + k, c + k };
}

pair t_tan_sec2 (series t, series s, const series u, int k, bool trig) {
    if (k) {
        _chain_(t + k, s, u, k, NULL, 1);
        _chain_(s + k, t, t, k, NULL, trig ? 2 : -2);
    } else {
        CHECK(t != s && t != u && s != u);
        CHECK(trig ? mpfr_cmpabs(u[0], PI_2) < 0 : true);
        trig ? mpfr_tan(t[k], u[k], RND) : mpfr_tanh(t[k], u[k], RND);
        trig ? mpfr_sec(s[k], u[k], RND) : mpfr_sech(s[k], u[k], RND);
        mpfr_sqr(s[k], s[k], RND);
    }
    return (pair){ t + k, s + k };
}

pair t_asin_cos (series u, series c, const series s, int k, bool trig) {
    if (k) {
        _chain_(u + k, c, u, k, s + k, 1);
        _chain_(c + k, s, u, k, NULL, trig ? -1 : 1);
    } else {
        CHECK(u != c && u != s && c != s);
        CHECK(trig ? mpfr_cmpabs_ui(s[0], 1) < 0 : true);
        trig ? mpfr_asin(u[k], s[k], RND) : mpfr_asinh(u[k], s[k], RND);
        trig ?  mpfr_cos(c[k], u[k], RND) :  mpfr_cosh(c[k], u[k], RND);
    }
    return (pair){ u + k, c + k };
}

pair t_acos_sin (series u, series s, const series c, int k, bool trig) {
    if (k) {
        _chain_(u + k, s, u, k, c + k, trig ? -1 : 1);
        _chain_(s + k, c, u, k, NULL, 1);
    } else {
        CHECK(u != s && u != c && s != c);
        CHECK(trig ? mpfr_cmpabs_ui(c[0], 1) < 0 : mpfr_cmp_si(c[0], 1) > 0);
        trig ? mpfr_acos(u[k], c[k], RND) : mpfr_acosh(u[k], c[k], RND);
        trig ?  mpfr_sin(s[k], u[k], RND) :  mpfr_sinh(s[k], u[k], RND);
        if (trig) mpfr_neg(s[k], s[k], RND);
    }
    return (pair){ u + k, s + k };
}

pair t_atan_sec2 (series u, series s, const series t, int k, bool trig) {
    if (k) {
        _chain_(u + k, s, u, k, t + k, 1);
        _chain_(s + k, t, t, k, NULL, trig ? 2 : -2);
    } else {
        CHECK(u != s && u != t && s != t);
        CHECK(trig ? true : mpfr_cmpabs_ui(t[0], 1) < 0);
        trig ? mpfr_atan(u[k], t[k], RND) : mpfr_atanh(u[k], t[k], RND);
        trig ?  mpfr_sec(s[k], u[k], RND) :  mpfr_sech(s[k], u[k], RND);
        mpfr_sqr(s[k], s[k], RND);
    }
    return (pair){ u + k, s + k };
}

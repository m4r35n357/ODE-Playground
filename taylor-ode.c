/*
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

const int BASE = 10;

const mpfr_rnd_t RND = MPFR_RNDN;

static long dp;

static mpfr_t D1, D_1, D2, D_2, _, __, ___;

void t_tempvars (void) {
    mpfr_inits(_, __, ___, NULL);
    mpfr_init_set_ui(D1, 1, RND);
    mpfr_init_set_si(D_1, -1, RND);
    mpfr_init_set_ui(D2, 2, RND);
    mpfr_init_set_si(D_2, -2, RND);
}

void t_output (mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t h, long step) {
    char s[42];
    dp == 0 ? sprintf(s, "%%.RNe %%.RNe %%.RNe %%.9RNe\n") : sprintf(s, "%%+.%luRNe %%+.%luRNe %%+.%luRNe %%+.9RNe\n", dp, dp, dp);
    mpfr_mul_ui(_, h, step, RND);
    mpfr_printf(s, x, y, z, _);
}

void t_stepper (char **argv, long *n, mpfr_t *h, long *nsteps) {
    dp = strtol(argv[1], NULL, BASE);
    mpfr_set_default_prec(strtod(argv[2], NULL) * 3.322);
    fprintf(stderr, " MPFR default precision: %lu bits\n", mpfr_get_default_prec());
    *n = strtol(argv[3], NULL, BASE);
    mpfr_init_set_str(*h, argv[4], BASE, RND);
    *nsteps = strtol(argv[5], NULL, BASE);
    t_tempvars();
}

void t_args (char **argv, int argc, ...) {
    va_list vars;
    va_start(vars, argc);
    for (int i = 6; i < argc; i++) {
        mpfr_init_set_str(*va_arg(vars, mpfr_t *), argv[i], BASE, RND);
    }
    va_end(vars);
}

series t_jet (int n) {
    assert(n > 0);
    mpfr_t *jet = malloc(sizeof (mpfr_t) * n);
    for (int i = 0; i < n; i++) {
        mpfr_init_set_ui(jet[i], 0, RND);
    }
    return (series){jet, n};
}

series t_jet_c (int n, mpfr_t value) {
    series s = t_jet(n);
    mpfr_set(s.a[0], value, RND);
    return s;
}

void t_next (series jet, mpfr_t dot, int k, sign s) {
    s == NEG ? mpfr_div_si(jet.a[k + 1], dot, - (k + 1), RND) : mpfr_div_ui(jet.a[k + 1], dot, k + 1, RND);
}

mpfr_t *t_horner (series jet, mpfr_t h) {
    mpfr_set_zero(_, 1);
    for (int i = jet.size - 1; i >= 0; i--) {
        mpfr_fma(_, _, h, jet.a[i], RND);
    }
    if (mpfr_number_p(_) == 0) exit(1);
    mpfr_swap(jet.a[0], _);
    return &jet.a[0];
}

mpfr_t *t_abs (series u, int k) {
    assert(k >= 0);
    mpfr_sgn(u.a[0]) < 0 ? mpfr_neg(_, u.a[k], RND) : mpfr_set(_, u.a[k], RND);
    return &_;
}

static mpfr_t *cauchy (mpfr_t *c, series a, series b, int k, int lower, int upper) {
    mpfr_set_zero(*c, 1);
    for (int j = lower; j <= upper; j++) {
        mpfr_fma(*c, a.a[j], b.a[k - j], *c, RND);
    }
    return c;
}

mpfr_t *t_prod (series u, series v, int k) {
    assert(k >= 0);
    return cauchy(&__, u, v, k, 0, k);
}

mpfr_t *t_quot (series q, series u, series v, int k) {
    assert(mpfr_zero_p(v.a[0]) == 0);
    assert(q.a != u.a && q.a != v.a);
    assert(k >= 0);
    mpfr_sub(_, u.a[k], *cauchy(&__, q, v, k, 0, k - 1), RND);
    mpfr_div(q.a[k], _, v.a[0], RND);
    return &q.a[k];
}

mpfr_t *t_sqr (series u, int k) {
    assert(k >= 0);
    mpfr_mul_2ui(_, *cauchy(&__, u, u, k, 0, (k - (k % 2 == 0 ? 2 : 1)) / 2), 1, RND);
    if (k % 2 == 0) mpfr_fma(_, u.a[k / 2], u.a[k / 2], _, RND);
    return &_;
}

mpfr_t *t_sqrt (series r, series u, int k) {
    assert(mpfr_sgn(u.a[0]) > 0);
    assert(r.a != u.a);
    assert(k >= 0);
    if (k == 0) {
        mpfr_sqrt(r.a[k], u.a[0], RND);
    } else {
        mpfr_mul_2ui(_, *cauchy(&__, r, r, k, 1, (k - (k % 2 == 0 ? 2 : 1)) / 2), 1, RND);
        if (k % 2 == 0) mpfr_fma(_, r.a[k / 2], r.a[k / 2], _, RND);
        mpfr_sub(_, u.a[k], _, RND);
        mpfr_div_2ui(_, _, 1, RND);
        mpfr_div(r.a[k], _, r.a[0], RND);
    }
    return &r.a[k];
}

static mpfr_t *d_cauchy (mpfr_t *f, series h, series u, int k, int lower, int upper, mpfr_t factor) {
    assert(f != &_ && h.a != &_ && u.a != &_);  // _ is used internally so it cannot be a parameter
    mpfr_set_zero(*f, 1);
    for (int j = lower; j <= upper; j++) {
        mpfr_mul_ui(_, u.a[k - j], k - j, RND);
        mpfr_fma(*f, h.a[j], _, *f, RND);
    }
    mpfr_div_ui(*f, *f, k, RND);
    mpfr_mul(*f, *f, factor, RND);
    return f;
}

mpfr_t *t_exp (series e, series u, int k) {
    assert(e.a != u.a);
    assert(k >= 0);
    if (k == 0) {
        mpfr_exp(e.a[k], u.a[0], RND);
    } else {
        d_cauchy(&e.a[k], e, u, k, 0, k - 1, D1);
    }
    return &e.a[k];
}

tuple t_sin_cos (series s, series c, series u, int k, geometry g) {
    assert(s.a != c.a && s.a != u.a && c.a != u.a);
    assert(k >= 0);
    if (k == 0) {
        g == TRIG ? mpfr_sin_cos(s.a[k], c.a[k], u.a[0], RND) : mpfr_sinh_cosh(s.a[k], c.a[k], u.a[0], RND);
    } else {
        d_cauchy(&s.a[k], c, u, k, 0, k - 1, D1);
        d_cauchy(&c.a[k], s, u, k, 0, k - 1, g == TRIG ? D_1 : D1);
    }
    return (tuple){&s.a[k], &c.a[k], u.size};
}

tuple t_tan_sec2 (series t, series s, series u, int k, geometry g) {
    assert(t.a != s.a && t.a != u.a && s.a != u.a);
    assert(k >= 0);
    if (k == 0) {
        g == TRIG ? mpfr_tan(t.a[k], u.a[0], RND) : mpfr_tanh(t.a[k], u.a[0], RND);
        g == TRIG ? mpfr_sec(s.a[k], u.a[0], RND) : mpfr_sech(s.a[k], u.a[0], RND);
        mpfr_sqr(s.a[k], s.a[k], RND);
    } else {
        d_cauchy(&t.a[k], s, u, k, 0, k - 1, D1);
        d_cauchy(&s.a[k], t, t, k, 0, k - 1, g == TRIG ? D2 : D_2);
    }
    return (tuple){&t.a[k], &s.a[k], u.size};
}

mpfr_t *t_pwr (series p, series u, mpfr_t a, int k) {
    assert(mpfr_sgn(u.a[0]) > 0);
    assert(p.a != u.a);
    assert(k >= 0);
    if (k == 0) {
        mpfr_pow(p.a[k], u.a[0], a, RND);
    } else {
        mpfr_sub(_, *d_cauchy(&__, p, u, k, 0, k - 1, a), *d_cauchy(&___, u, p, k, 1, k - 1, D1), RND);
        mpfr_div(p.a[k], _, u.a[0], RND);
    }
    return &p.a[k];
}

mpfr_t *t_ln (series l, series u, int k) {
    assert(mpfr_sgn(u.a[0]) > 0);
    assert(l.a != u.a);
    assert(k >= 0);
    if (k == 0) {
        mpfr_log(l.a[k], u.a[0], RND);
    } else {
        mpfr_sub(_, u.a[k], *d_cauchy(&__, u, l, k, 1, k - 1, D1), RND);
        mpfr_div(l.a[k], _, u.a[0], RND);
    }
    return &l.a[k];
}

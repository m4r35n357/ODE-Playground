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

static char fs[42];

static mpfr_t D1, D_1, D2, D_2, _, __, ___;

void t_tempvars (void) {
    mpfr_inits(_, __, ___, NULL);
    mpfr_init_set_ui(D1, 1, RND);
    mpfr_init_set_si(D_1, -1, RND);
    mpfr_init_set_ui(D2, 2, RND);
    mpfr_init_set_si(D_2, -2, RND);
}

void t_output (mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t h, long step) {
    mpfr_mul_ui(_, h, step, RND);
    mpfr_printf(fs, x, y, z, _);
}

void t_stepper (char **argv, long *n, mpfr_t *h, long *nsteps) {
    long d = strtol(argv[1], NULL, BASE);
    d == 0 ? sprintf(fs, "%%.RNe %%.RNe %%.RNe %%.9RNe\n") : sprintf(fs, "%%+.%luRNe %%+.%luRNe %%+.%luRNe %%+.9RNe\n", d, d, d);
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
        mpfr_set_str(*va_arg(vars, mpfr_t *), argv[i], BASE, RND);
    }
    va_end(vars);
}

series t_series (int n) {
    assert(n > 0);
    mpfr_t *jet = malloc(sizeof (mpfr_t) * n);
    for (int i = 0; i < n; i++) {
        mpfr_init_set_ui(jet[i], 0, RND);
    }
    return (series){jet, n};
}

void t_next (series s, mpfr_t dot, int k, sign sgn) {
    sgn == NEG ? mpfr_div_si(s.jet[k + 1], dot, - (k + 1), RND) : mpfr_div_ui(s.jet[k + 1], dot, k + 1, RND);
}

mpfr_t *t_horner (series s, mpfr_t h) {
    mpfr_set_zero(_, 1);
    for (int i = s.size - 1; i >= 0; i--) {
        mpfr_fma(_, _, h, s.jet[i], RND);
    }
    if (mpfr_number_p(_) == 0) exit(1);
    mpfr_swap(s.jet[0], _);
    return s.jet;
}

mpfr_t *t_abs (series u, int k) {
    assert(k >= 0);
    mpfr_sgn(u.jet[0]) < 0 ? mpfr_neg(_, u.jet[k], RND) : mpfr_set(_, u.jet[k], RND);
    return &_;
}

static mpfr_t *cauchy (mpfr_t *ck, series a, series b, int k, int lower, int upper) {
    mpfr_set_zero(*ck, 1);
    for (int j = lower; j <= upper; j++) {
        mpfr_fma(*ck, a.jet[j], b.jet[k - j], *ck, RND);
    }
    return ck;
}

mpfr_t *t_prod (series u, series v, int k) {
    assert(k >= 0);
    return cauchy(&__, u, v, k, 0, k);
}

mpfr_t *t_quot (series q, series u, series v, int k) {
    assert(mpfr_zero_p(v.jet[0]) == 0);
    assert(q.jet != u.jet && q.jet != v.jet);
    assert(k >= 0);
    k == 0 ? mpfr_set(_, u.jet[0], RND) : mpfr_sub(_, u.jet[k], *cauchy(&__, q, v, k, 0, k - 1), RND);
    mpfr_div(q.jet[k], _, v.jet[0], RND);
    return &q.jet[k];
}

mpfr_t *t_inv (series i, series v, int k) {
    assert(mpfr_zero_p(v.jet[0]) == 0);
    assert(i.jet != v.jet);
    assert(k >= 0);
    k == 0 ? mpfr_set(_, D1, RND) : mpfr_neg(_, *cauchy(&__, i, v, k, 0, k - 1), RND);
    mpfr_div(i.jet[k], _, v.jet[0], RND);
    return &i.jet[k];
}

mpfr_t *t_sqr (series u, int k) {
    assert(k >= 0);
    return cauchy(&__, u, u, k, 0, k);
}

mpfr_t *t_sqrt (series r, series u, int k) {
    assert(mpfr_sgn(u.jet[0]) > 0);
    assert(r.jet != u.jet);
    assert(k >= 0);
    if (k == 0) {
        mpfr_sqrt(r.jet[k], u.jet[0], RND);
    } else {
        mpfr_sub(_, u.jet[k], *cauchy(&__, r, r, k, 1, k - 1), RND);
        mpfr_div_2ui(_, _, 1, RND);
        mpfr_div(r.jet[k], _, r.jet[0], RND);
    }
    return &r.jet[k];
}

static mpfr_t *d_cauchy (mpfr_t *fk, series h, series u, int k, int lower, int upper, mpfr_t factor) {
    assert(fk != &_ && h.jet != &_ && u.jet != &_);  // _ is used internally so it cannot be a parameter
    mpfr_set_zero(*fk, 1);
    for (int j = lower; j <= upper; j++) {
        mpfr_mul_ui(_, u.jet[k - j], k - j, RND);
        mpfr_fma(*fk, h.jet[j], _, *fk, RND);
    }
    mpfr_div_ui(*fk, *fk, k, RND);
    mpfr_mul(*fk, *fk, factor, RND);
    return fk;
}

mpfr_t *t_exp (series e, series u, int k) {
    assert(e.jet != u.jet);
    assert(k >= 0);
    if (k == 0) {
        mpfr_exp(e.jet[k], u.jet[0], RND);
    } else {
        d_cauchy(&e.jet[k], e, u, k, 0, k - 1, D1);
    }
    return &e.jet[k];
}

tuple t_sin_cos (series s, series c, series u, int k, geometry g) {
    assert(s.jet != c.jet && s.jet != u.jet && c.jet != u.jet);
    assert(k >= 0);
    if (k == 0) {
        g == TRIG ? mpfr_sin_cos(s.jet[k], c.jet[k], u.jet[0], RND) : mpfr_sinh_cosh(s.jet[k], c.jet[k], u.jet[0], RND);
    } else {
        d_cauchy(&s.jet[k], c, u, k, 0, k - 1, D1);
        d_cauchy(&c.jet[k], s, u, k, 0, k - 1, g == TRIG ? D_1 : D1);
    }
    return (tuple){&s.jet[k], &c.jet[k], u.size};
}

tuple t_tan_sec2 (series t, series s, series u, int k, geometry g) {
    assert(t.jet != s.jet && t.jet != u.jet && s.jet != u.jet);
    assert(k >= 0);
    if (k == 0) {
        g == TRIG ? mpfr_tan(t.jet[k], u.jet[0], RND) : mpfr_tanh(t.jet[k], u.jet[0], RND);
        g == TRIG ? mpfr_sec(s.jet[k], u.jet[0], RND) : mpfr_sech(s.jet[k], u.jet[0], RND);
        mpfr_sqr(s.jet[k], s.jet[k], RND);
    } else {
        d_cauchy(&t.jet[k], s, u, k, 0, k - 1, D1);
        d_cauchy(&s.jet[k], t, t, k, 0, k - 1, g == TRIG ? D2 : D_2);
    }
    return (tuple){&t.jet[k], &s.jet[k], u.size};
}

mpfr_t *t_pwr (series p, series u, mpfr_t a, int k) {
    assert(mpfr_sgn(u.jet[0]) > 0);
    assert(p.jet != u.jet);
    assert(k >= 0);
    if (k == 0) {
        mpfr_pow(p.jet[k], u.jet[0], a, RND);
    } else {
        mpfr_sub(_, *d_cauchy(&__, p, u, k, 0, k - 1, a), *d_cauchy(&___, u, p, k, 1, k - 1, D1), RND);
        mpfr_div(p.jet[k], _, u.jet[0], RND);
    }
    return &p.jet[k];
}

mpfr_t *t_ln (series l, series u, int k) {
    assert(mpfr_sgn(u.jet[0]) > 0);
    assert(l.jet != u.jet);
    assert(k >= 0);
    if (k == 0) {
        mpfr_log(l.jet[k], u.jet[0], RND);
    } else {
        mpfr_sub(_, u.jet[k], *d_cauchy(&__, u, l, k, 1, k - 1, D1), RND);
        mpfr_div(l.jet[k], _, u.jet[0], RND);
    }
    return &l.jet[k];
}

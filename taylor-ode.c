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

static mpfr_t _, __, ___;

void t_tempvars (void) {
    mpfr_inits(_, __, ___, NULL);
}

void t_xyz_output (mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t t) {
    mpfr_printf("%+.12RNe %+.12RNe %+.12RNe %+.6RNe\n", x, y, z, t);
}

void t_stepper (char **argv, long *n, mpfr_t *t, mpfr_t *h, long *nsteps) {
    mpfr_set_default_prec(strtod(argv[1], NULL) * 3.322);
    fprintf(stderr, " MPFR default precision: %lu bits\n", mpfr_get_default_prec());
    *n = strtol(argv[2], NULL, BASE);
    mpfr_init_set_ui(*t, 0, RND);
    mpfr_init_set_str(*h, argv[3], BASE, RND);
    *nsteps = strtol(argv[4], NULL, BASE);
    t_tempvars();
}

void t_args (char **argv, int argc, ...) {
    va_list vars;
    va_start(vars, argc);
    for (int i = 5; i < argc; i++) {
        mpfr_init_set_str(*va_arg(vars, mpfr_t *), argv[i], BASE, RND);
    }
    va_end(vars);
}

mpfr_t *t_jet (int n) {
    assert(n > 0);
    mpfr_t *jet = malloc(sizeof (mpfr_t) * n);
    for (int i = 0; i < n; i++) {
        mpfr_init(jet[i]);
        mpfr_set_zero(jet[i], 1);
    }
    return jet;
}

mpfr_t *t_jet_c (int n, mpfr_t value) {
    mpfr_t *jet = t_jet(n);
    mpfr_set(jet[0], value, RND);
    return jet;
}

void t_horner (mpfr_t *jet, int n, mpfr_t h) {
    assert(n > 0);
    mpfr_set_zero(_, 1);
    for (int i = n; i > - 1; i--) {
        mpfr_fma(_, _, h, jet[i], RND);
    }
    assert(mpfr_number_p(_));
    mpfr_swap(jet[0], _);
}

mpfr_t *t_abs (mpfr_t *u, int k) {
    assert(k >= 0);
    mpfr_sgn(u[0]) < 0 ? mpfr_neg(_, u[k], RND) : mpfr_set(_, u[k], RND);
    return &_;
}

static mpfr_t *cauchy (mpfr_t *c, mpfr_t *a, mpfr_t *b, int k, int lower, int upper) {
    mpfr_set_zero(*c, 1);
    for (int j = lower; j <= upper; j++) {
        mpfr_fma(*c, a[j], b[k - j], *c, RND);
    }
    return c;
}

mpfr_t *t_prod (mpfr_t *u, mpfr_t *v, int k) {
    assert(k >= 0);
    return cauchy(&__, u, v, k, 0, k);
}

mpfr_t *t_quot (mpfr_t *q, mpfr_t *u, mpfr_t *v, int k) {
    assert(mpfr_zero_p(v[0]) == 0);
    assert(q != u && q != v && u != v);
    assert(k >= 0);
    mpfr_sub(_, u[k], *cauchy(&__, q, v, k, 0, k - 1), RND);
    mpfr_div(q[k], _, v[0], RND);
    return &q[k];
}

mpfr_t *t_sqr (mpfr_t *u, int k) {
    assert(k >= 0);
    mpfr_mul_2ui(_, *cauchy(&__, u, u, k, 0, (k - (k % 2 == 0 ? 2 : 1)) / 2), 1, RND);
    if (k % 2 == 0) mpfr_fma(_, u[k / 2], u[k / 2], _, RND);
    return &_;
}

mpfr_t *t_sqrt (mpfr_t *r, mpfr_t *u, int k) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(r != u);
    assert(k >= 0);
    if (k == 0) {
        mpfr_sqrt(r[k], u[0], RND);
    } else {
        mpfr_mul_2ui(_, *cauchy(&__, r, r, k, 1, (k - (k % 2 == 0 ? 2 : 1)) / 2), 1, RND);
        if (k % 2 == 0) mpfr_fma(_, r[k / 2], r[k / 2], _, RND);
        mpfr_sub(_, u[k], _, RND);
        mpfr_div_2ui(_, _, 1, RND);
        mpfr_div(r[k], _, r[0], RND);
    }
    return &r[k];
}

static mpfr_t *d_cauchy (mpfr_t *f, mpfr_t *h, mpfr_t *u, int k, int lower, int upper, double factor) {
    assert(f != &_);
    mpfr_set_zero(*f, 1);
    for (int j = lower; j <= upper; j++) {
        mpfr_mul_ui(_, u[k - j], k - j, RND);
        mpfr_fma(*f, h[j], _, *f, RND);
    }
    mpfr_div_ui(*f, *f, k, RND);
    mpfr_mul_d(*f, *f, factor, RND);
    return f;
}

mpfr_t *t_exp (mpfr_t *e, mpfr_t *u, int k) {
    assert(e != u);
    assert(k >= 0);
    if (k == 0) {
        mpfr_exp(e[k], u[0], RND);
    } else {
        d_cauchy(&e[k], e, u, k, 0, k - 1, 1.0);
    }
    return &e[k];
}

tuple t_sin_cos (mpfr_t *s, mpfr_t *c, mpfr_t *u, int k, geometry g) {
    assert(s != c && s != u && c != u);
    assert(k >= 0);
    if (k == 0) {
        g == TRIG ? mpfr_sin_cos(s[k], c[k], u[0], RND) : mpfr_sinh_cosh(s[k], c[k], u[0], RND);
    } else {
        d_cauchy(&s[k], c, u, k, 0, k - 1, 1.0);
        d_cauchy(&c[k], s, u, k, 0, k - 1, g == TRIG ? -1.0 : 1.0);
    }
    return (tuple){&s[k], &c[k]};
}

tuple t_tan_sec2 (mpfr_t *t, mpfr_t *s, mpfr_t *u, int k, geometry g) {
    assert(t != s && t != u && s != u);
    assert(k >= 0);
    if (k == 0) {
        g == TRIG ? mpfr_tan(t[k], u[0], RND) : mpfr_tanh(t[k], u[0], RND);
        g == TRIG ? mpfr_sec(s[k], u[0], RND) : mpfr_sech(s[k], u[0], RND);
        mpfr_sqr(s[0], s[0], RND);
    } else {
        d_cauchy(&t[k], s, u, k, 0, k - 1, 1.0);
        d_cauchy(&s[k], t, t, k, 0, k - 1, g == TRIG ? 2.0 : -2.0);
    }
    return (tuple){&t[k], &s[k]};
}

mpfr_t *t_pwr (mpfr_t *p, mpfr_t *u, double a, int k) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(p != u);
    assert(k >= 0);
    if (k == 0) {
        mpfr_set_d(_, a, RND);
        mpfr_pow(p[k], u[0], _, RND);
    } else {
        mpfr_sub(_, *d_cauchy(&__, p, u, k, 0, k - 1, a), *d_cauchy(&___, u, p, k, 1, k - 1, 1.0), RND);
        mpfr_div(p[k], _, u[0], RND);
    }
    return &p[k];
}

mpfr_t *t_ln (mpfr_t *l, mpfr_t *u, int k) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(l != u);
    assert(k >= 0);
    if (k == 0) {
        mpfr_log(l[k], u[0], RND);
    } else {
        mpfr_sub(_, u[k], *d_cauchy(&__, u, l, k, 1, k - 1, 1.0), RND);
        mpfr_div(l[k], _, u[0], RND);
    }
    return &l[k];
}

/*
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

const int BASE = 10;

const mpfr_rnd_t RND = MPFR_RNDN;

void t_arg (char **argv, int arg, mpfr_t *dest) {
    mpfr_init_set_str(*dest, argv[arg], BASE, RND);
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
}

mpfr_t *t_jet (int n) {
    assert(n > 0);
    mpfr_t *jet = malloc(sizeof (mpfr_t) * n);
    for (int i = 0; i < n; i++) {
        mpfr_init(jet[i]);
    }
    return jet;
}

mpfr_t *t_jet_c (int n, mpfr_t value) {
    mpfr_t *jet = t_jet(n);
    mpfr_set(jet[0], value, RND);
    for (int i = 1; i < n; i++) {
         mpfr_set_zero(jet[i], 1);
    }
    return jet;
}

void t_horner (mpfr_t *sum, mpfr_t *jet, int n, mpfr_t h) {
    assert(n > 0);
    mpfr_set_zero(*sum, 1);
    for (int i = n; i > - 1; i--) {
        mpfr_fma(*sum, *sum, h, jet[i], RND);
    }
    assert(mpfr_number_p(*sum));
}

mpfr_t *t_abs (mpfr_t *a, mpfr_t *u, int k) {
    assert(a != u);
    assert(k >= 0);
    mpfr_mul_si(*a, u[k], mpfr_sgn(u[0]), RND);
    return a;
}

static mpfr_t *cauchy (mpfr_t *c, mpfr_t *a, mpfr_t *b, int k, int lower, int upper) {
    mpfr_set_zero(*c, 1);
    for (int j = lower; j < upper; j++) {
        mpfr_fma(*c, a[j], b[k - j], *c, RND);
    }
    return c;
}

mpfr_t *t_sqr (mpfr_t *s, mpfr_t *u, int k) {
    assert(s != u);
    assert(k >= 0);
    mpfr_mul_2ui(*s, *cauchy(s, u, u, k, 0, (k - (k % 2 == 0 ? 2 : 1)) / 2 + 1), 1, RND);
    if (k % 2 == 0) mpfr_fma(*s, u[k / 2], u[k / 2], *s, RND);
    return s;
}

mpfr_t *t_prod (mpfr_t *p, mpfr_t *u, mpfr_t *v, int k) {
    assert(p != u && p != v);
    assert(k >= 0);
    return cauchy(p, u, v, k, 0, k + 1);
}

mpfr_t *t_quot (mpfr_t *q, mpfr_t *u, mpfr_t *v, int k) {
    assert(mpfr_sgn(v[0]) != 0);
    assert(q != u && q != v && u != v);
    assert(k >= 0);
    mpfr_sub(q[k], u[k], *cauchy(&q[k], q, v, k, 0, k), RND);
    mpfr_div(q[k], q[k], v[0], RND);
    return &q[k];
}

mpfr_t *t_sqrt (mpfr_t *r, mpfr_t *u, int k) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(r != u);
    assert(k >= 0);
    if (k == 0) {
        mpfr_sqrt(r[0], u[0], RND);
    } else {
        mpfr_mul_2ui(r[k], *cauchy(&r[k], r, r, k, 0, (k - (k % 2 == 0 ? 2 : 1)) / 2 + 1), 1, RND);
        if (k % 2 == 0) mpfr_fma(r[k], r[k / 2], r[k / 2], r[k], RND);
        mpfr_sub(r[k], u[k], r[k], RND);
        mpfr_div_2ui(r[k], r[k], 1, RND);
        mpfr_div(r[k], r[k], r[0], RND);
    }
    return &r[k];
}

static mpfr_t *d_cauchy (mpfr_t *f, mpfr_t *h, mpfr_t *u, int k, double factor, int lower, int upper, mpfr_t *_) {
    mpfr_set_zero(*f, 1);
    for (int j = lower; j < upper; j++) {
        mpfr_mul_d(*_, u[k - j], factor * (k - j) / k, RND);
        mpfr_fma(*f, h[j], *_, *f, RND);
    }
    return f;
}

mpfr_t *t_exp (mpfr_t *e, mpfr_t *u, int k, mpfr_t *_) {
    assert(e != u);
    assert(_ != e && _ != u);
    assert(k >= 0);
    if (k == 0) {
        mpfr_exp(e[0], u[0], RND);
        return &e[0];
    }
    return d_cauchy(&e[k], e, u, k, 1.0, 0, k, _);
}

tuple t_sin_cos (mpfr_t *s, mpfr_t *c, mpfr_t *u, int k, mpfr_t *_, geometry g) {
    assert(s != c && s != u && c != u);
    assert(_ != s && _ != c && _ != u);
    assert(k >= 0);
    if (k == 0) {
        (g == TRIG) ? mpfr_sin_cos(s[0], c[0], u[0], RND) : mpfr_sinh_cosh(s[0], c[0], u[0], RND);
        return (tuple){&s[0], &c[0]};
    }
    return (tuple){d_cauchy(&s[k], c, u, k, 1.0, 0, k, _), d_cauchy(&c[k], s, u, k, (g == TRIG) ? - 1.0 : 1.0, 0, k, _)};
}

tuple t_tan_sec2 (mpfr_t *t, mpfr_t *s2, mpfr_t *u, int k, mpfr_t *_, geometry g) {
    assert(t != s2 && t != u && s2 != u);
    assert(_ != t && _ != s2 && _ != u);
    assert(k >= 0);
    if (k == 0) {
        if (g == TRIG) {
            mpfr_tan(t[0], u[0], RND); mpfr_sec(s2[0], u[0], RND);
        } else {
            mpfr_tanh(t[0], u[0], RND); mpfr_sech(s2[0], u[0], RND);
        }
        mpfr_sqr(s2[0], s2[0], RND);
        return (tuple){&t[0], &s2[0]};
    }
    return (tuple){d_cauchy(&t[k], s2, u, k, 1.0, 0, k, _), d_cauchy(&s2[k], t, t, k, (g == TRIG) ? 2.0 : - 2.0, 0, k, _)};
}

mpfr_t *t_pwr (mpfr_t *p, mpfr_t *u, double a, int k, mpfr_t *_) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(p != u);
    assert(_ != p && _ != u);
    assert(k >= 0);
    if (k == 0) {
        mpfr_set_d(*_, a, RND);
        mpfr_pow(p[0], u[0], *_, RND);
    } else {
        mpfr_set_zero(p[k], 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_d(*_, u[k - j], a * (k - j) - j, RND);
            mpfr_fma(p[k], p[j], *_, p[k], RND);
        }
        mpfr_div_ui(p[k], p[k], k, RND);
        mpfr_div(p[k], p[k], u[0], RND);
    }
    return &p[k];
}

mpfr_t *t_ln (mpfr_t *l, mpfr_t *u, int k, mpfr_t *_) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(l != u);
    assert(_ != l && _ != u);
    assert(k >= 0);
    if (k == 0) {
        mpfr_log(l[0], u[0], RND);
    } else {
        mpfr_sub(*_, u[k], *d_cauchy(&l[k], u, l, k, 1.0, 1, k, _), RND);
        mpfr_div(l[k], *_, u[0], RND);
    }
    return &l[k];
}

/*
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include <mpfr.h>
#include "taylor-ode.h"

const int BASE = 10;

const mpfr_rnd_t RND = MPFR_RNDN;

void t_line_output (mpfr_t t, int count, ...) {
    va_list data;
    va_start(data, count);
    for (int i = 0; i < count; i++) {
        mpfr_vprintf("%.9RNe ", data);
    }
    va_end(data);
    mpfr_printf("%.5RNe\n", t);
}

void t_stepper (char **argv, long *n, mpfr_t *t, mpfr_t *h, long *nsteps) {
    mpfr_set_default_prec(strtod(argv[1], NULL) * 3.32);
    *n = strtol(argv[2], NULL, BASE);
    assert(*n > 1);
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
    assert(n > 0);
    mpfr_t *jet = t_jet(n);
    mpfr_set(jet[0], value, RND);
    for (int i = 1; i < n; i++) {
         mpfr_set_zero(jet[i], 1);
    }
    return jet;
}

void t_horner (mpfr_t *sum, const mpfr_t *jet, long n, mpfr_t h) {
    assert(sizeof *sum == sizeof (mpfr_t));
    assert(n > 0);
    mpfr_set(*sum, jet[n], RND);
    for (int j = n - 1; j > - 1; j--) {
        mpfr_fma(*sum, *sum, h, jet[j], RND);
    }
}

void t_square (mpfr_t *s, const mpfr_t *u, int k) {
    assert(sizeof *s == sizeof (mpfr_t));
    assert(k >= 0);
    if (k == 0) {
        mpfr_sqr(*s, u[0], RND);
    } else {
        mpfr_set_zero(*s, 1);
        if (k % 2 == 1) {
            for (int j = 0; j < (k - 1) / 2 + 1; j++) {
                mpfr_fma(*s, u[j], u[k - j], *s, RND);
            }
            mpfr_mul_2ui(*s, *s, 1, RND);
        } else {
            for (int j = 0; j < (k - 2) / 2 + 1; j++) {
                mpfr_fma(*s, u[j], u[k - j], *s, RND);
            }
            mpfr_mul_2ui(*s, *s, 1, RND);
            mpfr_fma(*s, u[k / 2], u[k / 2], *s, RND);
        }
    }
}

void t_product (mpfr_t *p, const mpfr_t *u, const mpfr_t *v, int k) {
    assert(sizeof *p == sizeof (mpfr_t));
    assert(sizeof *u == sizeof *v);
    assert(k >= 0);
    mpfr_set_zero(*p, 1);
    for (int j = 0; j < k + 1; j++) {
        mpfr_fma(*p, u[j], v[k - j], *p, RND);
    }
}

void t_quotient (mpfr_t *q, const mpfr_t *u, const mpfr_t *v, int k) {
    assert(mpfr_sgn(v[0]) != 0);
    assert(q != u && q != v && u != v);
    assert(sizeof *u == sizeof *q && sizeof *v == sizeof *q);
    assert(k >= 0);
    mpfr_set_zero(q[k], 1);
    for (int j = 0; j < k + 1; j++) {
        mpfr_fma(q[k], v[j], q[k - j], q[k], RND);
    }
    mpfr_sub(q[k], u[k], q[k], RND);
    mpfr_div(q[k], q[k], v[0], RND);
}

void t_sqrt (mpfr_t *r, const mpfr_t *u, int k) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(r != u);
    assert(sizeof *u == sizeof *r);
    assert(k >= 0);
    if (k == 0) {
        mpfr_sqrt(r[0], u[0], RND);
    } else {
        mpfr_set_zero(r[k], RND);
        if (k % 2 == 1) {
            for (int j = 0; j < (k - 1) / 2 + 1; j++) {
                mpfr_fma(r[k], r[j], r[k - j], r[k], RND);
            }
            mpfr_mul_ui(r[k], r[k], 2, RND);
        } else {
            for (int j = 0; j < (k - 2) / 2 + 1; j++) {
                mpfr_fma(r[k], r[j], r[k - j], r[k], RND);
            }
            mpfr_mul_ui(r[k], r[k], 2, RND);
            mpfr_fma(r[k], r[k / 2], r[k / 2], r[k], RND);
        }
        mpfr_sub(r[k], u[k], r[k], RND);
        mpfr_div_ui(r[k], r[k], 2, RND);
        mpfr_div(r[k], r[k], r[0], RND);
    }
}

static void _ddot (mpfr_t *d, const mpfr_t *v, const mpfr_t *u, int k, mpfr_t *_) {
    assert(d != _);
    assert(sizeof *d == sizeof (mpfr_t));
    assert(sizeof *u == sizeof *v);
    assert(sizeof *_ == sizeof (mpfr_t));
    assert(k > 0);
    mpfr_set_zero(*d, 1);
    for (int j = 1; j < k; j++) {
        mpfr_mul_ui(*_, u[j], j, RND);
        mpfr_fma(*d, *_, v[k - j], *d, RND);
    }
    mpfr_div_ui(*d, *d, k, RND);
}

void t_exp (mpfr_t *e, const mpfr_t *u, int k, mpfr_t *_) {
    assert(e != u);
    assert(sizeof *u == sizeof *e);
    assert(sizeof *_ == sizeof (mpfr_t));
    assert(k >= 0);
    if (k == 0) {
        mpfr_exp(e[0], u[0], RND);
    } else {
        _ddot(&e[k], e, u, k, _);
        mpfr_fma(e[k], e[0], u[k], e[k], RND);
    }
}

void t_sin_cos (mpfr_t *s, mpfr_t *c, const mpfr_t *u, int k, mpfr_t *_, geometry g) {
    assert(s != c && s != u && c != u);
    assert(sizeof *u == sizeof *s && sizeof *u == sizeof *c);
    assert(sizeof *_ == sizeof (mpfr_t));
    assert(k >= 0);
    if (k == 0) {
        if (g == TRIG) {
            mpfr_sin_cos(s[0], c[0], u[0], RND);
        } else {
            mpfr_sinh_cosh(s[0], c[0], u[0], RND);
        }
    } else {
        _ddot(&s[k], c, u, k, _);
        mpfr_fma(s[k], c[0], u[k], s[k], RND);
        _ddot(&c[k], s, u, k, _);
        mpfr_fma(c[k], s[0], u[k], c[k], RND);
        if (g == TRIG) {
            mpfr_neg(c[k], c[k], RND);
        }
    }
}

void t_tan_sec2 (mpfr_t *t, mpfr_t *s2, const mpfr_t *u, int k, mpfr_t *_, geometry g) {
    assert(t != s2 && t != u && s2 != u);
    assert(sizeof *u == sizeof *t && sizeof *u == sizeof *s2);
    assert(sizeof *_ == sizeof (mpfr_t));
    assert(k >= 0);
    if (k == 0) {
        if (g == TRIG) {
            mpfr_tan(t[0], u[0], RND);
            mpfr_sqr(*_, t[0], RND);
            mpfr_add_ui(s2[0], *_, 1, RND);
        } else {
            mpfr_tanh(t[0], u[0], RND);
            mpfr_sqr(*_, t[0], RND);
            mpfr_ui_sub(s2[0], 1, *_, RND);
        }
    } else {
        _ddot(&t[k], s2, u, k, _);
        mpfr_fma(t[k], s2[0], u[k], t[k], RND);
        _ddot(&s2[k], t, t, k, _);
        mpfr_fma(s2[k], t[0], t[k], s2[k], RND);
        mpfr_mul_2ui(s2[k], s2[k], 1, RND);
        if (g == HYP) {
            mpfr_neg(s2[k], s2[k], RND);
        }
    }
}

void t_power (mpfr_t *p, const mpfr_t *u, mpfr_t a, int k, mpfr_t *_, mpfr_t *__) {
    assert(mpfr_sgn(u[0]) != 0);
    assert(p != u);
    assert(sizeof *u == sizeof *p);
    assert(sizeof *_ == sizeof (mpfr_t));
    assert(sizeof *__ == sizeof (mpfr_t));
    assert(k >= 0);
    mpfr_set_zero(p[k], 1);
    if (k == 0) {
        mpfr_pow(p[0], u[0], a, RND);
    } else {
        _ddot(__, p, u, k, _);
        mpfr_fma(*__, p[0], u[k], *__, RND);
        mpfr_mul(p[k], *__, a, RND);
        _ddot(__, u, p, k, _);
        mpfr_sub(p[k], p[k], *__, RND);
        mpfr_div(p[k], p[k], u[0], RND);
    }
}

void t_ln (mpfr_t *l, const mpfr_t *u, int k, mpfr_t *_) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(_ != l && _ != u);
    assert(l != u);
    assert(sizeof *_ == sizeof (mpfr_t));
    assert(k >= 0);
    if (k == 0) {
        mpfr_log(l[0], u[0], RND);
    } else {
        _ddot(&l[k], u, l, k, _);
        mpfr_sub(l[k], u[k], l[k], RND);
        mpfr_div(l[k], l[k], u[0], RND);
    }
}

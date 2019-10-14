/*
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

const int BASE = 10;

const mpfr_rnd_t RND = MPFR_RNDN;

void t_xyz_output (mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t t) {
    mpfr_printf("%.9RNe %.9RNe %.9RNe %.5RNe\n", x, y, z, t);
}

void t_stepper (char **argv, long *n, mpfr_t *t, mpfr_t *h, long *nsteps) {
    mpfr_set_default_prec(strtod(argv[1], NULL) * 3.32);
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

void t_horner (mpfr_t *sum, const mpfr_t *jet, int n, mpfr_t h) {
    assert(n > 0);
    mpfr_set_zero(*sum, 1);
    for (int i = n; i > - 1; i--) {
        mpfr_fma(*sum, *sum, h, jet[i], RND);
    }
}

void t_abs (mpfr_t *a, const mpfr_t *u, int k) {
    assert(a != u);
    assert(k >= 0);
    mpfr_mul_si(*a, u[k], mpfr_sgn(u[0]), RND);
}

void t_sqr (mpfr_t *s, const mpfr_t *u, int k) {
    assert(s != u);
    assert(k >= 0);
    int upper = (k - (k % 2 == 0 ? 2 : 1)) / 2 + 1;
    mpfr_set_zero(*s, 1);
    for (int j = 0; j < upper; j++) {
        mpfr_fma(*s, u[j], u[k - j], *s, RND);
    }
    mpfr_mul_2ui(*s, *s, 1, RND);
    if (k % 2 == 0) {
        mpfr_fma(*s, u[k / 2], u[k / 2], *s, RND);
    }
}

void t_prod (mpfr_t *p, const mpfr_t *u, const mpfr_t *v, int k) {
    assert(p != u && p != v);
    assert(k >= 0);
    mpfr_set_zero(*p, 1);
    for (int j = 0; j < k + 1; j++) {
        mpfr_fma(*p, u[j], v[k - j], *p, RND);
    }
}

void t_quot (mpfr_t *q, const mpfr_t *u, const mpfr_t *v, int k) {
    assert(mpfr_sgn(v[0]) != 0);
    assert(q != u && q != v && u != v);
    assert(k >= 0);
    mpfr_set_zero(q[k], 1);
    for (int j = 0; j < k; j++) {
        mpfr_fma(q[k], q[j], v[k - j], q[k], RND);
    }
    mpfr_sub(q[k], u[k], q[k], RND);
    mpfr_div(q[k], q[k], v[0], RND);
}

void t_sqrt (mpfr_t *r, const mpfr_t *u, int k) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(r != u);
    assert(k >= 0);
    if (k == 0) {
        mpfr_sqrt(r[0], u[0], RND);
    } else {
        int upper = (k - (k % 2 == 0 ? 2 : 1)) / 2 + 1;
        mpfr_set_zero(r[k], RND);
        for (int j = 1; j < upper; j++) {
            mpfr_fma(r[k], r[j], r[k - j], r[k], RND);
        }
        mpfr_mul_2ui(r[k], r[k], 1, RND);
        if (k % 2 == 0) {
            mpfr_fma(r[k], r[k / 2], r[k / 2], r[k], RND);
        }
        mpfr_sub(r[k], u[k], r[k], RND);
        mpfr_div_2ui(r[k], r[k], 1, RND);
        mpfr_div(r[k], r[k], r[0], RND);
    }
}

void t_pwr (mpfr_t *p, const mpfr_t *u, mpfr_t a, int k, mpfr_t *_) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(p != u);
    assert(_ != p && _ != u);
    assert(k >= 0);
    if (k == 0) {
        mpfr_pow(p[0], u[0], a, RND);
    } else {
        mpfr_set_zero(p[k], 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_ui(*_, a, k - j, RND);
            mpfr_sub_ui(*_, *_, j, RND);
            mpfr_mul(*_, *_, p[j], RND);
            mpfr_fma(p[k], *_, u[k - j], p[k], RND);
        }
        mpfr_div_ui(p[k], p[k], k, RND);
        mpfr_div(p[k], p[k], u[0], RND);
    }
}

void t_exp (mpfr_t *e, const mpfr_t *u, int k, mpfr_t *_) {
    assert(e != u);
    assert(_ != e && _ != u);
    assert(k >= 0);
    if (k == 0) {
        mpfr_exp(e[0], u[0], RND);
    } else {
        mpfr_set_zero(e[k], 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_d(*_, e[j], (k - j) / (double)k, RND);
            mpfr_fma(e[k], *_, u[k - j], e[k], RND);
        }
    }
}

void t_ln (mpfr_t *l, const mpfr_t *u, int k, mpfr_t *_) {
    assert(mpfr_sgn(u[0]) > 0);
    assert(l != u);
    assert(_ != l && _ != u);
    assert(k >= 0);
    if (k == 0) {
        mpfr_log(l[0], u[0], RND);
    } else {
        mpfr_set_zero(l[k], 1);
        for (int j = 1; j < k; j++) {
            mpfr_mul_d(*_, l[j], j / (double)k, RND);
            mpfr_fma(l[k], *_, u[k - j], l[k], RND);
        }
        mpfr_sub(l[k], u[k], l[k], RND);
        mpfr_div(l[k], l[k], u[0], RND);
    }
}

void t_sin_cos (mpfr_t *s, mpfr_t *c, const mpfr_t *u, int k, mpfr_t *_, geometry g) {
    assert(s != c && s != u && c != u);
    assert(_ != s && _ != c && _ != u);
    assert(k >= 0);
    if (k == 0) {
        if (g == TRIG) {
            mpfr_sin_cos(s[0], c[0], u[0], RND);
        } else {
            mpfr_sinh_cosh(s[0], c[0], u[0], RND);
        }
    } else {
        mpfr_set_zero(s[k], 1);
        mpfr_set_zero(c[k], 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_d(*_, c[j], (k - j) / (double)k, RND);
            mpfr_fma(s[k], *_, u[k - j], s[k], RND);
            mpfr_mul_d(*_, s[j], (k - j) / (double)k, RND);
            mpfr_fma(c[k], *_, u[k - j], c[k], RND);
        }
        if (g == TRIG) {
            mpfr_neg(c[k], c[k], RND);
        }
    }
}

void t_tan_sec2 (mpfr_t *t, mpfr_t *s2, const mpfr_t *u, int k, mpfr_t *_, geometry g) {
    assert(t != s2 && t != u && s2 != u);
    assert(_ != t && _ != s2 && _ != u);
    assert(k >= 0);
    if (k == 0) {
        if (g == TRIG) {
            mpfr_tan(t[0], u[0], RND);
            mpfr_sec(s2[0], u[0], RND);
        } else {
            mpfr_tanh(t[0], u[0], RND);
            mpfr_sech(s2[0], u[0], RND);
        }
        mpfr_sqr(s2[0], s2[0], RND);
    } else {
        mpfr_set_zero(t[k], 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_d(*_, s2[j], (k - j) / (double)k, RND);
            mpfr_fma(t[k], *_, u[k - j], t[k], RND);
        }
        mpfr_mul(s2[k], t[0], t[k], RND);
        for (int j = 1; j < k; j++) {
            mpfr_mul_d(*_, t[j], (k - j) / (double)k, RND);
            mpfr_fma(s2[k], *_, t[k - j], s2[k], RND);
        }
        mpfr_mul_2ui(s2[k], s2[k], 1, RND);
        if (g == HYP) {
            mpfr_neg(s2[k], s2[k], RND);
        }
    }
}


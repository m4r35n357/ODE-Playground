/*
 * Automatic Differentiation of Taylor Series, Newton's method
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"
#include "ad.h"

long n_max = 13, steps;

mpfr_t x0, x1, x_step, x_prev, f_prev, f_value, tmp, tmp1, D_1, D_05, D0, D1, D2, D3, D4, D5, D6, D7, *w1, *w2, *w3, *w4, *w5, *w6, *w7, *w_value, *_1, *_2, *_3, *wx, *wf, f_tol, x_tol, int_tol;

solver method;

model m;

void test_sqr (mpfr_t *f, mpfr_t *x, int n) {
    ad_square(f, x, n);
    ad_minus(f, f, w_value, n);
}

void trig (mpfr_t *f, mpfr_t *x, int n) {
    ad_tan_sec2(f, _1, x, n);
}

void cosx_x3 (mpfr_t *f, mpfr_t *x, int n) {
    //  Example: ./ad-test-newton-dbg 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    ad_square(_1, x, n);
    ad_product(_2, _1, x, n);
    ad_sin_cos(_1, _3, x, n);
    ad_minus(_3, _1, _2, n);
    ad_minus(f, _3, w_value, n);
}

void test_polynomial (mpfr_t *f, mpfr_t *x, int n) {
    ad_square(f, x, n);
    ad_product(_2, f, x, n);
    ad_scale(_1, x, D2, n);
    ad_minus(f, _2, _1, n);
    ad_minus(f, f, w5, n);
}

void septic (mpfr_t *f, mpfr_t *x, int n) {
    //  Example: ./ad-test-newton-dbg 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 50000 >/dev/null
    ad_minus(_2, x, w1, n);
    ad_plus(_1, x, w2, n);
    ad_product(_3, _2, _1, n);
    ad_minus(_1, x, w3, n);
    ad_product(_2, _3, _1, n);
    ad_plus(_1, x, w5, n);
    ad_product(_3, _2, _1, n);
    ad_minus(_1, x, w6, n);
    ad_product(_2, _3, _1, n);
    ad_plus(_1, x, w7, n);
    ad_product(_3, _2, _1, n);
    ad_product(f, _3, x, n);
}

void composite1 (mpfr_t *f, mpfr_t *x, int n) {
    //  Example: ./ad-test-newton-dbg 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    ad_square(_1, x, n);
    ad_minus(_2, _1, w4, n);
    ad_exp(_1, _2, n);
    ad_exp(_2, x, n);
    ad_plus(_3, _1, _2, n);
    ad_ln(f, _3, n);
}

void composite2 (mpfr_t *f, mpfr_t *x, int n) {
    //  Example: ./ad-test-newton-dbg 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    ad_exp(_1, x, n);
    ad_minus(_2, _1, w4, n);
    ad_square(_1, _2, n);
    ad_square(_2, x, n);
    ad_plus(_3, _1, _2, n);
    ad_sqrt(f, _3, n);
}

void lorentz (mpfr_t *f, mpfr_t *x, int n) {
    ad_square(_1, x, n);
    ad_minus(_2, w1, _1, n);
    ad_power(f, _2, D_05, n);
}

int main (int argc, char **argv) {
    assert(argc == 8);

    mpfr_set_default_prec(236);
    mpfr_inits(x_step, x_prev, f_prev, tmp, tmp1, NULL);

    method = strtol(argv[1], NULL, BASE);
    mpfr_init_set_str(x0, argv[2], BASE, RND);
    mpfr_init_set_str(x1, argv[3], BASE, RND);
    steps = strtol(argv[4], NULL, BASE);
    mpfr_init_set_str(f_value, argv[5], BASE, RND);
    mpfr_init_set_str(f_tol, argv[6], BASE, RND);
    mpfr_init_set_str(x_tol, argv[7], BASE, RND);

    mpfr_init_set_si(D_1, -1, RND);
    mpfr_init_set_d(D_05, -0.5, RND);
    mpfr_init_set_si(D0, 0, RND);
    mpfr_init_set_si(D1, 1, RND);
    mpfr_init_set_si(D2, 2, RND);
    mpfr_init_set_si(D3, 3, RND);
    mpfr_init_set_si(D4, 4, RND);
    mpfr_init_set_si(D5, 5, RND);
    mpfr_init_set_si(D6, 6, RND);
    mpfr_init_set_si(D7, 7, RND);
    w1 = t_jet_c(n_max, D1);
    w2 = t_jet_c(n_max, D2);
    w3 = t_jet_c(n_max, D3);
    w4 = t_jet_c(n_max, D4);
    w5 = t_jet_c(n_max, D5);
    w6 = t_jet_c(n_max, D6);
    w7 = t_jet_c(n_max, D7);

    w_value = t_jet_c(n_max, f_value);
    _1 = t_jet(n_max);
    _2 = t_jet(n_max);
    _3 = t_jet(n_max);

    wx = t_jet_c(n_max, x0);
    set_ad_status(wx, VARIABLE);
    wf = t_jet(n_max);

    m = septic;

    mpfr_sub(tmp, x1, x0, RND);
    mpfr_div_ui(x_step, tmp, steps, RND);
    for (int k = 0; k < steps + 1; k++) {
        mpfr_mul_ui(tmp, x_step, k, RND);
        mpfr_add(wx[0], x0, tmp, RND);
        m(wf, wx, n_max);
        mpfr_printf("%.3RNe ", wx[0]);
        jet_to_derivs(wf, n_max);
        for (int i = 0; i < n_max; i++) {
            mpfr_printf("%.6RNe ", wf[i]);
        }
        printf("\n");
        if (method != NONE) {
            if(k > 0) {
                mpfr_mul(tmp, f_prev, wf[0], RND);
                if (mpfr_sgn(tmp) < 0) {
                    fprintf(stderr, "Bracketed root, solving ");
                    if (method == BISECT) {
                        fprintf(stderr, "using bisection\n");
                        ad_bisect(m, t_jet_c(3, x_prev), wx, 100, f_tol, x_tol, t_jet(3), t_jet(3), t_jet(3));
                    } else if (method == NEWTON) {
                        fprintf(stderr, "using Newton's method\n");
                        ad_newton(m, t_jet_c(method + 2L, D1), wx, 100, f_tol, x_tol);
                    } else {
                        fprintf(stderr, "using Householder's method of degree %lu\n", method - 1L);
                        ad_householder(m, t_jet_c(method + 2L, D1), wx, method, 100, f_tol, x_tol, t_jet(method + 2L), t_jet_c(method + 2L, D1));
                    }
                    fprintf(stderr, "\n");
                }
            }
        }
        mpfr_set(x_prev, wx[0], RND);
        mpfr_set(f_prev, wf[0], RND);
    }

    return 0;
}

/*
 * Automatic Differentiation of Taylor Series, Newton's s
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"
#include "ad.h"

mpfr_t _, D_1, D_05, D0, D1, D2, D3, D4, D5, D6, D7, *w1, *w2, *w3, *w4, *w5, *w6, *w7, *target, *_1, *_2, *_3;

void test_sqr (mpfr_t *f, mpfr_t *x, int n) {
    ad_square(f, x, n);
    ad_minus(f, f, target, n);
}

void trig (mpfr_t *f, mpfr_t *x, int n) {
    ad_tan_sec2(f, _1, x, n, TRIG);
}

void cosx_x3 (mpfr_t *f, mpfr_t *x, int n) {
    //  Example: ./ad-test-newton-dbg 13 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    ad_square(_1, x, n);
    ad_product(_2, _1, x, n);
    ad_sin_cos(_1, _3, x, n, TRIG);
    ad_minus(_3, _1, _2, n);
    ad_minus(f, _3, target, n);
}

void test_polynomial (mpfr_t *f, mpfr_t *x, int n) {
    ad_square(f, x, n);
    ad_product(_2, f, x, n);
    ad_scale(_1, x, D2, n);
    ad_minus(f, _2, _1, n);
    ad_minus(f, f, w5, n);
}

void septic (mpfr_t *f, mpfr_t *x, int n) {
    //  Example: ./ad-test-newton-dbg 13 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 50000 >/dev/null
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
    //  Example: ./ad-test-newton-dbg 13 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    ad_square(_1, x, n);
    ad_minus(_2, _1, w4, n);
    ad_exp(_1, _2, n);
    ad_exp(_2, x, n);
    ad_plus(_3, _1, _2, n);
    ad_ln(f, _3, n);
}

void composite2 (mpfr_t *f, mpfr_t *x, int n) {
    //  Example: ./ad-test-newton-dbg 13 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
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
    ad_power(f, _2, - 0.5, n);
}

int main (int argc, char **argv) {
    mpfr_t x0, x1, x_step, x_prev, f_prev, f_target, f_tol, x_tol;

    assert(argc == 9);
    mpfr_set_default_prec(236);
    mpfr_inits(x_step, x_prev, f_prev, _, NULL);

    long order = strtol(argv[1], NULL, BASE);
    solver s = strtol(argv[2], NULL, BASE);
    mpfr_init_set_str(x0, argv[3], BASE, RND);
    mpfr_init_set_str(x1, argv[4], BASE, RND);
    long steps = strtol(argv[5], NULL, BASE);
    mpfr_init_set_str(f_target, argv[6], BASE, RND);
    mpfr_init_set_str(f_tol, argv[7], BASE, RND);
    mpfr_init_set_str(x_tol, argv[8], BASE, RND);

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

    w1 = t_jet_c(order, D1);
    w2 = t_jet_c(order, D2);
    w3 = t_jet_c(order, D3);
    w4 = t_jet_c(order, D4);
    w5 = t_jet_c(order, D5);
    w6 = t_jet_c(order, D6);
    w7 = t_jet_c(order, D7);

    target = t_jet_c(order, f_target);
    _1 = t_jet(order);
    _2 = t_jet(order);
    _3 = t_jet(order);

    model func = composite1;
    mpfr_t *f = t_jet(order);
    mpfr_t *x = t_jet_c(order, x0);
    set_ad_status(x, VARIABLE);

    mpfr_sub(_, x1, x0, RND);
    mpfr_div_ui(x_step, _, steps, RND);
    for (int k = 0; k < steps + 1; k++) {
        mpfr_mul_ui(_, x_step, k, RND);
        mpfr_add(x[0], x0, _, RND);
        func(f, x, order);
        mpfr_printf("%.3RNe ", x[0]);
        jet_to_derivs(f, order);
        for (int i = 0; i < order; i++) {
            mpfr_printf("%.6RNe ", f[i]);
        }
        printf("\n");
        if (s != NONE) {
            if(k > 0) {
                mpfr_mul(_, f_prev, f[0], RND);
                if (mpfr_sgn(_) < 0) {
                    fprintf(stderr, "Bracketed a root, using ");
                    if (s == BISECT) {
                        fprintf(stderr, "Bisection\n");
                        ad_bisect(func, t_jet_c(1, x_prev), x, 100, f_tol, x_tol, t_jet(1), t_jet(1), t_jet(1));
                    } else if (s == NEWTON) {
                        fprintf(stderr, "Newton's method\n");
                        ad_newton(func, t_jet_c(s, D1), x, 100, f_tol, x_tol);
                    } else {
                        fprintf(stderr, "Householder's method of degree %d\n", s - 1);
                        ad_householder(func, t_jet_c(s, D1), x, s, 100, f_tol, x_tol, t_jet(s), t_jet_c(s, D1));
                    }
                    fprintf(stderr, "\n");
                }
            }
        }
        mpfr_set(x_prev, x[0], RND);
        mpfr_set(f_prev, f[0], RND);
    }

    return 0;
}

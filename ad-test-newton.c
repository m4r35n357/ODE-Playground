/*
 * Automatic Differentiation of Taylor Series, Newton's method
 *
 * (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 *
 * Example: ./ad-test-newton-dbg 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 8 8000 &
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"
#include "ad.h"

long n_max = 7, n, steps;
mpfr_t x0, x1, x_step, x_prev, f_prev, f_value, tmp, D0, D1, D2, D3, D5, D6, D7, *w1, *w2, *w3, *w5, *w6, *w7, *w_value, *w_tmp, *w_tmp1, *w_tmp2, *w_tmp3, *wx, *wf, *wxa, *wxb, f_tol, x_tol, int_tol;
model m;

void test_sqr (mpfr_t *f, mpfr_t *x, int n) {
    ad_square(f, x, n, &tmp);
    ad_minus(f, f, w_value, n);
}

void cosx_x3 (mpfr_t *f, mpfr_t *x, int n) {
    ad_square(w_tmp1, x, n, &tmp);
    ad_product(w_tmp2, w_tmp1, x, n, &tmp);
    ad_sin_cos(w_tmp1, f, x, n, &tmp);
    ad_minus(f, w_tmp1, w_tmp2, n);
    ad_minus(f, f, w_value, n);
}

void test_polynomial (mpfr_t *f, mpfr_t *x, int n) {
    ad_square(f, x, n, &tmp);
    ad_product(w_tmp2, f, x, n, &tmp);
    ad_scale(w_tmp1, x, D2, n);
    ad_minus(f, w_tmp2, w_tmp1, n);
    ad_minus(f, f, w5, n);
}

void quintic (mpfr_t *f, mpfr_t *x, int n) {
    ad_minus(w_tmp2, x, w1, n);
    ad_plus(w_tmp1, x, w2, n);
    ad_product(w_tmp3, w_tmp2, w_tmp1, n, &tmp);
    ad_minus(w_tmp1, x, w3, n);
    ad_product(w_tmp2, w_tmp3, w_tmp1, n, &tmp);
    ad_plus(w_tmp1, x, w5, n);
    ad_product(w_tmp3, w_tmp2, w_tmp1, n, &tmp);
    ad_minus(w_tmp1, x, w6, n);
    ad_product(f, w_tmp3, w_tmp1, n, &tmp);
}

void sextic (mpfr_t *f, mpfr_t *x, int n) {
    ad_minus(w_tmp2, x, w1, n);
    ad_plus(w_tmp1, x, w2, n);
    ad_product(w_tmp3, w_tmp2, w_tmp1, n, &tmp);
    ad_minus(w_tmp1, x, w3, n);
    ad_product(w_tmp2, w_tmp3, w_tmp1, n, &tmp);
    ad_plus(w_tmp1, x, w5, n);
    ad_product(w_tmp3, w_tmp2, w_tmp1, n, &tmp);
    ad_minus(w_tmp1, x, w6, n);
    ad_product(w_tmp2, w_tmp3, w_tmp1, n, &tmp);
    ad_plus(w_tmp1, x, w7, n);
    ad_product(f, w_tmp2, w_tmp1, n, &tmp);
}

int main (int argc, char **argv) {
    assert(argc == 8);

    mpfr_set_default_prec(236);
    mpfr_inits(x_step, x_prev, f_prev, tmp, NULL);

    n = strtol(argv[1], NULL, BASE);
    assert(n > 0 && n <= n_max);
    mpfr_init_set_str(x0, argv[2], BASE, RND);
    mpfr_init_set_str(x1, argv[3], BASE, RND);
    steps = strtol(argv[4], NULL, BASE);
    mpfr_init_set_str(f_value, argv[5], BASE, RND);
    mpfr_init_set_str(f_tol, argv[6], BASE, RND);
    mpfr_init_set_str(x_tol, argv[7], BASE, RND);

    mpfr_init_set_si(D0, 0, RND);
    mpfr_init_set_si(D1, 1, RND);
    mpfr_init_set_si(D2, 2, RND);
    mpfr_init_set_si(D3, 3, RND);
    mpfr_init_set_si(D5, 5, RND);
    mpfr_init_set_si(D6, 6, RND);
    mpfr_init_set_si(D7, 7, RND);
    w1 = t_jet_constant(n_max, D1);
    w2 = t_jet_constant(n_max, D2);
    w3 = t_jet_constant(n_max, D3);
    w5 = t_jet_constant(n_max, D5);
    w6 = t_jet_constant(n_max, D6);
    w7 = t_jet_constant(n_max, D7);

    w_value = t_jet_constant(n_max, f_value);
    w_tmp1 = t_jet(n_max);
    w_tmp2 = t_jet(n_max);
    w_tmp3 = t_jet(n_max);

    wx = t_jet_constant(n_max, x0);
    mpfr_set_si(wx[1], 1, RND);
    wf = t_jet(n_max);

    m = sextic;

    mpfr_sub(tmp, x1, x0, RND);
    mpfr_div_si(x_step, tmp, steps, RND);
    for (int k = 0; k < steps + 1; k++) {
        mpfr_mul_si(tmp, x_step, k, RND);
        mpfr_add(wx[0], x0, tmp, RND);
        m(wf, wx, n_max);
        mpfr_printf("%.3RNe ", wx[0]);
        jet_to_derivs(wf, n_max);
        for (int i = 0; i < n_max; i++) {
            mpfr_printf("%.6RNe ", wf[i]);
        }
        printf("\n");
        if(k > 0) {
            mpfr_mul(tmp, f_prev, wf[0], RND);
            if (mpfr_sgn(tmp) < 0) {
                fprintf(stderr, "Bracketed root, solving . . .\n");
                if (n == 1) {
                    ad_bisect(m, t_jet_constant(n, x_prev), t_jet_constant(n, wx[0]), 100, f_tol, x_tol);
                } else if (n == 2) {
                    ad_newton(m, t_jet(n), t_jet_constant(n, wx[0]), 100, f_tol, x_tol);
                } else {
                    ad_householder(m, t_jet(n), t_jet_constant(n, wx[0]), n, 100, f_tol, x_tol);
                }
                fprintf(stderr, "\n");
            }
        }
        mpfr_set(x_prev, wx[0], RND);
        mpfr_set(f_prev, wf[0], RND);
    }

    return 0;
}

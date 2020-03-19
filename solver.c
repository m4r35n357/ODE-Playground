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

mpfr_t _, D_1, D_05, D0, D1, D2, D3, D4, D5, D6, D7, D54, D160, D641, D828, D1260;
mpfr_t *w1, *w2, *w3, *w4, *w5, *w6, *w7, *target, *_1, *_2, *_3, *__;

static void septic (mpfr_t *f, mpfr_t *x, int n) {
    //  Example: ./solver-dbg 0 13 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 50000 >/dev/null 2>&1
    // (x + 7)(x + 5)(x + 2)x(x - 1)(x - 3)(x - 6) = x^7 + 4x^6 - 54x^5 - 160x^4 + 641x^3 + 828x^2 - 1260x
    ad_minus(_2, x, w1, n);
    ad_plus(_1, x, w2, n);
    ad_prod(_3, _2, _1, n);
    ad_minus(_1, x, w3, n);
    ad_prod(_2, _3, _1, n);
    ad_plus(_1, x, w5, n);
    ad_prod(_3, _2, _1, n);
    ad_minus(_1, x, w6, n);
    ad_prod(_2, _3, _1, n);
    ad_plus(_1, x, w7, n);
    ad_prod(_3, _2, _1, n);
    ad_prod(f, _3, x, n);
}

static void septic2 (mpfr_t *f, mpfr_t *x, int n) {
    //  Example: ./solver-dbg 4 13 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 50000 >/dev/null 2>&1
    // x^7 + 4x^6 - 54x^5 - 160x^4 + 641x^3 + 828x^2 - 1260x = ((((((x + 4)x - 54)x - 160)x + 641)x + 828)x - 1260)x
    ad_set(f, x, n);
    mpfr_add(f[0], f[0], D4, RND);
    ad_prod(__, x, f, n);
    mpfr_sub(__[0], __[0], D54, RND);
    ad_prod(f, x, __, n);
    mpfr_sub(f[0], f[0], D160, RND);
    ad_prod(__, x, f, n);
    mpfr_add(__[0], __[0], D641, RND);
    ad_prod(f, x, __, n);
    mpfr_add(f[0], f[0], D828, RND);
    ad_prod(__, x, f, n);
    mpfr_sub(__[0], __[0], D1260, RND);
    ad_prod(f, x, __, n);
}

static void composite1 (mpfr_t *f, mpfr_t *x, int n) {
    //  Example: ./solver-dbg 1 13 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null 2>&1
    ad_sqr(_1, x, n);
    ad_minus(_2, _1, w4, n);
    ad_exp(_1, _2, n);
    ad_exp(_2, x, n);
    ad_plus(_3, _1, _2, n);
    ad_ln(f, _3, n);
}

static void composite2 (mpfr_t *f, mpfr_t *x, int n) {
    //  Example: ./solver-dbg 2 13 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null 2>&1
    ad_exp(_1, x, n);
    ad_minus(_2, _1, w4, n);
    ad_sqr(_1, _2, n);
    ad_sqr(_2, x, n);
    ad_plus(_3, _1, _2, n);
    ad_sqrt(f, _3, n);
}

static void cosx_x3 (mpfr_t *f, mpfr_t *x, int n) {
    //  Example: ./solver-dbg 3 13 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null 2>&1
    ad_sqr(_1, x, n);
    ad_prod(_2, _1, x, n);
    ad_sin_cos(_1, _3, x, n, TRIG);
    ad_minus(_3, _1, _2, n);
    ad_minus(f, _3, target, n);
}

int main (int argc, char **argv) {
    mpfr_t x0, x1, x_step, x_prev, f0_prev, f1_prev, f2_prev, f_target, f_tol, x_tol;
    model function;

    assert(argc == 10);
    mpfr_set_default_prec(236);
    ad_tempvars();
    mpfr_inits(x_step, x_prev, f0_prev, f1_prev, f2_prev, _, NULL);

    long model_code = strtol(argv[1], NULL, BASE);
    switch (model_code) {
        case 0 :
            function = septic;
            break;
        case 1 :
            function = composite1;
            break;
        case 2 :
            function = composite2;
            break;
        case 3 :
            function = cosx_x3;
            break;
        case 4 :
            function = septic2;
            break;
        default :
            fprintf(stderr, "Invalid model ID %ld, parameter 1 must be between %d and %d\n", model_code, 0, 4);
            return 1;
    }
    long order = strtol(argv[2], NULL, BASE);
    solver s = strtol(argv[3], NULL, BASE);
    switch (s) {
        case NONE :
        case NEWTON :
            fprintf(stderr, "\n");
            break;
        default :
            fprintf(stderr, "Invalid solver ID %d, parameter 3 must be %d or %d\n", s, NONE, NEWTON);
            return 1;
    }
    mpfr_init_set_str(x0, argv[4], BASE, RND);
    mpfr_init_set_str(x1, argv[5], BASE, RND);
    long steps = strtol(argv[6], NULL, BASE);
    mpfr_init_set_str(f_target, argv[7], BASE, RND);
    mpfr_init_set_str(f_tol, argv[8], BASE, RND);
    mpfr_init_set_str(x_tol, argv[9], BASE, RND);

    mpfr_init_set_si(D_1, -1, RND);
    mpfr_init_set_d(D_05, -0.5, RND);
    mpfr_init_set_ui(D0, 0, RND);
    mpfr_init_set_ui(D1, 1, RND);
    mpfr_init_set_ui(D2, 2, RND);
    mpfr_init_set_ui(D3, 3, RND);
    mpfr_init_set_ui(D4, 4, RND);
    mpfr_init_set_ui(D5, 5, RND);
    mpfr_init_set_ui(D6, 6, RND);
    mpfr_init_set_ui(D7, 7, RND);
    mpfr_init_set_ui(D54, 54, RND);
    mpfr_init_set_ui(D160, 160, RND);
    mpfr_init_set_ui(D641, 641, RND);
    mpfr_init_set_ui(D828, 828, RND);
    mpfr_init_set_ui(D1260, 1260, RND);

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
    __ = t_jet(order);

    mpfr_t *f = t_jet_c(order, D0);
    mpfr_t *f_ = t_jet_c(NEWTON + INFLECTION, D0);
    mpfr_t *x = t_jet_c(order, x0);
    set_ad_status(x, VARIABLE);

    mpfr_sub(_, x1, x0, RND);
    mpfr_div_ui(x_step, _, steps, RND);
    for (int k = 0; k < steps + 1; k++) {
        mpfr_mul_ui(_, x_step, k, RND);
        mpfr_add(x[0], x0, _, RND);
        function(f, x, order);
        mpfr_printf("%.3RNe ", x[0]);
        jet_to_derivs(f, order);
        for (int i = 0; i < order; i++) {
            mpfr_printf("%.6RNe ", f[i]);
        }
        printf("\n");
        if (s != NONE) {
            if(k > 0) {
                mpfr_mul(_, f0_prev, f[ROOT], RND);
                if (mpfr_sgn(_) < 0) {
                    ad_newton(function, f_, x, 100, f_tol, x_tol, ROOT);
                    mpfr_cmp(f0_prev, f[ROOT]) > 0 ? fprintf(stderr, "\\") : fprintf(stderr, "/");
                    fprintf(stderr, " ROOT\n");
                }
                mpfr_mul(_, f1_prev, f[MIN_MAX], RND);
                if (mpfr_sgn(_) < 0) {
                    ad_newton(function, f_, x, 100, f_tol, x_tol, MIN_MAX);
                    mpfr_cmp(f1_prev, f[MIN_MAX]) > 0 ? fprintf(stderr, "\\ MAXIMUM\n") : fprintf(stderr, "/ MINIMUM\n");
                }
                mpfr_mul(_, f2_prev, f[INFLECTION], RND);
                if (mpfr_sgn(_) < 0) {
                    ad_newton(function, f_, x, 100, f_tol, x_tol, INFLECTION);
                    mpfr_cmp(f2_prev, f[INFLECTION]) > 0 ? fprintf(stderr, "\\") : fprintf(stderr, "/");
                    fprintf(stderr, " INFLECTION\n");
                }
            }
        }
        mpfr_set(x_prev, x[0], RND);
        mpfr_set(f0_prev, f[ROOT], RND);
        mpfr_set(f1_prev, f[MIN_MAX], RND);
        mpfr_set(f2_prev, f[INFLECTION], RND);
    }
    fprintf(stderr, "\n");
    return 0;
}

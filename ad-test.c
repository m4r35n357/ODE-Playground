/*
 * Automatic Differentiation of Taylor Series, validation checks
 *
 * Example: ./ad-test-dbg 7 2 1
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"
#include "ad.h"

#define KNRM "\x1B[0;37m"
#define KCYN "\x1B[36m"
#define KGRY "\x1B[2;37m"
#define KBLD "\x1B[1;37m"

long order, n;
mpfr_t x, y, _, D_1, D_05, D0, D05, D1, D2, D3, D5, D6, D7, *cx, *cy, *cx0, *c1, *c2, *c3, *c5, *c6, *c7, *PI_3, *PI_4, *we, *wl, *ws, *wc, *wt, *ws2, *wsqr, *wsqrt, *wsum, *wprod, *wquot, *wpwr, *__, *_1, *_2, *_3;

void septic (mpfr_t *f, const mpfr_t *x, int n) {
    ad_minus(_2, x, c1, n);
    ad_plus(_1, x, c2, n);
    ad_product(_3, _2, _1, n);
    ad_minus(_1, x, c3, n);
    ad_product(_2, _3, _1, n);
    ad_plus(_1, x, c5, n);
    ad_product(_3, _2, _1, n);
    ad_minus(_1, x, c6, n);
    ad_product(_2, _3, _1, n);
    ad_plus(_1, x, c7, n);
    ad_product(_3, _2, _1, n);
    ad_product(f, _3, x, n);
}

int main (int argc, char **argv) {
    assert(argc == 4);

    mpfr_set_default_prec(113);

    mpfr_inits(_, NULL);
    n = strtol(argv[1], NULL, BASE);
    assert(n > 1);
    mpfr_init_set_str(x, argv[2], BASE, RND);
    mpfr_init_set_str(y, argv[3], BASE, RND);

    mpfr_init_set_si(D_1, -1, RND);
    mpfr_init_set_d(D_05, -0.5, RND);
    mpfr_init_set_ui(D0, 0, RND);
    mpfr_init_set_d(D05, 0.5, RND);
    mpfr_init_set_ui(D1, 1, RND);
    mpfr_init_set_ui(D2, 2, RND);
    mpfr_init_set_ui(D3, 3, RND);
    mpfr_init_set_ui(D5, 5, RND);
    mpfr_init_set_ui(D6, 6, RND);
    mpfr_init_set_ui(D7, 7, RND);

    c1 = t_jet_c(n, D1);
    c2 = t_jet_c(n, D2);
    c3 = t_jet_c(n, D3);
    c5 = t_jet_c(n, D5);
    c6 = t_jet_c(n, D6);
    c7 = t_jet_c(n, D7);

    cx = t_jet_c(n, x);
    cx0 = t_jet_c(n, D0);
    cy = t_jet_c(n, y);
    set_ad_status(cx, VARIABLE);
    set_ad_status(cx0, VARIABLE);
    set_ad_status(cy, CONSTANT);

    mpfr_const_pi(_, RND);
    mpfr_div_ui(_, _, 3, RND);
    PI_3 = t_jet_c(n, _);
    set_ad_status(PI_3, VARIABLE);

    mpfr_const_pi(_, RND);
    mpfr_div_ui(_, _, 4, RND);
    PI_4 = t_jet_c(n, _);
    set_ad_status(PI_4, VARIABLE);

    wsqr = t_jet(n);
    wsqrt = t_jet(n);
    wsum = t_jet(n);
    wprod = t_jet(n);
    wquot = t_jet(n);
    wpwr = t_jet(n);
    we = t_jet(n);
    wl = t_jet(n);
    ws = t_jet(n);
    wc = t_jet(n);
    wt = t_jet(n);
    ws2 = t_jet(n);
    __ = t_jet(n);
    _1 = t_jet(n);
    _2 = t_jet(n);
    _3 = t_jet(n);

    printf("\n%sx = %s, y = %s, order = %ld%s\n\n", KBLD, argv[2], argv[3], n - 1, KNRM);

//    printf("%s%s%s\n", KCYN, "Horner", KNRM);
//    t_horner(&x, (int[]){-19, 7, -4, 6}, n, D3);
//    t_horner(&x, (float[]){-19.0, 7.0, -4.0, 6.0}, n, D3);
//    mpfr_printf("int: %9.6RNf, float: %9.6RNf\n", *x, *y);

    printf("%s%s%s\n", KCYN, "f(x) = x^2", KNRM);
    ad_square(wsqr, cx, n);
    jet_output(wsqr, n, KNRM, KGRY);
    ad_product(wprod, cx, cx, n);
    jet_output(wprod, n, KNRM, KGRY);
    derivative_output(wsqr, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^4", KNRM);
    ad_square(__, cx, n);
    ad_square(wsqr, __, n);
    jet_output(wsqr, n, KNRM, KGRY);
    ad_product(__, cx, cx, n);
    ad_product(wprod, __, __, n);
    jet_output(wprod, n, KNRM, KGRY);
    derivative_output(wsqr, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = sqrt(x)", KNRM);
    ad_power(wpwr, cx, D05, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    ad_sqrt(wsqrt, cx, n);
    jet_output(wsqrt, n, KNRM, KGRY);
    derivative_output(wsqrt, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = sqrt(x^2)", KNRM);
    ad_square(wsqr, cx, n);
    ad_sqrt(wsqrt, wsqr, n);
    jet_output(wsqrt, n, KNRM, KGRY);
    derivative_output(wsqrt, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = 1 / x", KNRM);
    ad_power(wpwr, cx, D_1, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    ad_quotient(wquot, c1, cx, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x / x", KNRM);
    ad_scale(__, cx, D1, n);
    ad_quotient(wquot, cx, __, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = 1 / septic(x)", KNRM);
    septic(__, cx, n);
    ad_power(wpwr, __, D_1, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    ad_quotient(wquot, c1, __, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^-0.5", KNRM);
    ad_power(wpwr, cx, D_05, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    ad_sqrt(wsqrt, cx, n);
    ad_quotient(wquot, c1, wsqrt, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = e^x", KNRM);
    ad_exp(we, cx, n);
    jet_output(we, n, KNRM, KGRY);
    derivative_output(we, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = log(x)", KNRM);
    ad_ln(wl, cx, n);
    jet_output(wl, n, KNRM, KGRY);
    derivative_output(wl, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = log(e^x)", KNRM);
    ad_exp(we, cx, n);
    ad_ln(wl, we, n);
    jet_output(wl, n, KNRM, KGRY);
    derivative_output(wl, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = sin(0), f(x) = cos(0)", KNRM);
    ad_sin_cos(ws, wc, cx0, n);
    jet_output(ws, n, KNRM, KGRY);
    derivative_output(ws, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(wc, n, KNRM, KGRY);
    derivative_output(wc, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = sin(pi/3), f(x) = cos(pi/3)", KNRM);
    ad_sin_cos(ws, wc, PI_3, n);
    jet_output(ws, n, KNRM, KGRY);
    derivative_output(ws, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(wc, n, KNRM, KGRY);
    derivative_output(wc, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = sinh(pi/3), f(x) = cosh(pi/3)", KNRM);
    ad_sinh_cosh(ws, wc, PI_3, n);
    jet_output(ws, n, KNRM, KGRY);
    derivative_output(ws, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(wc, n, KNRM, KGRY);
    derivative_output(wc, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = tan(pi/4), f(x) = sec(pi/4)^2", KNRM);
    ad_tan_sec2(wt, ws2, PI_4, n);
    jet_output(wt, n, KNRM, KGRY);
    derivative_output(wt, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(ws2, n, KNRM, KGRY);
    derivative_output(ws2, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = tanh(pi/4), f(x) = sech(pi/4)^2", KNRM);
    ad_tanh_sech2(wt, ws2, PI_4, n);
    jet_output(wt, n, KNRM, KGRY);
    derivative_output(wt, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(ws2, n, KNRM, KGRY);
    derivative_output(ws2, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x, y) = x * y, d/dx", KNRM);
    ad_product(wprod, cx, cy, n);
    jet_output(wprod, n, KNRM, KGRY);
    derivative_output(wprod, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x, y) = e^x / (x - y * sin(x^2), d/dx", KNRM);
    ad_square(wsqr, cx, n);
    ad_sin_cos(ws, wc, wsqr, n);
    ad_product(wprod, cy, ws, n);
    ad_minus(wsum, cx, wprod, n);
    ad_exp(we, cx, n);
    ad_quotient(wquot, we, wsum, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    set_ad_status(cx, CONSTANT);
    set_ad_status(cy, VARIABLE);

    printf("%s%s%s\n", KCYN, "f(x, y) = x * y, d/dy", KNRM);
    ad_product(wprod, cx, cy, n);
    jet_output(wprod, n, KNRM, KGRY);
    derivative_output(wprod, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x, y) = e^x / (x - y * sin(x^2), d/dy", KNRM);
    ad_square(wsqr, cx, n);
    ad_sin_cos(ws, wc, wsqr, n);
    ad_product(wprod, cy, ws, n);
    ad_minus(wsum, cx, wprod, n);
    ad_exp(we, cx, n);
    ad_quotient(wquot, we, wsum, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    return 0;
}

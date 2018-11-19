/*
 * Automatic Differentiation of Taylor Series, validation checks
 *
 * Example: ./ad-test-dbg 7 2 1
 *
 * (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
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
mpfr_t x, y, _, D_1, D_05, D0, D1, D2, D3, D5, D6, D7, D9, *cx, *cy, *cx0, *cx1, *c1, *c2, *c3, *c5, *c6, *c7, *cx9, *PI_3, *PI_4, *we, *wl, *ws, *wc, *wt, *ws2, *wsq, *wsum, *wprod, *wnum, *wdenom, *wquot, *wpwr, *__, *_1, *_2, *_3;

void septic (mpfr_t *f, mpfr_t *x, int n) {
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
    mpfr_init_set_ui(D1, 1, RND);
    mpfr_init_set_ui(D2, 2, RND);
    mpfr_init_set_ui(D3, 3, RND);
    mpfr_init_set_ui(D5, 5, RND);
    mpfr_init_set_ui(D6, 6, RND);
    mpfr_init_set_ui(D7, 7, RND);
    mpfr_init_set_ui(D9, 9, RND);

    cx = t_jet_constant(n, x);
    cx0 = t_jet_constant(n, D0);
    cx1 = t_jet_constant(n, D1);
    cy = t_jet_constant(n, y);

    set_ad_status(cx, VARIABLE);
    set_ad_status(cx0, VARIABLE);
    set_ad_status(cx1, VARIABLE);
    set_ad_status(cy, CONSTANT);

    c1 = t_jet_constant(n, D1);
    c2 = t_jet_constant(n, D2);
    c3 = t_jet_constant(n, D3);
    c5 = t_jet_constant(n, D5);
    c6 = t_jet_constant(n, D6);
    c7 = t_jet_constant(n, D7);
    cx9 = t_jet_constant(n, D9);
    set_ad_status(cx9, VARIABLE);

    mpfr_const_pi(_, RND);
    mpfr_div_ui(_, _, 3, RND);
    PI_3 = t_jet_constant(n, _);
    set_ad_status(PI_3, VARIABLE);

    mpfr_const_pi(_, RND);
    mpfr_div_ui(_, _, 4, RND);
    PI_4 = t_jet_constant(n, _);
    set_ad_status(PI_4, VARIABLE);

    wsq = t_jet(n);
    wsum = t_jet(n);
    wprod = t_jet(n);
    wnum = t_jet(n);
    wdenom = t_jet(n);
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

    printf("%s%s%s\n", KCYN, "f(x) = x", KNRM);
    jet_output(cx, n, KNRM, KGRY);
    derivative_output(cx, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^2", KNRM);
    ad_product(wprod, cx, cx, n);
    jet_output(wprod, n, KNRM, KGRY);
    ad_square(wsq, cx, n);
    jet_output(wsq, n, KNRM, KGRY);
    derivative_output(wsq, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^3", KNRM);
    ad_square(wsq, cx, n);
    ad_product(wprod, wsq, cx, n);
    jet_output(wprod, n, KNRM, KGRY);
    derivative_output(wprod, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^4", KNRM);
    ad_square(wsq, cx, n);
    ad_product(wprod, wsq, wsq, n);
    jet_output(wprod, n, KNRM, KGRY);
    ad_product(wprod, cx, cx, n);
    ad_square(wsq, wprod, n);
    jet_output(wsq, n, KNRM, KGRY);
    derivative_output(wsq, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = 1 / x", KNRM);
    ad_power(wpwr, cx, D_1, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    ad_quotient(wquot, c1, cx, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = 1 / (x + 3)", KNRM);
    ad_plus(wdenom, cx, c3, n);
    ad_power(wpwr, wdenom, D_1, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    ad_quotient(wquot, c1, wdenom, n);
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

    printf("%s%s%s\n", KCYN, "f(x) = (2x + 3) / (x + 2)", KNRM);
    ad_scale(__, cx, D2, n);
    ad_plus(wnum, __, c3, n);
    ad_plus(wdenom, cx, c2, n);
    ad_quotient(wquot, wnum, wdenom, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^a, x = 9, a = -0.5", KNRM);
    ad_power(wpwr, cx9, D_05, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = e^x", KNRM);
    ad_exp(we, cx1, n);
    jet_output(we, n, KNRM, KGRY);
    derivative_output(we, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = log(x)", KNRM);
    ad_ln(wl, cx, n);
    jet_output(wl, n, KNRM, KGRY);
    derivative_output(wl, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = sin(x), f(x) = cos(x), x = 0", KNRM);
    ad_sin_cos(ws, wc, cx0, n);
    jet_output(ws, n, KNRM, KGRY);
    derivative_output(ws, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(wc, n, KNRM, KGRY);
    derivative_output(wc, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = sin(x), f(x) = cos(x), x = pi / 3", KNRM);
    ad_sin_cos(ws, wc, PI_3, n);
    jet_output(ws, n, KNRM, KGRY);
    derivative_output(ws, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(wc, n, KNRM, KGRY);
    derivative_output(wc, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = sinh(x), f(x) = cosh(x), x = pi / 3", KNRM);
    ad_sinh_cosh(ws, wc, PI_3, n);
    jet_output(ws, n, KNRM, KGRY);
    derivative_output(ws, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(wc, n, KNRM, KGRY);
    derivative_output(wc, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = tan(x), f(x) = sec(x)^2, x = pi / 4", KNRM);
    ad_tan_sec2(wt, ws2, PI_4, n);
    jet_output(wt, n, KNRM, KGRY);
    derivative_output(wt, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(ws2, n, KNRM, KGRY);
    derivative_output(ws2, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = tanh(x), f(x) = sech(x)^2, x = pi / 4", KNRM);
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
    ad_square(wsq, cx, n);
    ad_sin_cos(ws, wc, wsq, n);
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
    ad_square(wsq, cx, n);
    ad_sin_cos(ws, wc, wsq, n);
    ad_product(wprod, cy, ws, n);
    ad_minus(wsum, cx, wprod, n);
    ad_exp(we, cx, n);
    ad_quotient(wquot, we, wsum, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    return 0;
}

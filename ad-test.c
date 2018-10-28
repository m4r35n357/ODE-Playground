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
mpfr_t x, y, tmp, D_1, D_05, D0, D1, D2, D3, D9, *cx, *cy, *cx0, *cx1, *c1, *c2, *c3, *cx9, *PI_3, *PI_4, *we, *wl, *ws, *wc, *wt, *ws2,
        *wsq, *wsum, *wprod, *wnum, *wdenom, *wquot, *wpwr, *wtmp;

int main (int argc, char **argv) {
    assert(argc == 4);

    mpfr_set_default_prec(113);

    mpfr_inits(tmp, NULL);
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
    mpfr_init_set_ui(D9, 9, RND);

    cx = t_jet_constant(n, x);
    cy = t_jet_constant(n, y);

    cx0 = t_jet_constant(n, D0);
    mpfr_set_ui(cx0[1], 1, RND);
    cx1 = t_jet_constant(n, D1);
    mpfr_set_ui(cx1[1], 1, RND);

    c1 = t_jet_constant(n, D1);
    c2 = t_jet_constant(n, D2);
    c3 = t_jet_constant(n, D3);

    cx9 = t_jet_constant(n, D9);
    mpfr_set_ui(cx9[1], 1, RND);

    mpfr_const_pi(tmp, RND);
    mpfr_div_ui(tmp, tmp, 3, RND);
    PI_3 = t_jet_constant(n, tmp);
    mpfr_set_ui(PI_3[1], 1, RND);

    mpfr_const_pi(tmp, RND);
    mpfr_div_ui(tmp, tmp, 4, RND);
    PI_4 = t_jet_constant(n, tmp);
    mpfr_set_ui(PI_4[1], 1, RND);

    wtmp = t_jet(n);
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

    printf("\n%sx = %s, y = %s, order = %ld%s\n\n", KBLD, argv[2], argv[3], n - 1, KNRM);

    mpfr_set_ui(cx[1], 1, RND);
    mpfr_set_ui(cy[1], 0, RND);

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

    printf("%s%s%s\n", KCYN, "f(x) = 1 / (x + 3)", KNRM);
    ad_plus(wdenom, cx, c3, n);
    ad_power(wpwr, wdenom, D_1, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    ad_quotient(wquot, c1, wdenom, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = (2x + 3) / (x + 2)", KNRM);
    ad_scale(wtmp, cx, D2, n);
    ad_plus(wnum, wtmp, c3, n);
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
    ad_exp(we, cx1, n, &tmp);
    jet_output(we, n, KNRM, KGRY);
    derivative_output(we, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = sin(x), f(x) = cos(x), x = 0", KNRM);
    ad_sin_cos(ws, wc, cx0, n, &tmp);
    jet_output(ws, n, KNRM, KGRY);
    derivative_output(ws, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(wc, n, KNRM, KGRY);
    derivative_output(wc, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = sin(x), f(x) = cos(x), x = pi / 3", KNRM);
    ad_sin_cos(ws, wc, PI_3, n, &tmp);
    jet_output(ws, n, KNRM, KGRY);
    derivative_output(ws, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(wc, n, KNRM, KGRY);
    derivative_output(wc, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = sinh(x), f(x) = cosh(x), x = pi / 3", KNRM);
    ad_sinh_cosh(ws, wc, PI_3, n, &tmp);
    jet_output(ws, n, KNRM, KGRY);
    derivative_output(ws, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(wc, n, KNRM, KGRY);
    derivative_output(wc, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = tan(x), f(x) = sec(x)^2, x = pi / 4", KNRM);
    ad_tan_sec2(wt, ws2, PI_4, n, &tmp);
    jet_output(wt, n, KNRM, KGRY);
    derivative_output(wt, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(ws2, n, KNRM, KGRY);
    derivative_output(ws2, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = tanh(x), f(x) = sech(x)^2, x = pi / 4", KNRM);
    ad_tanh_sech2(wt, ws2, PI_4, n, &tmp);
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
    ad_sin_cos(ws, wc, wsq, n, &tmp);
    ad_product(wprod, cy, ws, n);
    ad_minus(wsum, cx, wprod, n);
    ad_exp(we, cx, n, &tmp);
    ad_quotient(wquot, we, wsum, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    mpfr_set_ui(cx[1], 0, RND);
    mpfr_set_ui(cy[1], 1, RND);

    printf("%s%s%s\n", KCYN, "f(x, y) = x * y, d/dy", KNRM);
    ad_product(wprod, cx, cy, n);
    jet_output(wprod, n, KNRM, KGRY);
    derivative_output(wprod, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x, y) = e^x / (x - y * sin(x^2), d/dy", KNRM);
    ad_square(wsq, cx, n);
    ad_sin_cos(ws, wc, wsq, n, &tmp);
    ad_product(wprod, cy, ws, n);
    ad_minus(wsum, cx, wprod, n);
    ad_exp(we, cx, n, &tmp);
    ad_quotient(wquot, we, wsum, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    return 0;
}

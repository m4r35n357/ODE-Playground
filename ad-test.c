/*
 * Automatic Differentiation of Taylor Series, validation checks
 *
 * Example: ./ad-test-dbg 7 2 1 >/tmp/ad-test.txt; diff --context=1 /tmp/ad-test.txt ad-test.txt
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
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

static mpfr_t *___;

static void septic (mpfr_t *f, mpfr_t *value, int order) {
    ad_set(f, value, order);
    mpfr_add_ui(f[0], f[0], 4, RND);
    ad_prod(___, value, f, order);
    mpfr_sub_ui(___[0], ___[0], 54, RND);
    ad_prod(f, value, ___, order);
    mpfr_sub_ui(f[0], f[0], 160, RND);
    ad_prod(___, value, f, order);
    mpfr_add_ui(___[0], ___[0], 641, RND);
    ad_prod(f, value, ___, order);
    mpfr_add_ui(f[0], f[0], 828, RND);
    ad_prod(___, value, f, order);
    mpfr_sub_ui(___[0], ___[0], 1260, RND);
    ad_prod(f, value, ___, order);
}

int main (int argc, char **argv) {
    long n;
    mpfr_t x, y, _, D0, D1, D3;

    assert(argc == 4);
    mpfr_set_default_prec(113);
    ad_tempvars();
    mpfr_inits(_, NULL);
    n = strtol(argv[1], NULL, BASE);
    assert(n > 1);
    mpfr_init_set_str(x, argv[2], BASE, RND);
    mpfr_init_set_str(y, argv[3], BASE, RND);

    mpfr_init_set_ui(D0, 0, RND);
    mpfr_init_set_ui(D1, 1, RND);
    mpfr_init_set_ui(D3, 3, RND);

    ___ = t_jet(n);

    mpfr_t *c1 = t_jet_c(n, D1);

    mpfr_t *cx = t_jet_c(n, x);
    mpfr_t *cx0 = t_jet_c(n, D0);
    mpfr_t *cy = t_jet_c(n, y);
    set_ad_status(cx, VARIABLE);
    set_ad_status(cx0, VARIABLE);
    set_ad_status(cy, CONSTANT);

    mpfr_const_pi(_, RND);
    mpfr_div_ui(_, _, 3, RND);
    mpfr_t *PI_3 = t_jet_c(n, _);
    set_ad_status(PI_3, VARIABLE);

    mpfr_const_pi(_, RND);
    mpfr_div_ui(_, _, 4, RND);
    mpfr_t *PI_4 = t_jet_c(n, _);
    set_ad_status(PI_4, VARIABLE);

    mpfr_t *wsqr = t_jet(n);
    mpfr_t *wabs = t_jet(n);
    mpfr_t *wsqrt = t_jet(n);
    mpfr_t *wsum = t_jet(n);
    mpfr_t *wprod = t_jet(n);
    mpfr_t *wquot = t_jet(n);
    mpfr_t *wpwr = t_jet(n);
    mpfr_t *we = t_jet(n);
    mpfr_t *wl = t_jet(n);
    mpfr_t *ws = t_jet(n);
    mpfr_t *wc = t_jet(n);
    mpfr_t *wt = t_jet(n);
    mpfr_t *ws2 = t_jet(n);
    mpfr_t *__ = t_jet_c(n + 1, D0);

    printf("%s%s%s\n", KCYN, "Horner", KNRM);
    mpfr_set_si(__[0], -19, RND);
    mpfr_set_ui(__[1], 7, RND);
    mpfr_set_si(__[2], -4, RND);
    mpfr_set_ui(__[3], 6, RND);
    t_horner(__, n, D3);
    mpfr_printf(" %7.3RNf\n", __[0]);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "x", KNRM);
    derivative_output(cx, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "y", KNRM);
    derivative_output(cy, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^2.0", KNRM);
    ad_pwr(wpwr, cx, 2.0, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = sqr(x)", KNRM);
    ad_sqr(wsqr, cx, n);
    jet_output(wsqr, n, KNRM, KGRY);
    derivative_output(wsqr, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = x * x", KNRM);
    ad_prod(wprod, cx, cx, n);
    jet_output(wprod, n, KNRM, KGRY);
    derivative_output(wprod, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^4.0", KNRM);
    ad_pwr(wpwr, cx, 4.0, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = sqr(sqr(x))", KNRM);
    ad_sqr(__, cx, n);
    ad_sqr(wsqr, __, n);
    jet_output(wsqr, n, KNRM, KGRY);
    derivative_output(wsqr, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = x * x * x * x", KNRM);
    ad_prod(__, cx, cx, n);
    ad_prod(wprod, __, __, n);
    jet_output(wprod, n, KNRM, KGRY);
    derivative_output(wprod, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^0.5", KNRM);
    ad_pwr(wpwr, cx, 0.5, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = sqrt(x)", KNRM);
    ad_sqrt(wsqrt, cx, n);
    jet_output(wsqrt, n, KNRM, KGRY);
    derivative_output(wsqrt, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^-0.5", KNRM);
    ad_pwr(wpwr, cx, - 0.5, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = 1 / sqrt(x)", KNRM);
    ad_sqrt(wsqrt, cx, n);
    ad_quot(wquot, c1, wsqrt, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^-1", KNRM);
    ad_pwr(wpwr, cx, - 1.0, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = 1 / x", KNRM);
    ad_quot(wquot, c1, cx, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^0", KNRM);
    ad_pwr(wpwr, cx, 0.0, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = x / x", KNRM);
    ad_scale(__, cx, D1, n);
    ad_quot(wquot, cx, __, n);
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

    printf("%s%s%s\n", KCYN, "f(x) = sqrt(x * x)", KNRM);
    ad_sqrt(wsqrt, ad_prod(wprod, cx, cx, n), n);
    jet_output(wsqrt, n, KNRM, KGRY);
    derivative_output(wsqrt, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = sqrt(x * x) - |x|", KNRM);
    ad_minus(__, wsqrt, ad_abs(wabs, cx, n), n);
    jet_output(__, n, KNRM, KGRY);
    derivative_output(__, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = log(e^x)", KNRM);
    ad_exp(we, cx, n);
    ad_ln(wl, we, n);
    jet_output(wl, n, KNRM, KGRY);
    derivative_output(wl, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    ad_sin_cos(ws, wc, cx0, n, TRIG);
    printf("%s%s%s\n", KCYN, "f(x) = sin(0)", KNRM);
    jet_output(ws, n, KNRM, KGRY);
    derivative_output(ws, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = cos(0)", KNRM);
    jet_output(wc, n, KNRM, KGRY);
    derivative_output(wc, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    ad_sin_cos(ws, wc, PI_3, n, TRIG);
    printf("%s%s%s\n", KCYN, "f(x) = sin(pi/3)", KNRM);
    jet_output(ws, n, KNRM, KGRY);
    derivative_output(ws, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = cos(pi/3)", KNRM);
    jet_output(wc, n, KNRM, KGRY);
    derivative_output(wc, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    ad_sin_cos(ws, wc, PI_3, n, HYP);
    printf("%s%s%s\n", KCYN, "f(x) = sinh(pi/3)", KNRM);
    jet_output(ws, n, KNRM, KGRY);
    derivative_output(ws, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = cosh(pi/3)", KNRM);
    jet_output(wc, n, KNRM, KGRY);
    derivative_output(wc, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    ad_tan_sec2(wt, ws2, PI_4, n, TRIG);
    printf("%s%s%s\n", KCYN, "f(x) = tan(pi/4)", KNRM);
    jet_output(wt, n, KNRM, KGRY);
    derivative_output(wt, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = sec(pi/4)^2", KNRM);
    jet_output(ws2, n, KNRM, KGRY);
    derivative_output(ws2, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    ad_tan_sec2(wt, ws2, PI_4, n, HYP);
    printf("%s%s%s\n", KCYN, "f(x) = tanh(pi/4)", KNRM);
    jet_output(wt, n, KNRM, KGRY);
    derivative_output(wt, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = sech(pi/4)^2", KNRM);
    jet_output(ws2, n, KNRM, KGRY);
    derivative_output(ws2, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x, y) = x * y, d/dx", KNRM);
    ad_prod(wprod, cx, cy, n);
    jet_output(wprod, n, KNRM, KGRY);
    derivative_output(wprod, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x, y) = sqr(x^2 - y), d/dx", KNRM);
    ad_sqr(wsqr, cx, n);
    ad_minus(wsum, wsqr, cy, n);
    ad_sqr(wsqr, wsum, n);
    jet_output(wsqr, n, KNRM, KGRY);
    derivative_output(wsqr, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x, y) = e^x / (x - y * sin(x^2), d/dx", KNRM);
    ad_sqr(wsqr, cx, n);
    ad_sin_cos(ws, wc, wsqr, n, TRIG);
    ad_prod(wprod, cy, ws, n);
    ad_minus(wsum, cx, wprod, n);
    ad_exp(we, cx, n);
    ad_quot(wquot, we, wsum, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = septic(x)^-1", KNRM);
    septic(__, cx, n);
    ad_pwr(wpwr, __, - 1.0, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = 1 / septic(x)", KNRM);
    ad_quot(wquot, c1, __, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    set_ad_status(cx, CONSTANT);
    set_ad_status(cy, VARIABLE);

    printf("%s%s%s\n", KCYN, "f(x, y) = x * y, d/dy", KNRM);
    ad_prod(wprod, cx, cy, n);
    jet_output(wprod, n, KNRM, KGRY);
    derivative_output(wprod, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x, y) = sqr(x^2 - y), d/dy", KNRM);
    ad_sqr(wsqr, cx, n);
    ad_minus(wsum, wsqr, cy, n);
    ad_sqr(wsqr, wsum, n);
    jet_output(wsqr, n, KNRM, KGRY);
    derivative_output(wsqr, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x, y) = e^x / (x - y * sin(x^2), d/dy", KNRM);
    ad_sqr(wsqr, cx, n);
    ad_sin_cos(ws, wc, wsqr, n, TRIG);
    ad_prod(wprod, cy, ws, n);
    ad_minus(wsum, cx, wprod, n);
    ad_exp(we, cx, n);
    ad_quot(wquot, we, wsum, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    return 0;
}

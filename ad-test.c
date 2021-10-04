/*
 * Automatic Differentiation of Taylor Series, legacy validation checks
 *
 * Example: ./ad-test-dbg 6 2 1 >/tmp/ad-test.txt; diff --context=1 /tmp/ad-test.txt ad-test.txt
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
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

static series ___;

static mpfr_t D0, D05, D_05, D1, D_1, D2, D_2, D3, D4;

static void ad_test_tempvars (void) {
    ad_tempvars();
    mpfr_init_set_ui(D0, 0, RND);
    mpfr_init_set_str(D05, "0.5", BASE, RND);
    mpfr_init_set_str(D_05, "-0.5", BASE, RND);
    mpfr_init_set_ui(D1, 1, RND);
    mpfr_init_set_si(D_1, -1, RND);
    mpfr_init_set_ui(D2, 2, RND);
    mpfr_init_set_si(D_2, -2, RND);
    mpfr_init_set_ui(D3, 3, RND);
    mpfr_init_set_ui(D4, 4, RND);
}

static void septic (series f, series value, int n) {
    ad_set(f, value, n);
    mpfr_add_ui(f[0], f[0], 4, RND);
    ad_prod(___, value, f, n);
    mpfr_sub_ui(___[0], ___[0], 54, RND);
    ad_prod(f, value, ___, n);
    mpfr_sub_ui(f[0], f[0], 160, RND);
    ad_prod(___, value, f, n);
    mpfr_add_ui(___[0], ___[0], 641, RND);
    ad_prod(f, value, ___, n);
    mpfr_add_ui(f[0], f[0], 828, RND);
    ad_prod(___, value, f, n);
    mpfr_sub_ui(___[0], ___[0], 1260, RND);
    ad_prod(f, value, ___, n);
}

typedef struct { mpfr_t a; } parameters;

void *get_p (int argc, char **argv, long n) {
    (void)argc, (void)argv, (void)n;
    parameters *p = malloc(sizeof (parameters));
    return p;
}

void ode (series x, series y, series z, components *c, void *params, int k) {
    (void)x, (void)y, (void)z, (void)c, (void)params, (void)k;
}

int main (int argc, char **argv) {
    int n;
    mpfr_t x, y, _;

    assert(argc == 4);
    mpfr_set_default_prec(113);
    ad_test_tempvars();
    n = (int)strtol(argv[1], NULL, BASE) + 1;
    assert(n > 1);
    mpfr_init_set_str(x, argv[2], BASE, RND);
    mpfr_init_set_str(y, argv[3], BASE, RND);
    mpfr_init(_);

    ___ = t_jet(n);

    series cx = ad_series_v(n, x);
    series cx0 = ad_series_v(n, D0);
    series cy = ad_series_c(n, y);

    mpfr_const_pi(_, RND);
    mpfr_div_ui(_, _, 3, RND);
    series PI_3 = ad_series_v(n, _);

    mpfr_const_pi(_, RND);
    mpfr_div_ui(_, _, 4, RND);
    series PI_4 = ad_series_v(n, _);

    series wsqr = t_jet(n);
    series wabs = t_jet(n);
    series wsqrt = t_jet(n);
    series wsum = t_jet(n);
    series wprod = t_jet(n);
    series wquot = t_jet(n);
    series wpwr = t_jet(n);
    series we = t_jet(n);
    series wl = t_jet(n);
    series ws = t_jet(n);
    series wc = t_jet(n);
    series wt = t_jet(n);
    series ws2 = t_jet(n);
    series __ = ad_series_c(n, D0);

    printf("%s%s%s\n", KCYN, "Horner", KNRM);
    mpfr_set_si(__[0], -19, RND);
    mpfr_set_ui(__[1], 7, RND);
    mpfr_set_si(__[2], -4, RND);
    mpfr_set_ui(__[3], 6, RND);
    t_horner(__, 3, D3);
    mpfr_printf(" %7.3RNf\n", *t_horner(__, 3, D3));
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "x", KNRM);
    derivative_output(cx, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "y", KNRM);
    derivative_output(cy, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^2.0", KNRM);
    ad_pwr(wpwr, cx, D2, n);
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
    ad_pwr(wpwr, cx, D4, n);
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
    ad_pwr(wpwr, cx, D05, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = sqrt(x)", KNRM);
    ad_sqrt(wsqrt, cx, n);
    jet_output(wsqrt, n, KNRM, KGRY);
    derivative_output(wsqrt, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^-0.5", KNRM);
    ad_pwr(wpwr, cx, D_05, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = 1 / sqrt(x)", KNRM);
    ad_sqrt(wsqrt, cx, n);
    ad_inv(wquot, wsqrt, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^-1", KNRM);
    ad_pwr(wpwr, cx, D_1, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = 1 / x", KNRM);
    ad_inv(wquot, cx, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^0", KNRM);
    ad_pwr(wpwr, cx, D0, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = x / x", KNRM);
    ad_quot(wquot, cx, cx, n);
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

    ad_sin_cos(ws, wc, cx0, TRIG, n);
    printf("%s%s%s\n", KCYN, "f(x) = sin(0)", KNRM);
    jet_output(ws, n, KNRM, KGRY);
    derivative_output(ws, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = cos(0)", KNRM);
    jet_output(wc, n, KNRM, KGRY);
    derivative_output(wc, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    ad_sin_cos(ws, wc, PI_3, TRIG, n);
    printf("%s%s%s\n", KCYN, "f(x) = sin(pi/3)", KNRM);
    jet_output(ws, n, KNRM, KGRY);
    derivative_output(ws, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = cos(pi/3)", KNRM);
    jet_output(wc, n, KNRM, KGRY);
    derivative_output(wc, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    ad_sin_cos(ws, wc, PI_3, HYP, n);
    printf("%s%s%s\n", KCYN, "f(x) = sinh(pi/3)", KNRM);
    jet_output(ws, n, KNRM, KGRY);
    derivative_output(ws, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = cosh(pi/3)", KNRM);
    jet_output(wc, n, KNRM, KGRY);
    derivative_output(wc, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    ad_tan_sec2(wt, ws2, PI_4, TRIG, n);
    printf("%s%s%s\n", KCYN, "f(x) = tan(pi/4)", KNRM);
    jet_output(wt, n, KNRM, KGRY);
    derivative_output(wt, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = sec(pi/4)^2", KNRM);
    jet_output(ws2, n, KNRM, KGRY);
    derivative_output(ws2, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    ad_tan_sec2(wt, ws2, PI_4, HYP, n);
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
    ad_sin_cos(ws, wc, wsqr, TRIG, n);
    ad_prod(wprod, cy, ws, n);
    ad_minus(wsum, cx, wprod, n);
    ad_exp(we, cx, n);
    ad_quot(wquot, we, wsum, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = septic(x)^-1", KNRM);
    septic(__, cx, n);
    ad_pwr(wpwr, __, D_1, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    printf("%s", KNRM);
    printf("%s%s%s\n", KCYN, "f(x) = 1 / septic(x)", KNRM);
    ad_inv(wquot, __, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    cx = ad_series_c(n, x);
    cy = ad_series_v(n, y);

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
    ad_sin_cos(ws, wc, wsqr, TRIG, n);
    ad_prod(wprod, cy, ws, n);
    ad_minus(wsum, cx, wprod, n);
    ad_exp(we, cx, n);
    ad_quot(wquot, we, wsum, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    return 0;
}

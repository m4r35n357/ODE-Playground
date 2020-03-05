/*
 * Automatic Differentiation of Taylor Series, validation checks
 *
 * Example: ./libad-test-dbg 32 20 2 1e-18
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
#define KGRN "\x1B[1;32m"
#define KCYN "\x1B[1;36m"
#define KRED "\x1B[1;31m"

typedef enum {PASS, FAIL} result;

int total = 0, passed = 0, skipped = 0;

static result check_jets (char* name, mpfr_t *a, mpfr_t *b, int size, mpfr_t threshold, mpfr_t *_) {
    total++;
    for (int k = 0; k < size; k++) {
        mpfr_sub(*_, a[k], b[k], RND);
        if (mpfr_cmp_abs(*_, threshold) > 0) {
            double da = mpfr_get_d(a[k], RND);
            double db = mpfr_get_d(b[k], RND);
            double difference = mpfr_get_d(*_, RND);
            printf("%sFAILED%s %s,  k: %d diff: %+.3e a: %+.3e b: %+.3e\n", KRED, KNRM, name, k, difference, da, db);
            return FAIL;
        }
    }
    printf("%sPASSED%s %s\n", KGRN, KNRM, name);
    passed++;
    return PASS;
}

int main (int argc, char **argv) {
    mpfr_t x0, DA, D05, D0, D1, D2, tmp, tol;

    assert(argc == 5);

    mpfr_init(tmp);
    mpfr_set_default_prec(strtod(argv[1], NULL) * 3.322);
    long n = strtol(argv[2], NULL, BASE) + 1;
    assert(n > 1);
    mpfr_init_set_str(x0, argv[3], BASE, RND);
    mpfr_init_set_str(tol, argv[4], BASE, RND);

    double a = -3.7;
    mpfr_init_set_d(DA, a, RND);
    mpfr_init_set_d(D05, 0.5, RND);
    mpfr_init_set_ui(D0, 0, RND);
    mpfr_init_set_ui(D1, 1, RND);
    mpfr_init_set_ui(D2, 2, RND);

    mpfr_t *wplus = t_jet(n);
    mpfr_t *wminus = t_jet(n);
    mpfr_t *wscale = t_jet(n);
    mpfr_t *wsqr1 = t_jet(n);
    mpfr_t *wsqr2 = t_jet(n);
    mpfr_t *wabs = t_jet(n);
    mpfr_t *wsqrt = t_jet(n);
    mpfr_t *wprod1 = t_jet(n);
    mpfr_t *wprod2 = t_jet(n);
    mpfr_t *wquot = t_jet(n);
    mpfr_t *wpwr1 = t_jet(n);
    mpfr_t *wpwr2 = t_jet(n);
    mpfr_t *we1 = t_jet(n);
    mpfr_t *we2 = t_jet(n);
    mpfr_t *wl1 = t_jet(n);
    mpfr_t *wl2 = t_jet(n);
    mpfr_t *ws = t_jet(n);
    mpfr_t *wc = t_jet(n);
    mpfr_t *wt = t_jet(n);
    mpfr_t *ws2 = t_jet(n);

    printf("\n");

    mpfr_t *x = t_jet_c(n, x0);
    set_ad_status(x, VARIABLE);

    mpfr_t *c0 = t_jet_c(n, D0);
    mpfr_t *c1 = t_jet_c(n, D1);

    check_jets("x * x == sqr(x)", ad_product(wprod1, x, x, n), ad_square(wsqr1, x, n), n, tol, &tmp);
    if (mpfr_sgn(x[0]) > 0) {
        check_jets("x * x == x^2", wprod1, ad_power(wpwr1, x, 2.0, n), n, tol, &tmp);
    } else skipped++;

    check_jets("x * x * x * x == sqr(sqr(x))", ad_product(wprod2, wprod1, wprod1, n), ad_square(wsqr2, wsqr1, n), n, tol, &tmp);
    if (mpfr_sgn(x[0]) > 0) {
        check_jets("x * x * x * x == x^4", wprod2, ad_power(wpwr2, wpwr1, 2.0, n), n, tol, &tmp);
    } else skipped++;

    if (mpfr_sgn(x[0]) > 0) {
        check_jets("x^0.5 == sqrt(x)", ad_power(wpwr1, x, 0.5, n), ad_sqrt(wsqrt, x, n), n, tol, &tmp);
    } else skipped++;
    if (mpfr_sgn(x[0]) > 0) {
        check_jets("x^-0.5 == 1 / sqrt(x)", ad_power(wpwr1, x, -0.5, n), ad_quotient(wquot, c1, ad_sqrt(wsqrt, x, n), n), n, tol, &tmp);
    } else skipped++;
    if (mpfr_sgn(x[0]) > 0) {
        check_jets("x^-1 == 1 / x", ad_power(wpwr1, x, -1.0, n), ad_quotient(wquot, c1, x, n), n, tol, &tmp);
    } else skipped++;
    if (mpfr_sgn(x[0]) > 0) {
        check_jets("x^0 == x / x", ad_power(wpwr1, x, 0.0, n), ad_quotient(wquot, x, ad_scale(wscale, x, D1, n), n), n, tol, &tmp);
    } else skipped++;

    if (mpfr_zero_p(x[0]) == 0) {
        check_jets("sqrt(x * x) == |x|", ad_sqrt(wsqrt, ad_product(wprod1, x, x, n), n), ad_abs(wabs, x, n), n, tol, &tmp);
    } else skipped++;

    if (mpfr_sgn(x[0]) > 0) {
        check_jets("log(x^a) == a * log(x)", ad_ln(wl1, ad_power(wpwr1, x, a, n), n), ad_scale(wscale, ad_ln(wl2, x, n), DA, n), n, tol, &tmp);
    } else skipped++;
    check_jets("log(e^x) == x", ad_ln(wl1, ad_exp(we1, x, n), n), x, n, tol, &tmp);

    ad_sin_cos(ws, wc, x, n, TRIG);
    check_jets("cos^2(x) + sin^2(x) == 1", ad_plus(wplus, ad_square(wsqr1, wc, n), ad_square(wsqr2, ws, n), n), c1, n, tol, &tmp);

    ad_sin_cos(ws, wc, x, n, HYP);
    check_jets("cosh^2(x) - sinh^2(x) == 1", ad_minus(wminus, ad_square(wsqr1, wc, n), ad_square(wsqr2, ws, n), n), c1, n, tol, &tmp);

    ad_exp(we1, x, n);
    ad_exp(we2, ad_minus(wminus, c0, x, n), n);
    check_jets("sinh(x) == 0.5 * (e^x - e^-x)", ws, ad_scale(wscale, ad_minus(wminus, we1, we2, n), D05, n), n, tol, &tmp);
    check_jets("cosh(x) == 0.5 * (e^x + e^-x)", wc, ad_scale(wscale, ad_plus(wplus, we1, we2, n), D05, n), n, tol, &tmp);

    ad_tan_sec2(wt, ws2, x, n, TRIG);
    check_jets("sec^2(x) - tan^2(x) == 1", ad_minus(wminus, ws2, ad_square(wsqr1, wt, n), n), c1, n, tol, &tmp);

    ad_sin_cos(ws, wc, x, n, TRIG);
    check_jets("tan(x) == sin(x) / cos(x)", wt, ad_quotient(wquot, ws, wc, n), n, tol, &tmp);
    check_jets("sec^2(x) == 1 / cos^2(x)", ws2, ad_quotient(wquot, c1, ad_square(wsqr1, wc, n), n), n, tol, &tmp);

    ad_tan_sec2(wt, ws2, x, n, HYP);
    check_jets("sech^2(x) + tanh^2(x) == 1", ad_plus(wplus, ws2, ad_square(wsqr1, wt, n), n), c1, n, tol, &tmp);

    ad_sin_cos(ws, wc, x, n, HYP);
    check_jets("tanh(x) == sinh(x) / cosh(x)", wt, ad_quotient(wquot, ws, wc, n), n, tol, &tmp);
    check_jets("sech^2(x) == 1 / cosh^2(x)", ws2, ad_quotient(wquot, c1, ad_square(wsqr1, wc, n), n), n, tol, &tmp);

    printf("Total: %d, Passed: %d", total, passed);
    if (skipped > 0) {
        printf(", %sSKIPPED%s %d", KCYN, KNRM, skipped);
    }
    if (passed == total) {
        printf("\n\n");
        return 0;
    } else {
        printf(", %sFAILED%s %d\n\n", KRED, KNRM, total - passed);
        return 1;
    }
}

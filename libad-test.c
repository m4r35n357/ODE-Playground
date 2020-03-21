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
#define KYLW "\x1B[1;33m"
#define KRED "\x1B[1;31m"

typedef enum {PASS, FAIL} result;

static int total = 0, passed = 0, skipped = 0;

static mpfr_t delta;

static result compare (char* name, mpfr_t *a, mpfr_t *b, int size, mpfr_t threshold) {
    total++;
    for (int k = 0; k < size; k++) {
        mpfr_sub(delta, a[k], b[k], RND);
        if (mpfr_cmp_abs(delta, threshold) > 0) {
            printf("%sFAILED%s %s  k: %d  diff: %.3e  a: %.3e  b: %.3e\n",
                    KRED, KNRM, name, k, mpfr_get_d(delta, RND), mpfr_get_d(a[k], RND), mpfr_get_d(b[k], RND));
            return FAIL;
        }
    }
    printf("%sPASSED%s %s\n", KGRN, KNRM, name);
    passed++;
    return PASS;
}

int main (int argc, char **argv) {
    mpfr_t x0, DA, D05, D0, D1, D2, tol;

    assert(argc == 5);

    mpfr_init(delta);
    mpfr_set_default_prec(strtod(argv[1], NULL) * 3.322);
    ad_tempvars();
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

    mpfr_t *plus = t_jet(n);
    mpfr_t *minus = t_jet(n);
    mpfr_t *neg = t_jet(n);
    mpfr_t *scale = t_jet(n);
    mpfr_t *sqr1 = t_jet(n);
    mpfr_t *sqr2 = t_jet(n);
    mpfr_t *abs = t_jet(n);
    mpfr_t *sqrt = t_jet(n);
    mpfr_t *prod = t_jet(n);
    mpfr_t *quot = t_jet(n);
    mpfr_t *pwr = t_jet(n);
    mpfr_t *e1 = t_jet(n);
    mpfr_t *e2 = t_jet(n);
    mpfr_t *l1 = t_jet(n);
    mpfr_t *l2 = t_jet(n);
    mpfr_t *s = t_jet(n);
    mpfr_t *c = t_jet(n);
    mpfr_t *t = t_jet(n);
    mpfr_t *s2 = t_jet(n);

    printf("\n");

    mpfr_t *c1 = t_jet_c(n, D1);
    mpfr_t *x = t_jet_c(n, x0);
    set_ad_status(x, VARIABLE);

    compare("x * x == sqr(x)", ad_prod(prod, x, x, n), ad_sqr(sqr1, x, n), n, tol);
    if (mpfr_sgn(x[0]) > 0) {
        compare("x * x == x^2", prod, ad_pwr(pwr, x, 2.0, n), n, tol);
    } else skipped++;

    if (mpfr_sgn(x[0]) > 0) {
        compare("x^1 == x", ad_pwr(pwr, x, 1.0, n), x, n, tol);
    } else skipped++;
    if (mpfr_sgn(x[0]) > 0) {
        compare("x^0.5 == sqrt(x)", ad_pwr(pwr, x, 0.5, n), ad_sqrt(sqrt, x, n), n, tol);
    } else skipped++;
    if (mpfr_sgn(x[0]) > 0) {
        compare("x^0 == x / x", ad_pwr(pwr, x, 0.0, n), ad_quot(quot, x, ad_scale(scale, x, D1, n), n), n, tol);
    } else skipped++;
    if (mpfr_sgn(x[0]) > 0) {
        compare("x^-0.5 == 1 / sqrt(x)", ad_pwr(pwr, x, -0.5, n), ad_quot(quot, c1, ad_sqrt(sqrt, x, n), n), n, tol);
    } else skipped++;
    if (mpfr_sgn(x[0]) > 0) {
        compare("x^-1 == 1 / x", ad_pwr(pwr, x, -1.0, n), ad_quot(quot, c1, x, n), n, tol);
    } else skipped++;

    if (mpfr_zero_p(x[0]) == 0) {
        compare("sqrt(x * x) == |x|", ad_sqrt(sqrt, ad_prod(prod, x, x, n), n), ad_abs(abs, x, n), n, tol);
    } else skipped++;

    if (mpfr_sgn(x[0]) > 0) {
        compare("log(x^a) == a * log(x)", ad_ln(l1, ad_pwr(pwr, x, a, n), n), ad_scale(scale, ad_ln(l2, x, n), DA, n), n, tol);
    } else skipped++;
    compare("log(e^x) == x", ad_ln(l1, ad_exp(e1, x, n), n), x, n, tol);

    ad_sin_cos(s, c, x, n, TRIG);
    compare("cos^2(x) + sin^2(x) == 1", ad_plus(plus, ad_sqr(sqr1, c, n), ad_sqr(sqr2, s, n), n), c1, n, tol);

    ad_sin_cos(s, c, x, n, HYP);
    compare("cosh^2(x) - sinh^2(x) == 1", ad_minus(minus, ad_sqr(sqr1, c, n), ad_sqr(sqr2, s, n), n), c1, n, tol);

    ad_exp(e1, x, n);
    ad_exp(e2, ad_neg(neg, x, n), n);
    compare("sinh(x) == 0.5 * (e^x - e^-x)", s, ad_scale(scale, ad_minus(minus, e1, e2, n), D05, n), n, tol);
    compare("cosh(x) == 0.5 * (e^x + e^-x)", c, ad_scale(scale, ad_plus(plus, e1, e2, n), D05, n), n, tol);

    ad_tan_sec2(t, s2, x, n, TRIG);
    compare("sec^2(x) - tan^2(x) == 1", ad_minus(minus, s2, ad_sqr(sqr1, t, n), n), c1, n, tol);

    ad_sin_cos(s, c, x, n, TRIG);
    compare("tan(x) == sin(x) / cos(x)", t, ad_quot(quot, s, c, n), n, tol);
    compare("sec^2(x) == 1 / cos^2(x)", s2, ad_quot(quot, c1, ad_sqr(sqr1, c, n), n), n, tol);

    ad_tan_sec2(t, s2, x, n, HYP);
    compare("sech^2(x) + tanh^2(x) == 1", ad_plus(plus, s2, ad_sqr(sqr1, t, n), n), c1, n, tol);

    ad_sin_cos(s, c, x, n, HYP);
    compare("tanh(x) == sinh(x) / cosh(x)", t, ad_quot(quot, s, c, n), n, tol);
    compare("sech^2(x) == 1 / cosh^2(x)", s2, ad_quot(quot, c1, ad_sqr(sqr1, c, n), n), n, tol);

    printf("Total: %d, Passed: %d", total, passed);
    if (skipped > 0) {
        printf(", %sSKIPPED%s %d", KYLW, KNRM, skipped);
    }
    if (passed == total) {
        printf("\n\n");
        return 0;
    } else {
        printf(", %sFAILED%s %d\n\n", KRED, KNRM, total - passed);
        return 1;
    }
}

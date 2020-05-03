/*
 * Automatic Differentiation of Taylor Seriesewest validation checks
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

static mpfr_t delta, D0, D05, D_05, D1, D_1, D2, D_2, D3, D_3, D_5, DA;

static void ad_lib_test_tempvars (void) {
    ad_tempvars();
    mpfr_init(delta);
    mpfr_init_set_ui(D0, 0, RND);
    mpfr_init_set_str(D05, "0.5", BASE, RND);
    mpfr_init_set_str(D_05, "-0.5", BASE, RND);
    mpfr_init_set_ui(D1, 1, RND);
    mpfr_init_set_si(D_1, -1, RND);
    mpfr_init_set_ui(D2, 2, RND);
    mpfr_init_set_si(D_2, -2, RND);
    mpfr_init_set_ui(D3, 3, RND);
    mpfr_init_set_si(D_3, -3, RND);
    mpfr_init_set_si(D_5, -5, RND);
    mpfr_init_set_str(DA, "-3.7", BASE, RND);
}

static result compare (char* name, series a, series b, mpfr_t threshold) {
    assert(b.size == a.size);
    total++;
    for (int k = 0; k < b.size; k++) {
        mpfr_sub(delta, a.a[k], b.a[k], RND);
        if (mpfr_cmp_abs(delta, threshold) > 0) {
            printf("%sFAILED%s %s  k: %d  diff: %.3e  a: %.3e  b: %.3e\n",
                    KRED, KNRM, name, k, mpfr_get_d(delta, RND), mpfr_get_d(a.a[k], RND), mpfr_get_d(b.a[k], RND));
            return FAIL;
        }
    }
    printf("%sPASSED%s %s\n", KGRN, KNRM, name);
    passed++;
    return PASS;
}

int main (int argc, char **argv) {
    mpfr_t x0, tol;

    assert(argc == 5);

    mpfr_set_default_prec(strtod(argv[1], NULL) * 3.322);
    ad_lib_test_tempvars();
    long n = strtol(argv[2], NULL, BASE) + 1;
    assert(n > 1);
    mpfr_init_set_str(x0, argv[3], BASE, RND);
    mpfr_init_set_str(tol, argv[4], BASE, RND);

    series a = t_jet(n);
    series p = t_jet(n);
    series q = t_jet(n);
    series m = t_jet(n);
    series s1 = t_jet(n);
    series s2 = t_jet(n);
    series p1 = t_jet(n);
    series p2 = t_jet(n);
    series p3 = t_jet(n);
    series e1 = t_jet(n);
    series e2 = t_jet(n);
    series l1 = t_jet(n);
    series l2 = t_jet(n);
    series s = t_jet(n);
    series c = t_jet(n);
    series t = t_jet(n);

    series c1 = t_jet_c(n, D1);
    series x = t_jet_c(n, x0);
    set_ad_status(x, VARIABLE);

    printf("\n");
    compare("x * x == sqr(x)", ad_prod(p, x, x), ad_sqr(s1, x), tol);
    if (mpfr_sgn(x.a[0]) > 0) {
        compare("x * x == x^2", p, ad_pwr(p1, x, D2), tol);
    } else skipped++;

    if (mpfr_sgn(x.a[0]) > 0) {
        compare("x^1 == x", ad_pwr(p1, x, D1), x, tol);
    } else skipped++;
    if (mpfr_sgn(x.a[0]) > 0) {
        compare("x^0.5 == sqrt(x)", ad_pwr(p1, x, D05), ad_sqrt(s, x), tol);
    } else skipped++;
    if (mpfr_sgn(x.a[0]) > 0) {
        compare("x^0 == x / x", ad_pwr(p1, x, D0), ad_quot(q, x, x), tol);
    } else skipped++;
    if (mpfr_sgn(x.a[0]) > 0) {
        compare("x^-0.5 == 1 / sqrt(x)", ad_pwr(p1, x, D_05), ad_quot(q, c1, ad_sqrt(s, x)), tol);
    } else skipped++;
    if (mpfr_sgn(x.a[0]) > 0) {
        compare("x^-1 == 1 / x", ad_pwr(p1, x, D_1), ad_quot(q, c1, x), tol);
    } else skipped++;
    if (mpfr_sgn(x.a[0]) > 0) {
        compare("x^2 * x^-5 == x^-3", ad_prod(p, ad_pwr(p1, x, D2), ad_pwr(p2, x, D_5)), ad_pwr(p3, x, D_3), tol);
    } else skipped++;

    if (mpfr_zero_p(x.a[0]) == 0) {
        compare("sqrt(x * x) == |x|", ad_sqrt(s, ad_prod(p, x, x)), ad_abs(a, x), tol);
    } else skipped++;

    if (mpfr_sgn(x.a[0]) > 0) {
        compare("log(x^a) == a * log(x)", ad_ln(l1, ad_pwr(p1, x, DA)), ad_scale(s, ad_ln(l2, x), DA), tol);
    } else skipped++;
    compare("log(e^x) == x", ad_ln(l1, ad_exp(e1, x)), x, tol);

    ad_sin_cos(s, c, x, TRIG);
    compare("cos^2(x) + sin^2(x) == 1", ad_plus(p, ad_sqr(s1, c), ad_sqr(s2, s)), c1, tol);

    ad_sin_cos(s, c, x, HYP);
    compare("cosh^2(x) - sinh^2(x) == 1", ad_minus(m, ad_sqr(s1, c), ad_sqr(s2, s)), c1, tol);

    ad_exp(e1, x);
    ad_exp(e2, ad_neg(m, x));
    compare("sinh(x) == 0.5 * (e^x - e^-x)", s, ad_scale(s, ad_minus(m, e1, e2), D05), tol);
    compare("cosh(x) == 0.5 * (e^x + e^-x)", c, ad_scale(s, ad_plus(p, e1, e2), D05), tol);

    ad_tan_sec2(t, s2, x, TRIG);
    compare("sec^2(x) - tan^2(x) == 1", ad_minus(m, s2, ad_sqr(s1, t)), c1, tol);

    ad_sin_cos(s, c, x, TRIG);
    compare("tan(x) == sin(x) / cos(x)", t, ad_quot(q, s, c), tol);
    compare("sec^2(x) == 1 / cos^2(x)", s2, ad_quot(q, c1, ad_sqr(s1, c)), tol);

    ad_tan_sec2(t, s2, x, HYP);
    compare("sech^2(x) + tanh^2(x) == 1", ad_plus(p, s2, ad_sqr(s1, t)), c1, tol);

    ad_sin_cos(s, c, x, HYP);
    compare("tanh(x) == sinh(x) / cosh(x)", t, ad_quot(q, s, c), tol);
    compare("sech^2(x) == 1 / cosh^2(x)", s2, ad_quot(q, c1, ad_sqr(s1, c)), tol);

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

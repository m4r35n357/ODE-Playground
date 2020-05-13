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

static mpfr_t delta, tolerance, D0, D05, D_05, D1, D_1, D2, D_2, D3, D_3, D_5, DA;

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

static result compare (char* name, series a, series b) {
    assert(b.size == a.size);
    total++;
    for (int k = 0; k < b.size; k++) {
        mpfr_sub(delta, a.a[k], b.a[k], RND);
        if (mpfr_cmp_abs(delta, tolerance) > 0) {
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
    mpfr_t x0;

    assert(argc == 5);
    mpfr_set_default_prec(strtod(argv[1], NULL) * 3.322);
    ad_lib_test_tempvars();
    long n = strtol(argv[2], NULL, BASE) + 1;
    assert(n > 1);
    mpfr_init_set_str(x0, argv[3], BASE, RND);
    mpfr_init_set_str(tolerance, argv[4], BASE, RND);

    series abs = t_jet(n);
    series scale = t_jet(n);
    series sum = t_jet(n);
    series diff = t_jet(n);
    series prod = t_jet(n);
    series quot = t_jet(n);
    series neg = t_jet(n);
    series sqr1 = t_jet(n);
    series sqr2 = t_jet(n);
    series sqrt = t_jet(n);
    series pow1 = t_jet(n);
    series pow2 = t_jet(n);
    series pow3 = t_jet(n);
    series exp1 = t_jet(n);
    series exp2 = t_jet(n);
    series ln1 = t_jet(n);
    series ln2 = t_jet(n);
    series sin = t_jet(n);
    series cos = t_jet(n);
    series tan = t_jet(n);
    series sec2 = t_jet(n);

    series c1 = t_jet_c(n, D1);
    series x = t_jet_c(n, x0);
    set_ad_status(x, VARIABLE);

    printf("\n");
    compare("x * x == sqr(x)", ad_prod(prod, x, x), ad_sqr(sqr1, x));
    if (mpfr_zero_p(x.a[0]) == 0) {
        compare("sqr(x) / x == x", ad_quot(quot, sqr1, x), x);
    } else skipped++;
    if (mpfr_sgn(x.a[0]) > 0) {
        ad_sqrt(sqrt, x);
        compare("sqrt(x) * sqrt(x) == x", ad_prod(prod, sqrt, sqrt), x);
    } else skipped++;
    if (mpfr_sgn(x.a[0]) > 0) {
        compare("x / sqrt(x) == sqrt(x)", ad_quot(quot, x, sqrt), sqrt);
    } else skipped++;

    if (mpfr_sgn(x.a[0]) > 0) {
        compare("x^2 == sqr(x)", ad_pwr(pow1, x, D2), sqr1);
    } else skipped++;
    if (mpfr_sgn(x.a[0]) > 0) {
        compare("x^1 == x", ad_pwr(pow1, x, D1), x);
    } else skipped++;
    if (mpfr_sgn(x.a[0]) > 0) {
        compare("x^0.5 == sqrt(x)", ad_pwr(pow1, x, D05), sqrt);
    } else skipped++;
    if (mpfr_sgn(x.a[0]) > 0) {
        compare("x^0 == x / x", ad_pwr(pow1, x, D0), ad_quot(quot, x, x));
    } else skipped++;
    if (mpfr_sgn(x.a[0]) > 0) {
        compare("x^-0.5 == 1 / sqrt(x)", ad_pwr(pow1, x, D_05), ad_quot(quot, c1, sqrt));
    } else skipped++;
    if (mpfr_sgn(x.a[0]) > 0) {
        compare("x^-1 == 1 / x", ad_pwr(pow1, x, D_1), ad_quot(quot, c1, x));
    } else skipped++;
    if (mpfr_sgn(x.a[0]) > 0) {
        compare("x^-2 == 1 / sqr(x)", ad_pwr(pow1, x, D_2), ad_quot(quot, c1, sqr1));
    } else skipped++;

    if (mpfr_sgn(x.a[0]) > 0) {
        compare("x^2 * x^-5 == x^-3", ad_prod(prod, ad_pwr(pow1, x, D2), ad_pwr(pow2, x, D_5)), ad_pwr(pow3, x, D_3));
    } else skipped++;

    if (mpfr_zero_p(x.a[0]) == 0) {
        compare("(x * x)^0.5 == |x|", ad_pwr(pow1, ad_prod(prod, x, x), D05), ad_abs(abs, x));
    } else skipped++;

    compare("log(e^x) == x", ad_ln(ln1, ad_exp(exp1, x)), x);
    if (mpfr_sgn(x.a[0]) > 0) {
        compare("log(x^a) == a * log(x)", ad_ln(ln1, ad_pwr(pow1, x, DA)), ad_scale(scale, ad_ln(ln2, x), DA));
    } else skipped++;

    ad_sin_cos(sin, cos, x, TRIG);
    compare("cos^2(x) + sin^2(x) == 1", ad_plus(sum, ad_sqr(sqr1, cos), ad_sqr(sqr2, sin)), c1);

    ad_sin_cos(sin, cos, x, HYP);
    compare("cosh^2(x) - sinh^2(x) == 1", ad_minus(diff, ad_sqr(sqr1, cos), ad_sqr(sqr2, sin)), c1);

    ad_exp(exp1, x);
    ad_exp(exp2, ad_neg(neg, x));
    compare("sinh(x) == 0.5 * (e^x - e^-x)", sin, ad_scale(scale, ad_minus(diff, exp1, exp2), D05));
    compare("cosh(x) == 0.5 * (e^x + e^-x)", cos, ad_scale(scale, ad_plus(sum, exp1, exp2), D05));

    ad_tan_sec2(tan, sec2, x, TRIG);
    compare("sec^2(x) - tan^2(x) == 1", ad_minus(diff, sec2, ad_sqr(sqr1, tan)), c1);

    ad_sin_cos(sin, cos, x, TRIG);
    compare("tan(x) == sin(x) / cos(x)", tan, ad_quot(quot, sin, cos));
    compare("sec^2(x) == 1 / cos^2(x)", sec2, ad_quot(quot, c1, ad_sqr(sqr1, cos)));

    ad_tan_sec2(tan, sec2, x, HYP);
    compare("sech^2(x) + tanh^2(x) == 1", ad_plus(sum, sec2, ad_sqr(sqr1, tan)), c1);

    ad_sin_cos(sin, cos, x, HYP);
    compare("tanh(x) == sinh(x) / cosh(x)", tan, ad_quot(quot, sin, cos));
    compare("sech^2(x) == 1 / cosh^2(x)", sec2, ad_quot(quot, c1, ad_sqr(sqr1, cos)));

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

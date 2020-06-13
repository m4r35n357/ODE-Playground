/*
 * Automatic Differentiation of Taylor Seriesewest validation checks
 *
 * Example: ./libad-test-dbg 9 32 20 2 1e-18
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

typedef enum {PASSED, SKIPPED, FAILED} result;

static int debug = 0, total = 0, passed = 0, skipped = 0;

static mpfr_t x0, delta, tolerance, D0, D05, D_05, D1, D_1, D2, D_2, D_3;

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
    mpfr_init_set_si(D_3, -3, RND);
}

static result skip (char* name) {
    total++;
    skipped++;
    printf("%sSKIPPED%s %s\n", KYLW, KNRM, name);
    return SKIPPED;
}

static result compare (char* name, series a, series b) {
    assert(b.size == a.size);
    total++;
    for (int k = 0; k < b.size; k++) {
        mpfr_sub(delta, a.jet[k], b.jet[k], RND);
        if (mpfr_cmp_abs(delta, tolerance) > 0) {
            printf("%s FAILED%s %s  k: %d  LHS: %.6e  RHS: %.6e  diff %.3e\n",
                    KRED, KNRM, name, k, mpfr_get_d(a.jet[k], RND), mpfr_get_d(b.jet[k], RND), mpfr_get_d(delta, RND));
            return FAILED;
        }
        if (debug >= 2 && k == 0) printf("\n");
        if (debug >= 2) printf("%s  DEBUG%s  k: %2d  %+.6e %+.6e  diff %+.3e\n",
                               KNRM, KNRM, k, mpfr_get_d(a.jet[k], RND), mpfr_get_d(b.jet[k], RND), mpfr_get_d(delta, RND));
    }
    if (debug >= 1) printf("%s PASSED%s %s\n", KGRN, KNRM, name);
    passed++;
    return PASSED;
}

int main (int argc, char **argv) {
    mpfr_t PI, PI_2;

    assert(argc == 5 || argc == 6);
    mpfr_set_default_prec(strtod(argv[1], NULL) * 3.322);
    ad_lib_test_tempvars();
    long n = strtol(argv[2], NULL, BASE) + 1;
    assert(n > 1);
    mpfr_init_set_str(x0, argv[3], BASE, RND);
    mpfr_init_set_str(tolerance, argv[4], BASE, RND);
    if (argc == 6) debug = strtol(argv[5], NULL, BASE);

    series abs = t_series(n);
    series scale = t_series(n);
    series sum = t_series(n);
    series diff = t_series(n);
    series prod = t_series(n);
    series quot = t_series(n);
    series inv = t_series(n);
    series neg = t_series(n);
    series sqr1 = t_series(n);
    series sqr2 = t_series(n);
    series sqrt = t_series(n);
    series pow1 = t_series(n);
    series pow2 = t_series(n);
    series exp1 = t_series(n);
    series exp2 = t_series(n);
    series ln1 = t_series(n);
    series ln2 = t_series(n);
    series sin = t_series(n);
    series sin2 = t_series(n);
    series cos = t_series(n);
    series cos2 = t_series(n);
    series tan = t_series(n);
    series sec2 = t_series(n);

    series c1 = ad_series_c(n, D1);
    series x = ad_series_v(n, x0);

    mpfr_inits(PI, PI_2, NULL);
    mpfr_const_pi(PI, RND);
    mpfr_div_2ui(PI_2, PI, 1, RND);

    int positive = mpfr_sgn(x.jet[0]) > 0;
    int non_zero = mpfr_zero_p(x.jet[0]) == 0;
    int lt_pi = mpfr_cmpabs(x.jet[0], PI) < 0;
    int lt_pi_2 = mpfr_cmpabs(x.jet[0], PI_2) < 0;

    printf("\n");
    ad_sqr(sqr1, x);
    char* name = "x * x == sqr(x)";
    compare(name, ad_prod(prod, x, x), sqr1);
    name = "sqr(x) / x == x";
    if (non_zero) {
        compare(name, ad_quot(quot, sqr1, x), x);
    } else skip(name);
    name = "x * 1 / x == 1";
    if (non_zero) {
        compare(name, ad_prod(prod, x, ad_inv(inv, x)), c1);
    } else skip(name);
    name = "sqrt(x) * sqrt(x) == x";
    if (positive) {
        ad_sqrt(sqrt, x);
        compare(name, ad_prod(prod, sqrt, sqrt), x);
    } else skip(name);
    name = "x / sqrt(x) == sqrt(x)";
    if (positive) {
        compare(name, ad_quot(quot, x, sqrt), sqrt);
    } else skip(name);

    name = "x^2 == sqr(x)";
    if (positive) {
        compare(name, ad_pwr(pow1, x, D2), sqr1);
    } else skip(name);
    name = "x^1 == x";
    if (positive) {
        compare(name, ad_pwr(pow1, x, D1), x);
    } else skip(name);
    name = "x^0.5 == sqrt(x)";
    if (positive) {
        compare(name, ad_pwr(pow1, x, D05), sqrt);
    } else skip(name);
    name = "x^0 == 1";
    if (positive) {
        compare(name, ad_pwr(pow1, x, D0), c1);
    } else skip(name);
    name = "x^-0.5 == 1 / sqrt(x)";
    if (positive) {
        compare(name, ad_pwr(pow1, x, D_05), ad_inv(inv, sqrt));
    } else skip(name);
    name = "x^-1 == 1 / x";
    if (positive) {
        compare(name, ad_pwr(pow1, x, D_1), ad_inv(inv, x));
    } else skip(name);
    name = "x^-2 == 1 / sqr(x)";
    if (positive) {
        compare(name, ad_pwr(pow1, x, D_2), ad_inv(inv, sqr1));
    } else skip(name);

    name = "sqr(x) * x^-3 == 1 / x";
    if (positive) {
        compare(name, ad_prod(prod, sqr1, ad_pwr(pow2, x, D_3)), ad_inv(inv, x));
    } else skip(name);

    name = "sqr(x)^0.5 == |x|";
    if (non_zero) {
        compare(name, ad_pwr(pow1, sqr1, D05), ad_abs(abs, x));
    } else skip(name);

    name = "log(e^x) == x";
    compare(name, ad_ln(ln1, ad_exp(exp1, x)), x);
    name = "log(sqr(x)) == 2 * log(x)";
    if (positive) {
        ad_ln(ln2, x);
        compare(name, ad_ln(ln1, sqr1), ad_scale(scale, ln2, D2));
    } else skip(name);
    name = "log(sqrt(x)) == 0.5 * log(x)";
    if (positive) {
        compare(name, ad_ln(ln1, sqrt), ad_scale(scale, ln2, D05));
    } else skip(name);
    name = "log(1 / x) == - log(x)";
    if (positive) {
        compare(name, ad_ln(ln1, ad_inv(inv, x)), ad_neg(neg, ln2));
    } else skip(name);
    name = "log(x^-3) == - 3 * log(x)";
    if (positive) {
        compare(name, ad_ln(ln1, ad_pwr(pow1, x, D_3)), ad_scale(scale, ln2, D_3));
    } else skip(name);

    ad_sin_cos(sin, cos, x, HYP);
    ad_tan_sec2(tan, sec2, x, HYP);
    ad_sqr(sqr1, cos);
    ad_sqr(sqr2, sin);
    compare("cosh^2(x) - sinh^2(x) == 1", ad_minus(diff, sqr1, sqr2), c1);
    compare("sech^2(x) + tanh^2(x) == 1", ad_plus(sum, sec2, ad_sqr(sqr1, tan)), c1);
    compare("tanh(x) == sinh(x) / cosh(x)", tan, ad_quot(quot, sin, cos));
    compare("sech^2(x) == 1 / cosh^2(x)", sec2, ad_inv(inv, ad_sqr(sqr1, cos)));
    ad_sin_cos(sin2, cos2, ad_scale(scale, x, D2), HYP);
    compare("sinh(2 * x) == 2 * sinh(x) * cosh(x)", sin2, ad_scale(scale, ad_prod(prod, sin, cos), D2));
    compare("cosh(2 * x) == cosh^2(x) + sinh^2(x)", cos2, ad_plus(sum, sqr1, sqr2));

    ad_exp(exp1, x);
    ad_exp(exp2, ad_neg(neg, x));
    compare("sinh(x) == 0.5 * (e^x - e^-x)", sin, ad_scale(scale, ad_minus(diff, exp1, exp2), D05));
    compare("cosh(x) == 0.5 * (e^x + e^-x)", cos, ad_scale(scale, ad_plus(sum, exp1, exp2), D05));

    ad_sin_cos(sin, cos, x, TRIG);
    ad_tan_sec2(tan, sec2, x, TRIG);
    ad_sqr(sqr1, cos);
    ad_sqr(sqr2, sin);
    name = "cos^2(x) + sin^2(x) == 1";
    if (lt_pi) {
        compare(name, ad_plus(sum, sqr1, sqr2), c1);
    } else skip(name);
    name = "sec^2(x) - tan^2(x) == 1";
    if (lt_pi) {
        compare(name, ad_minus(diff, sec2, ad_sqr(sqr1, tan)), c1);
    } else skip(name);
    name = "tan(x) == sin(x) / cos(x)";
    if (lt_pi_2) {
        compare(name, tan, ad_quot(quot, sin, cos));
    } else skip(name);
    name = "sec^2(x) == 1 / cos^2(x)";
    if (lt_pi_2) {
        compare(name, sec2, ad_inv(inv, ad_sqr(sqr1, cos)));
    } else skip(name);
    ad_sin_cos(sin2, cos2, ad_scale(scale, x, D2), TRIG);
    name = "sin(2 * x) == 2 * sin(x) * cos(x)";
    if (lt_pi_2) {
        compare(name, sin2, ad_scale(scale, ad_prod(prod, sin, cos), D2));
    } else skip(name);
    name = "cos(2 * x) == cos^2(x) + sin^2(x)";
    if (lt_pi_2) {
        compare(name, cos2, ad_minus(diff, sqr1, sqr2));
    } else skip(name);

    printf("Total: %d, %sPASSED%s %d", total, KGRN, KNRM, passed);
    if (skipped > 0) {
        printf(", %sSKIPPED%s %d", KYLW, KNRM, skipped);
    }
    if (passed == total) {
        printf("\n\n");
        return 0;
    } else {
        printf(", %sFAILED%s %d\n\n", KRED, KNRM, total - passed - skipped);
        return 1;
    }
}

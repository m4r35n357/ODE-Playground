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
#define KWHT "\x1B[1;37m"
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
    mpfr_t PI_2;

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
    series inv_x = t_series(n);
    series neg = t_series(n);
    series sqr_x = t_series(n);
    series sqr_sin_x = t_series(n);
    series sqr_cos_x = t_series(n);
    series sqr_tan_x = t_series(n);
    series sqrt_x = t_series(n);
    series pow = t_series(n);
    series exp_x = t_series(n);
    series neg_exp_x = t_series(n);
    series ln = t_series(n);
    series ln_x = t_series(n);
    series sin = t_series(n);
    series sin_2x = t_series(n);
    series cos = t_series(n);
    series cos_2x = t_series(n);
    series tan = t_series(n);
    series sec2 = t_series(n);

    series c1 = ad_series_c(n, D1);
    series x = ad_series_v(n, x0);

    mpfr_init(PI_2);
    mpfr_const_pi(PI_2, RND);
    mpfr_div_2ui(PI_2, PI_2, 1, RND);

    int x_positive = mpfr_sgn(x.jet[0]) > 0;
    int x_non_zero = mpfr_zero_p(x.jet[0]) == 0;
    int x_lt_pi_2 = mpfr_cmpabs(x.jet[0], PI_2) < 0;

    printf("\n");
    ad_sqr(sqr_x, x);
    if (x_non_zero) ad_inv(inv_x, x);
    if (x_positive) ad_sqrt(sqrt_x, x);
    char* name = "x * x == sqr(x)";
    compare(name, ad_prod(prod, x, x), sqr_x);
    name = "sqr(x) / x == x";
    x_non_zero ? compare(name, ad_quot(quot, sqr_x, x), x) : skip(name);
    name = "x * 1 / x == 1";
    x_non_zero ? compare(name, ad_prod(prod, x, inv_x), c1) : skip(name);
    name = "sqrt(x) * sqrt(x) == x";
    x_positive ? compare(name, ad_prod(prod, sqrt_x, sqrt_x), x) : skip(name);
    name = "x / sqrt(x) == sqrt(x)";
    x_positive ? compare(name, ad_quot(quot, x, sqrt_x), sqrt_x) : skip(name);

    name = "x^2 == sqr(x)";
    x_positive ? compare(name, ad_pwr(pow, x, D2), sqr_x) : skip(name);
    name = "x^1 == x";
    x_positive ? compare(name, ad_pwr(pow, x, D1), x) : skip(name);
    name = "x^0.5 == sqrt(x)";
    x_positive ? compare(name, ad_pwr(pow, x, D05), sqrt_x): skip(name);
    name = "x^0 == 1";
    x_positive ? compare(name, ad_pwr(pow, x, D0), c1) : skip(name);
    name = "x^-0.5 == 1 / sqrt(x)";
    x_positive ? compare(name, ad_pwr(pow, x, D_05), ad_inv(inv, sqrt_x)) : skip(name);
    name = "x^-1 == 1 / x";
    x_positive ? compare(name, ad_pwr(pow, x, D_1), inv_x) : skip(name);
    name = "x^-2 == 1 / sqr(x)";
    x_positive ? compare(name, ad_pwr(pow, x, D_2), ad_inv(inv, sqr_x)) : skip(name);

    name = "sqr(x) * x^-3 == 1 / x";
    x_positive ? compare(name, ad_prod(prod, sqr_x, ad_pwr(pow, x, D_3)), inv_x) : skip(name);

    name = "sqr(x)^0.5 == |x|";
    x_non_zero ? compare(name, ad_pwr(pow, sqr_x, D05), ad_abs(abs, x)) : skip(name);

    name = "log(e^x) == x";
    compare(name, ad_ln(ln, ad_exp(exp_x, x)), x);
    name = "log(sqr(x)) == 2 * log(x)";
    if (x_positive) ad_ln(ln_x, x);
    x_positive ? compare(name, ad_ln(ln, sqr_x), ad_scale(scale, ln_x, D2)) : skip(name);
    name = "log(sqrt(x)) == 0.5 * log(x)";
    x_positive ? compare(name, ad_ln(ln, sqrt_x), ad_scale(scale, ln_x, D05)) : skip(name);
    name = "log(1 / x) == - log(x)";
    x_positive ? compare(name, ad_ln(ln, inv_x), ad_neg(neg, ln_x)) : skip(name);
    name = "log(x^-3) == - 3 * log(x)";
    x_positive ? compare(name, ad_ln(ln, ad_pwr(pow, x, D_3)), ad_scale(scale, ln_x, D_3)) : skip(name);

    ad_sin_cos(sin, cos, x, HYP);
    ad_tan_sec2(tan, sec2, x, HYP);
    ad_sqr(sqr_sin_x, sin);
    ad_sqr(sqr_cos_x, cos);
    ad_sin_cos(sin_2x, cos_2x, ad_scale(scale, x, D2), HYP);
    name = "cosh^2(x) - sinh^2(x) == 1";
    compare(name, ad_minus(diff, sqr_cos_x, sqr_sin_x), c1);
    name = "sech^2(x) + tanh^2(x) == 1";
    compare(name, ad_plus(sum, sec2, ad_sqr(sqr_tan_x, tan)), c1);
    name = "tanh(x) == sinh(x) / cosh(x)";
    compare(name, tan, ad_quot(quot, sin, cos));
    name = "sech^2(x) == 1 / cosh^2(x)";
    compare(name, sec2, ad_inv(inv, sqr_cos_x));
    name = "sinh(2x) == 2 * sinh(x) * cosh(x)";
    compare(name, sin_2x, ad_scale(scale, ad_prod(prod, sin, cos), D2));
    name = "cosh(2x) == cosh^2(x) + sinh^2(x)";
    compare(name, cos_2x, ad_plus(sum, sqr_cos_x, sqr_sin_x));

    ad_exp(exp_x, x);
    ad_exp(neg_exp_x, ad_neg(neg, x));
    name = "sinh(x) == 0.5 * (e^x - e^-x)";
    compare(name, sin, ad_scale(scale, ad_minus(diff, exp_x, neg_exp_x), D05));
    name = "cosh(x) == 0.5 * (e^x + e^-x)";
    compare(name, cos, ad_scale(scale, ad_plus(sum, exp_x, neg_exp_x), D05));

    ad_sin_cos(sin, cos, x, TRIG);
    ad_tan_sec2(tan, sec2, x, TRIG);
    ad_sqr(sqr_sin_x, sin);
    ad_sqr(sqr_cos_x, cos);
    ad_sin_cos(sin_2x, cos_2x, ad_scale(scale, x, D2), TRIG);
    name = "cos^2(x) + sin^2(x) == 1";
    compare(name, ad_plus(sum, sqr_cos_x, sqr_sin_x), c1);
    name = "sec^2(x) - tan^2(x) == 1";
    x_lt_pi_2 ? compare(name, ad_minus(diff, sec2, ad_sqr(sqr_tan_x, tan)), c1) : skip(name);
    name = "tan(x) == sin(x) / cos(x)";
    x_lt_pi_2 ? compare(name, tan, ad_quot(quot, sin, cos)) : skip(name);
    name = "sec^2(x) == 1 / cos^2(x)";
    x_lt_pi_2 ? compare(name, sec2, ad_inv(inv, sqr_cos_x)) : skip(name);
    name = "sin(2x) == 2 * sin(x) * cos(x)";
    compare(name, sin_2x, ad_scale(scale, ad_prod(prod, sin, cos), D2));
    name = "cos(2x) == cos^2(x) - sin^2(x)";
    compare(name, cos_2x, ad_minus(diff, sqr_cos_x, sqr_sin_x));

    printf("%sTotal%s: %d, %sPASSED%s %d", KWHT, KNRM, total, KGRN, KNRM, passed);
    if (skipped > 0) {
        printf(", %sSKIPPED%s %d", KYLW, KNRM, skipped);
    }
    if (passed == total - skipped) {
        printf("\n\n");
        return 0;
    } else {
        printf(", %sFAILED%s %d\n\n", KRED, KNRM, total - passed - skipped);
        return 1;
    }
}

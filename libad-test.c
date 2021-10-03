/*
 * Automatic Differentiation of Taylor Series, newest validation checks
 *
 * Example: ./libad-test-dbg 32 20 1 1e-18
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

static mpfr_t x0, delta, tolerance, D0, D01, D05, D_05, D1, D_1, D2, D_2, D_3;

static void ad_lib_test_tempvars (void) {
    ad_tempvars();
    mpfr_init(delta);
    mpfr_init_set_ui(D0, 0, RND);
    mpfr_init_set_str(D01, "0.1", BASE, RND);
    mpfr_init_set_str(D05, "0.5", BASE, RND);
    mpfr_init_set_str(D_05, "-0.5", BASE, RND);
    mpfr_init_set_ui(D1, 1, RND);
    mpfr_init_set_si(D_1, -1, RND);
    mpfr_init_set_ui(D2, 2, RND);
    mpfr_init_set_si(D_2, -2, RND);
    mpfr_init_set_si(D_3, -3, RND);
}

typedef struct { mpfr_t a, b, c; } parameters;

void *get_p (int argc, char **argv, long order) {
    (void)argc; (void)argv; (void)order;
    parameters *p = malloc(sizeof (parameters));
    mpfr_init_set(p->a, D1, RND);
    mpfr_init_set(p->b, D0, RND);
    mpfr_init_set(p->c, D_1, RND);
    return p;
}

void ode (series x, series y, series z, components *c, void *params, int k) {
    parameters *p = (parameters *)params;
	mpfr_mul(c->x, p->a, x[k], RND);
	mpfr_mul(c->y, p->b, y[k], RND);
	mpfr_mul(c->z, p->c, z[k], RND);
}

static result skip (char* name) {
    total++;
    skipped++;
    printf("%sSKIPPED%s %s\n", KYLW, KNRM, name);
    return SKIPPED;
}

static result compare (char* name, series a, series b, int n) {
    total++;
    for (int k = 0; k < n; k++) {
        mpfr_sub(delta, a[k], b[k], RND);
        if (mpfr_cmp_abs(delta, tolerance) > 0) {
            printf("%s FAILED%s %s  k: %d  LHS: %.6e  RHS: %.6e  diff %.3e\n",
                    KRED, KNRM, name, k, mpfr_get_d(a[k], RND), mpfr_get_d(b[k], RND), mpfr_get_d(delta, RND));
            return FAILED;
        }
        if (debug >= 2 && k == 0) printf("\n");
        if (debug >= 2) printf("%s  DEBUG%s  k: %2d  %+.6e %+.6e  diff %+.3e\n",
                               KNRM, KNRM, k, mpfr_get_d(a[k], RND), mpfr_get_d(b[k], RND), mpfr_get_d(delta, RND));
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

    series abs = t_jet(n);
    series scale = t_jet(n);
    series sum = t_jet(n);
    series diff = t_jet(n);
    series prod = t_jet(n);
    series quot = t_jet(n);
    series inv = t_jet(n);
    series inv_x = t_jet(n);
    series neg = t_jet(n);
    series sqr_x = t_jet(n);
    series sqr_sin_x = t_jet(n);
    series sqr_cos_x = t_jet(n);
    series sqr_tan_x = t_jet(n);
    series sqrt_x = t_jet(n);
    series pow = t_jet(n);
    series exp_x = t_jet(n);
    series neg_exp_x = t_jet(n);
    series ln = t_jet(n);
    series ln_x = t_jet(n);
    series sin = t_jet(n);
    series sin_2x = t_jet(n);
    series cos = t_jet(n);
    series cos_2x = t_jet(n);
    series tan = t_jet(n);
    series sec2 = t_jet(n);

    series c1 = ad_series_c(n, D1);
    series x = ad_series_v(n, x0);

    mpfr_init(PI_2);
    mpfr_const_pi(PI_2, RND);
    mpfr_div_2ui(PI_2, PI_2, 1, RND);

    int x_positive = mpfr_sgn(x[0]) > 0;
    int x_non_zero = mpfr_zero_p(x[0]) == 0;
    int x_lt_pi_2 = mpfr_cmpabs(x[0], PI_2) < 0;

    tsm(argc, argv, 10, D01, 10, D1, D1, D1);

    printf("\n");
    ad_sqr(sqr_x, x, n);
    if (x_non_zero) ad_inv(inv_x, x, n);
    if (x_positive) ad_sqrt(sqrt_x, x, n);
    char* name = "x * x == sqr(x)";
    compare(name, ad_prod(prod, x, x, n), sqr_x, n);
    name = "sqr(x) / x == x";
    x_non_zero ? compare(name, ad_quot(quot, sqr_x, x, n), x, n) : skip(name);
    name = "x * 1 / x == 1";
    x_non_zero ? compare(name, ad_prod(prod, x, inv_x, n), c1, n) : skip(name);
    name = "sqrt(x) * sqrt(x) == x";
    x_positive ? compare(name, ad_prod(prod, sqrt_x, sqrt_x, n), x, n) : skip(name);
    name = "x / sqrt(x) == sqrt(x)";
    x_positive ? compare(name, ad_quot(quot, x, sqrt_x, n), sqrt_x, n) : skip(name);

    name = "x^2 == sqr(x)";
    x_positive ? compare(name, ad_pwr(pow, x, D2, n), sqr_x, n) : skip(name);
    name = "x^1 == x";
    x_positive ? compare(name, ad_pwr(pow, x, D1, n), x, n) : skip(name);
    name = "x^0.5 == sqrt(x)";
    x_positive ? compare(name, ad_pwr(pow, x, D05, n), sqrt_x, n): skip(name);
    name = "x^0 == 1";
    x_positive ? compare(name, ad_pwr(pow, x, D0, n), c1, n) : skip(name);
    name = "x^-0.5 == 1 / sqrt(x)";
    x_positive ? compare(name, ad_pwr(pow, x, D_05, n), ad_inv(inv, sqrt_x, n), n) : skip(name);
    name = "x^-1 == 1 / x";
    x_positive ? compare(name, ad_pwr(pow, x, D_1, n), inv_x, n) : skip(name);
    name = "x^-2 == 1 / sqr(x)";
    x_positive ? compare(name, ad_pwr(pow, x, D_2, n), ad_inv(inv, sqr_x, n), n) : skip(name);

    name = "sqr(x) * x^-3 == 1 / x";
    x_positive ? compare(name, ad_prod(prod, sqr_x, ad_pwr(pow, x, D_3, n), n), inv_x, n) : skip(name);

    name = "sqr(x)^0.5 == |x|";
    x_non_zero ? compare(name, ad_pwr(pow, sqr_x, D05, n), ad_abs(abs, x, n), n) : skip(name);

    name = "log(e^x) == x";
    compare(name, ad_ln(ln, ad_exp(exp_x, x, n), n), x, n);
    name = "log(sqr(x)) == 2 * log(x)";
    if (x_positive) ad_ln(ln_x, x, n);
    x_positive ? compare(name, ad_ln(ln, sqr_x, n), ad_scale(scale, ln_x, D2, n), n) : skip(name);
    name = "log(sqrt(x)) == 0.5 * log(x)";
    x_positive ? compare(name, ad_ln(ln, sqrt_x, n), ad_scale(scale, ln_x, D05, n), n) : skip(name);
    name = "log(1 / x) == - log(x)";
    x_positive ? compare(name, ad_ln(ln, inv_x, n), ad_neg(neg, ln_x, n), n) : skip(name);
    name = "log(x^-3) == - 3 * log(x)";
    x_positive ? compare(name, ad_ln(ln, ad_pwr(pow, x, D_3, n), n), ad_scale(scale, ln_x, D_3, n), n) : skip(name);

    ad_sin_cos(sin, cos, x, HYP, n);
    ad_tan_sec2(tan, sec2, x, HYP, n);
    ad_sqr(sqr_sin_x, sin, n);
    ad_sqr(sqr_cos_x, cos, n);
    ad_sin_cos(sin_2x, cos_2x, ad_scale(scale, x, D2, n), HYP, n);
    name = "cosh^2(x) - sinh^2(x) == 1";
    compare(name, ad_minus(diff, sqr_cos_x, sqr_sin_x, n), c1, n);
    name = "sech^2(x) + tanh^2(x) == 1";
    compare(name, ad_plus(sum, sec2, ad_sqr(sqr_tan_x, tan, n), n), c1, n);
    name = "tanh(x) == sinh(x) / cosh(x)";
    compare(name, tan, ad_quot(quot, sin, cos, n), n);
    name = "sech^2(x) == 1 / cosh^2(x)";
    compare(name, sec2, ad_inv(inv, sqr_cos_x, n), n);
    name = "sinh(2x) == 2 * sinh(x) * cosh(x)";
    compare(name, sin_2x, ad_scale(scale, ad_prod(prod, sin, cos, n), D2, n), n);
    name = "cosh(2x) == cosh^2(x) + sinh^2(x)";
    compare(name, cos_2x, ad_plus(sum, sqr_cos_x, sqr_sin_x, n), n);

    ad_exp(exp_x, x, n);
    ad_exp(neg_exp_x, ad_neg(neg, x, n), n);
    name = "sinh(x) == 0.5 * (e^x - e^-x)";
    compare(name, sin, ad_scale(scale, ad_minus(diff, exp_x, neg_exp_x, n), D05, n), n);
    name = "cosh(x) == 0.5 * (e^x + e^-x)";
    compare(name, cos, ad_scale(scale, ad_plus(sum, exp_x, neg_exp_x, n), D05, n), n);

    ad_sin_cos(sin, cos, x, TRIG, n);
    ad_tan_sec2(tan, sec2, x, TRIG, n);
    ad_sqr(sqr_sin_x, sin, n);
    ad_sqr(sqr_cos_x, cos, n);
    ad_sin_cos(sin_2x, cos_2x, ad_scale(scale, x, D2, n), TRIG, n);
    name = "cos^2(x) + sin^2(x) == 1";
    compare(name, ad_plus(sum, sqr_cos_x, sqr_sin_x, n), c1, n);
    name = "sec^2(x) - tan^2(x) == 1";
    x_lt_pi_2 ? compare(name, ad_minus(diff, sec2, ad_sqr(sqr_tan_x, tan, n), n), c1, n) : skip(name);
    name = "tan(x) == sin(x) / cos(x)";
    x_lt_pi_2 ? compare(name, tan, ad_quot(quot, sin, cos, n), n) : skip(name);
    name = "sec^2(x) == 1 / cos^2(x)";
    x_lt_pi_2 ? compare(name, sec2, ad_inv(inv, sqr_cos_x, n), n) : skip(name);
    name = "sin(2x) == 2 * sin(x) * cos(x)";
    compare(name, sin_2x, ad_scale(scale, ad_prod(prod, sin, cos, n), D2, n), n);
    name = "cos(2x) == cos^2(x) - sin^2(x)";
    compare(name, cos_2x, ad_minus(diff, sqr_cos_x, sqr_sin_x, n), n);

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

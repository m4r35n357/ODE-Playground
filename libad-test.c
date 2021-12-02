/*
 * Automatic Differentiation of Taylor Series, newest validation checks
 *
 * Example: ./libad-test-dbg 32 20 1 1e-21 [ 0 | 1 | 2 ]
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
#define KWHT "\x1B[1;37m"
#define KGRN "\x1B[1;32m"
#define KYLW "\x1B[1;33m"
#define KRED "\x1B[1;31m"

static int n, debug = 0, total = 0, passed = 0, skipped = 0;

static mpfr_t delta, tolerance, D0, D01, D05, D_05, D1, D_1, D2, D_2, D3, D_3;

static void libad_test_init (void) {
    ad_init(n);
    t_init(12);
    mpfr_init(delta);
    mpfr_init_set_si(D0, 0, RND);
    mpfr_init_set_str(D01, "0.1", BASE, RND);
    mpfr_init_set_str(D05, "0.5", BASE, RND);
    mpfr_init_set_str(D_05, "-0.5", BASE, RND);
    mpfr_init_set_si(D1, 1, RND);
    mpfr_init_set_si(D2, 2, RND);
    mpfr_init_set_si(D3, 3, RND);
    mpfr_init_set_si(D_1, -1, RND);
    mpfr_init_set_si(D_2, -2, RND);
    mpfr_init_set_si(D_3, -3, RND);
}

typedef struct { mpfr_t a, b, c; } parameters;

void *get_p (int argc, char **argv, int order) {
    (void)argc; (void)argv; (void)order;
    parameters *p = malloc(sizeof (parameters));
    mpfr_init_set(p->a, D1, RND);
    mpfr_init_set(p->b, D0, RND);
    mpfr_init_set(p->c, D_1, RND);
    return p;
}

void ode (components *v, series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    mpfr_mul(v->x, p->a, x[k], RND);
    mpfr_mul(v->y, p->b, y[k], RND);
    mpfr_mul(v->z, p->c, z[k], RND);
}

static void skip (char* name) {
    total++;
    skipped++;
    if (debug >= 1) fprintf(stderr, "%sSKIPPED%s %s\n", KYLW, KNRM, name);
}

static void compare (char* name, series a, series b) {
    total++;
    for (int k = 0; k < n; k++) {
        mpfr_sub(delta, a[k], b[k], RND);
        if (mpfr_number_p(delta) == 0 || mpfr_cmpabs(delta, tolerance) > 0) {
            fprintf(stderr, "%s FAILED%s %s  k: %d  LHS: %.6e  RHS: %.6e  diff %.3e\n",
                    KRED, KNRM, name, k, mpfr_get_d(a[k], RND), mpfr_get_d(b[k], RND), mpfr_get_d(delta, RND));
            return;
        }
        if (debug >= 2) {
            if (k == 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s  DEBUG%s  k: %2d  %+.6e %+.6e  diff %+.3e\n",
                    KNRM, KNRM, k, mpfr_get_d(a[k], RND), mpfr_get_d(b[k], RND), mpfr_get_d(delta, RND));
        }
    }
    if (debug >= 1) fprintf(stderr, "%s PASSED%s %s\n", KGRN, KNRM, name);
    passed++;
}

static series contaminate (series a) {
    for (int k = 0; k < n; k++) {
        mpfr_set_inf(a[k], 1);
    }
    return a;
}

int main (int argc, char **argv) {
    mpfr_t PI_2;

    assert(argc == 5 || argc == 6);
    long double precision = strtod(argv[1], NULL) * 3.322L;
    mpfr_set_default_prec((int)precision);
    n = (int)strtol(argv[2], NULL, BASE);
    assert(n > 1);
    libad_test_init();
    series x = t_jet(n + 1);
    mpfr_init_set_str(x[0], argv[3], BASE, RND);
    for (int k = 1; k <= n; k++) {
        mpfr_div_si(x[k], x[0], k * k, RND);
    }
    mpfr_init_set_str(tolerance, argv[4], BASE, RND);
    if (argc == 6) debug = (int)strtol(argv[5], NULL, BASE);

    series c1 = t_jet(n);
    mpfr_set(c1[0], D1, RND);
    for (int i = 1; i < n; i++) {
        mpfr_set_zero(c1[i], 1);
    }

    mpfr_init(PI_2);
    mpfr_const_pi(PI_2, RND);
    mpfr_div_2si(PI_2, PI_2, 1, RND);

    int x_positive = mpfr_sgn(x[0]) > 0;
    int x_non_zero = mpfr_zero_p(x[0]) == 0;
    int x_lt_pi_2 = mpfr_cmpabs(x[0], PI_2) < 0;

    series abs = contaminate(t_jet(n));
    series scale = contaminate(t_jet(n));
    series sum = contaminate(t_jet(n));
    series diff = contaminate(t_jet(n));
    series prod = contaminate(t_jet(n));
    series quot = contaminate(t_jet(n));
    series inv = contaminate(t_jet(n));
    series inv_x = contaminate(t_jet(n));
    series sqr_x = contaminate(t_jet(n));
    series sqr_sin_x = contaminate(t_jet(n));
    series sqr_cos_x = contaminate(t_jet(n));
    series sqr_tan_x = contaminate(t_jet(n));
    series sqrt_x = contaminate(t_jet(n));
    series pow = contaminate(t_jet(n));
    series exp_x = contaminate(t_jet(n));
    series neg_exp_x = contaminate(t_jet(n));
    series ln = contaminate(t_jet(n));
    series ln_x = contaminate(t_jet(n));
    series sin = contaminate(t_jet(n));
    series sin_2x = contaminate(t_jet(n));
    series cos = contaminate(t_jet(n));
    series cos_2x = contaminate(t_jet(n));
    series tan = contaminate(t_jet(n));
    series sec2 = contaminate(t_jet(n));

    fprintf(stdout, "\n");
    fprintf(stdout, "Horner\n");
    series p = t_jet(n >= 8 ? n : 8);
    mpfr_set_si(p[0], 1, RND);
    mpfr_set_si(p[1], 3, RND);
    mpfr_set_si(p[2], 0, RND);
    mpfr_set_si(p[3], 2, RND);
    mpfr_fprintf(stdout, " 23 %8.3RNf\n", *t_horner(p, 3, D2));
    mpfr_set_si(p[0], 3, RND);
    mpfr_set_si(p[1], -1, RND);
    mpfr_set_si(p[2], 2, RND);
    mpfr_set_si(p[3], -4, RND);
    mpfr_set_si(p[4], 0, RND);
    mpfr_set_si(p[5], 1, RND);
    mpfr_fprintf(stdout, "153 %8.3RNf\n", *t_horner(p, 5, D3));
    mpfr_set_si(p[0], 1, RND);
    mpfr_set_si(p[1], -4, RND);
    mpfr_set_si(p[2], 0, RND);
    mpfr_set_si(p[3], 0, RND);
    mpfr_set_si(p[4], 2, RND);
    mpfr_set_si(p[5], 3, RND);
    mpfr_set_si(p[6], 0, RND);
    mpfr_set_si(p[7], -2, RND);
    mpfr_fprintf(stdout, "201 %8.3RNf\n", *t_horner(p, 7, D_2));

    fprintf(stdout, "\n");
    fprintf(stdout, "TSM\n");
    tsm(argc, argv, n, D01, 10, D1, D1, D1);
    mpfr_t e1, e0, e_1;
    mpfr_inits(e1, e0, e_1, NULL);
    mpfr_exp(e1, D1, RND);
    mpfr_exp(e0, D0, RND);
    mpfr_exp(e_1, D_1, RND);
    fprintf(stdout, "Check\n");
    t_output(e1, e0, e_1, D01, 10);

    fprintf(stderr, "\n");
    ad_sqr(sqr_x, x);
    if (x_non_zero) ad_inv(inv_x, x);
    if (x_positive) ad_sqrt(sqrt_x, x);

    char* name = "x * x == sqr(x)";
    compare(name, ad_mul(prod, x, x), sqr_x);

    name = "sqr(x) / x == x";
    x_non_zero ? compare(name, ad_div(quot, sqr_x, x), x) : skip(name);

    name = "x * 1 / x == 1";
    x_non_zero ? compare(name, ad_mul(prod, x, inv_x), c1) : skip(name);

    name = "sqrt(x) * sqrt(x) == x";
    x_positive ? compare(name, ad_mul(prod, sqrt_x, sqrt_x), x) : skip(name);

    name = "x / sqrt(x) == sqrt(x)";
    x_positive ? compare(name, ad_div(quot, x, sqrt_x), sqrt_x) : skip(name);

    if (debug != 0) fprintf(stderr, "\n");

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

    if (debug != 0) fprintf(stderr, "\n");

    name = "sqr(x) * x^-3 == 1 / x";
    x_positive ? compare(name, ad_mul(prod, sqr_x, ad_pwr(pow, x, D_3)), inv_x) : skip(name);

    name = "sqr(x)^0.5 == |x|";
    x_non_zero ? compare(name, ad_pwr(pow, sqr_x, D05), ad_abs(abs, x)) : skip(name);

    if (debug != 0) fprintf(stderr, "\n");
    ad_exp(exp_x, x);
    if (x_positive) ad_ln(ln_x, x);

    name = "log(e^x) == x";
    compare(name, ad_ln(ln, exp_x), x);

    name = "log(sqr(x)) == log(x) * 2";
    x_positive ? compare(name, ad_ln(ln, sqr_x), ad_scale(scale, ln_x, D2)) : skip(name);

    name = "log(sqrt(x)) == log(x) / 2";
    x_positive ? compare(name, ad_ln(ln, sqrt_x), ad_scale(scale, ln_x, D05)) : skip(name);

    name = "log(1 / x) == - log(x)";
    x_positive ? compare(name, ad_ln(ln, inv_x), ad_scale(scale, ln_x, D_1)) : skip(name);

    name = "log(x^-3) == - 3 * log(x)";
    x_positive ? compare(name, ad_ln(ln, ad_pwr(pow, x, D_3)), ad_scale(scale, ln_x, D_3)) : skip(name);

    if (debug != 0) fprintf(stderr, "\n");
    ad_sin_cos(sin, cos, x, HYP);
    ad_tan_sec2(tan, sec2, x, HYP);
    ad_sqr(sqr_sin_x, sin);
    ad_sqr(sqr_cos_x, cos);
    ad_sin_cos(sin_2x, cos_2x, ad_scale(scale, x, D2), HYP);

    name = "cosh^2(x) - sinh^2(x) == 1";
    compare(name, ad_sub(diff, sqr_cos_x, sqr_sin_x), c1);

    name = "sech^2(x) + tanh^2(x) == 1";
    compare(name, ad_add(sum, sec2, ad_sqr(sqr_tan_x, tan)), c1);

    name = "tanh(x) == sinh(x) / cosh(x)";
    compare(name, tan, ad_div(quot, sin, cos));

    name = "sech^2(x) == 1 / cosh^2(x)";
    compare(name, sec2, ad_inv(inv, sqr_cos_x));

    name = "sinh(2x) == 2 * sinh(x) * cosh(x)";
    compare(name, sin_2x, ad_scale(scale, ad_mul(prod, sin, cos), D2));

    name = "cosh(2x) == cosh^2(x) + sinh^2(x)";
    compare(name, cos_2x, ad_add(sum, sqr_cos_x, sqr_sin_x));

    if (debug != 0) fprintf(stderr, "\n");
    ad_exp(neg_exp_x, ad_scale(scale, x, D_1));

    name = "cosh(x) == (e^x + e^-x) / 2";
    compare(name, cos, ad_scale(scale, ad_add(sum, exp_x, neg_exp_x), D05));

    name = "sinh(x) == (e^x - e^-x) / 2";
    compare(name, sin, ad_scale(scale, ad_sub(diff, exp_x, neg_exp_x), D05));

    if (debug != 0) fprintf(stderr, "\n");
    ad_sin_cos(sin, cos, x, TRIG);
    ad_tan_sec2(tan, sec2, x, TRIG);
    ad_sqr(sqr_sin_x, sin);
    ad_sqr(sqr_cos_x, cos);
    ad_sin_cos(sin_2x, cos_2x, ad_scale(scale, x, D2), TRIG);

    name = "cos^2(x) + sin^2(x) == 1";
    compare(name, ad_add(sum, sqr_cos_x, sqr_sin_x), c1);

    name = "sec^2(x) - tan^2(x) == 1";
    x_lt_pi_2 ? compare(name, ad_sub(diff, sec2, ad_sqr(sqr_tan_x, tan)), c1) : skip(name);

    name = "tan(x) == sin(x) / cos(x)";
    x_lt_pi_2 ? compare(name, tan, ad_div(quot, sin, cos)) : skip(name);

    name = "sec^2(x) == 1 / cos^2(x)";
    x_lt_pi_2 ? compare(name, sec2, ad_inv(inv, sqr_cos_x)) : skip(name);

    name = "sin(2x) == 2 * sin(x) * cos(x)";
    compare(name, sin_2x, ad_scale(scale, ad_mul(prod, sin, cos), D2));

    name = "cos(2x) == cos^2(x) - sin^2(x)";
    compare(name, cos_2x, ad_sub(diff, sqr_cos_x, sqr_sin_x));

    if (debug != 0) fprintf(stderr, "\n");
    fprintf(stderr, "%sTotal%s: %d, %sPASSED%s %d", KWHT, KNRM, total, KGRN, KNRM, passed);
    if (skipped > 0) {
        fprintf(stderr, ", %sSKIPPED%s %d", KYLW, KNRM, skipped);
    }
    if (passed == total - skipped) {
        fprintf(stderr, "\n\n");
        return 0;
    } else {
        fprintf(stderr, ", %sFAILED%s %d\n\n", KRED, KNRM, total - passed - skipped);
        return 1;
    }
}

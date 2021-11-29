/*
 * Automatic Differentiation of Taylor Series, newest validation checks
 *
 * Example: ./libad-test-dbg 20 1 1e-18 [ 0 | 1 | 2 ]
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "taylor-ode.h"
#include "ad.h"

#define KNRM "\x1B[0;37m"
#define KWHT "\x1B[1;37m"
#define KGRN "\x1B[1;32m"
#define KYLW "\x1B[1;33m"
#define KRED "\x1B[1;31m"

static int n, debug = 0, total = 0, passed = 0, skipped = 0;

static real delta, tolerance;

static const real PLUS1 = 1.0L, ZERO = 0.0L, MINUS1 = -1.0L;

typedef struct { real a, b, c; } parameters;

void *get_p (int argc, char **argv, int order) {
    (void)argc; (void)argv; (void)order;
    parameters *p = malloc(sizeof (parameters));
    p->a = PLUS1;
    p->b = ZERO;
    p->c = MINUS1;
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = p->a * x[k],
        .y = p->b * y[k],
        .z = p->c * z[k]
    };
}

static void skip (char* name) {
    total++;
    skipped++;
    if (debug >= 1) fprintf(stderr, "%sSKIPPED%s %s\n", KYLW, KNRM, name);
}

static void compare (char* name, series a, series b) {
    total++;
    for (int k = 0; k < n; k++) {
        delta = a[k] - b[k];
        if (isnan(delta) || isinf(delta) || fabsl(delta) > tolerance) {
            fprintf(stderr, "%s FAILED%s %s  k: %d  LHS: %.6Le  RHS: %.6Le  diff %.3Le\n", KRED, KNRM, name, k, a[k], b[k], delta);
            return;
        }
        if (debug >= 2) {
            if (k == 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s  DEBUG%s  k: %2d  %+.6Le %+.6Le  diff %+.3Le\n", KNRM, KNRM, k, a[k], b[k], delta);
        }
    }
    if (debug >= 1) fprintf(stderr, "%s PASSED%s %s\n", KGRN, KNRM, name);
    passed++;
}

static series contaminate (series a) {
    for (int k = 0; k < n; k++) {
        a[k] = INFINITY;
    }
    return a;
}

int main (int argc, char **argv) {
    real PI_2 = 0.5L * MY_PI;

    assert(argc == 4 || argc == 5);
    n = (int)strtol(argv[1], NULL, BASE);
    assert(n > 1);
    ad_init(n);
    series x = t_jet(n + 1);
    x[0] = strtold(argv[2], NULL);
    for (int k = 1; k <= n; k++) {
        x[k] = x[0] / (k * k);
    }
    tolerance = strtold(argv[3], NULL);
    if (argc == 5) debug = (int)strtol(argv[4], NULL, BASE);

    series c1 = t_jet(n);
    c1[0] = 1.0L;

    int x_positive = x[0] > 0.0L;
    int x_non_zero = x[0] != 0.0L;
    int x_lt_pi_2 = x[0] < PI_2;

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
    p[0] = 1.0L; p[1] = 3.0L; p[2] = 0.0L; p[3] = 2.0L;
    fprintf(stdout, " 23 %8.3Lf\n", t_horner(p, 3, 2.0L));
    p[0] = 3; p[1] = -1.0L; p[2] = 2.0L; p[3] = -4.0L; p[4] = 0.0L; p[5] = 1.0L;
    fprintf(stdout, "153 %8.3Lf\n", t_horner(p, 5, 3.0L));
    p[0] = 1.0L; p[1] = -4.0L; p[2] = 0.0L; p[3] = 0.0L; p[4] = 2.0L; p[5] = 3.0L; p[6] = 0.0L; p[7] = -2.0L;
    fprintf(stdout, "201 %8.3Lf\n", t_horner(p, 7, -2.0L));

    fprintf(stdout, "\n");
    fprintf(stdout, "TSM\n");
    int dp = 12, steps = 10;
    real step = 0.1L;
    tsm(argc, argv, dp, n, step, steps, 1.0L, 1.0L, 1.0L);
    fprintf(stdout, "Check\n");
    t_output(dp, expl(PLUS1), expl(ZERO), expl(MINUS1), step * steps, "_", "_", "_");

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
    x_positive ? compare(name, ad_pwr(pow, x, 2.0L), sqr_x) : skip(name);

    name = "x^1 == x";
    x_positive ? compare(name, ad_pwr(pow, x, 1.0L), x) : skip(name);

    name = "x^0.5 == sqrt(x)";
    x_positive ? compare(name, ad_pwr(pow, x, 0.5L), sqrt_x): skip(name);

    name = "x^0 == 1";
    x_positive ? compare(name, ad_pwr(pow, x, 0.0L), c1) : skip(name);

    name = "x^-0.5 == 1 / sqrt(x)";
    x_positive ? compare(name, ad_pwr(pow, x, -0.5L), ad_inv(inv, sqrt_x)) : skip(name);

    name = "x^-1 == 1 / x";
    x_positive ? compare(name, ad_pwr(pow, x, -1.0L), inv_x) : skip(name);

    name = "x^-2 == 1 / sqr(x)";
    x_positive ? compare(name, ad_pwr(pow, x, -2.0L), ad_inv(inv, sqr_x)) : skip(name);

    if (debug != 0) fprintf(stderr, "\n");

    name = "sqr(x) * x^-3 == 1 / x";
    x_positive ? compare(name, ad_mul(prod, sqr_x, ad_pwr(pow, x, -3.0L)), inv_x) : skip(name);

    name = "sqr(x)^0.5 == |x|";
    x_non_zero ? compare(name, ad_pwr(pow, sqr_x, 0.5L), ad_abs(abs, x)) : skip(name);

    if (debug != 0) fprintf(stderr, "\n");
    ad_exp(exp_x, x);
    if (x_positive) ad_ln(ln_x, x);

    name = "log(e^x) == x";
    compare(name, ad_ln(ln, exp_x), x);

    name = "log(sqr(x)) == log(x) * 2";
    x_positive ? compare(name, ad_ln(ln, sqr_x), ad_scale(scale, ln_x, 2.0L)) : skip(name);

    name = "log(sqrt(x)) == log(x) / 2";
    x_positive ? compare(name, ad_ln(ln, sqrt_x), ad_scale(scale, ln_x, 0.5L)) : skip(name);

    name = "log(1 / x) == - log(x)";
    x_positive ? compare(name, ad_ln(ln, inv_x), ad_scale(scale, ln_x, -1.0L)) : skip(name);

    name = "log(x^-3) == - 3 * log(x)";
    x_positive ? compare(name, ad_ln(ln, ad_pwr(pow, x, -3.0L)), ad_scale(scale, ln_x, -3.0L)) : skip(name);

    if (debug != 0) fprintf(stderr, "\n");
    ad_sin_cos(sin, cos, x, HYP);
    ad_tan_sec2(tan, sec2, x, HYP);
    ad_sqr(sqr_sin_x, sin);
    ad_sqr(sqr_cos_x, cos);
    ad_sin_cos(sin_2x, cos_2x, ad_scale(scale, x, 2.0L), HYP);

    name = "cosh^2(x) - sinh^2(x) == 1";
    compare(name, ad_sub(diff, sqr_cos_x, sqr_sin_x), c1);

    name = "sech^2(x) + tanh^2(x) == 1";
    compare(name, ad_add(sum, sec2, ad_sqr(sqr_tan_x, tan)), c1);

    name = "tanh(x) == sinh(x) / cosh(x)";
    compare(name, tan, ad_div(quot, sin, cos));

    name = "sech^2(x) == 1 / cosh^2(x)";
    compare(name, sec2, ad_inv(inv, sqr_cos_x));

    name = "sinh(2x) == 2 * sinh(x) * cosh(x)";
    compare(name, sin_2x, ad_scale(scale, ad_mul(prod, sin, cos), 2.0L));

    name = "cosh(2x) == cosh^2(x) + sinh^2(x)";
    compare(name, cos_2x, ad_add(sum, sqr_cos_x, sqr_sin_x));

    if (debug != 0) fprintf(stderr, "\n");
    ad_exp(neg_exp_x, ad_scale(scale, x, -1.0L));

    name = "cosh(x) == (e^x + e^-x) / 2";
    compare(name, cos, ad_scale(scale, ad_add(sum, exp_x, neg_exp_x), 0.5L));

    name = "sinh(x) == (e^x - e^-x) / 2";
    compare(name, sin, ad_scale(scale, ad_sub(diff, exp_x, neg_exp_x), 0.5L));

    if (debug != 0) fprintf(stderr, "\n");
    ad_sin_cos(sin, cos, x, TRIG);
    ad_tan_sec2(tan, sec2, x, TRIG);
    ad_sqr(sqr_sin_x, sin);
    ad_sqr(sqr_cos_x, cos);
    ad_sin_cos(sin_2x, cos_2x, ad_scale(scale, x, 2.0L), TRIG);

    name = "cos^2(x) + sin^2(x) == 1";
    compare(name, ad_add(sum, sqr_cos_x, sqr_sin_x), c1);

    name = "sec^2(x) - tan^2(x) == 1";
    x_lt_pi_2 ? compare(name, ad_sub(diff, sec2, ad_sqr(sqr_tan_x, tan)), c1) : skip(name);

    name = "tan(x) == sin(x) / cos(x)";
    x_lt_pi_2 ? compare(name, tan, ad_div(quot, sin, cos)) : skip(name);

    name = "sec^2(x) == 1 / cos^2(x)";
    x_lt_pi_2 ? compare(name, sec2, ad_inv(inv, sqr_cos_x)) : skip(name);

    name = "sin(2x) == 2 * sin(x) * cos(x)";
    compare(name, sin_2x, ad_scale(scale, ad_mul(prod, sin, cos), 2.0L));

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

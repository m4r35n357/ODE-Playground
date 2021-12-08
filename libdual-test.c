/*
 * Dual numbers, newest validation checks
 *
 * Example: ./libdual-test-dbg 1 1e-18 [ 0 | 1 | 2 ]
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "dual.h"

#define KNRM "\x1B[0;37m"
#define KWHT "\x1B[1;37m"
#define KGRN "\x1B[1;32m"
#define KYLW "\x1B[1;33m"
#define KRED "\x1B[1;31m"

static int debug = 0, total = 0, passed = 0, skipped = 0;

static real x0, delta_val, delta_dot, tolerance;

static void skip (char* name) {
    total++;
    skipped++;
    if (debug >= 1) fprintf(stderr, "%s SKIP%s %s\n", KYLW, KNRM, name);
}

static void compare (char* name, dual a, dual b) {
    total++;
    delta_val = a.val - b.val;
    if (fabsl(delta_val) > tolerance) {
        printf("%s FAIL%s %s  LHS: %.6Le  RHS: %.6Le  diff %.3Le\n", KRED, KNRM, name, a.val, b.val, delta_val);
        return;
    }
    delta_dot = a.dot - b.dot;
    if (fabsl(delta_dot) > tolerance) {
        printf("%s FAIL%s %s  LHS: %.6Le  RHS: %.6Le  diff %.3Le\n", KRED, KNRM, name, a.dot, b.dot, delta_dot);
        return;
    }
    if (debug >= 2) printf("\n");
    if (debug >= 2) printf("%s  DEBUG%s  %+.6Le %+.6Le  diff %+.3Le\n", KNRM, KNRM, a.val, b.val, delta_val);
    if (debug >= 2) printf("%s  DEBUG%s  %+.6Le %+.6Le  diff %+.3Le\n", KNRM, KNRM, a.dot, b.dot, delta_dot);
    if (debug >= 1) printf("%s PASS%s %s\n", KGRN, KNRM, name);
    passed++;
}

int main (int argc, char **argv) {
    real PI_2 = 0.5L * MY_PI;

    assert(argc == 3 || argc == 4);
    //ad_init(n);
    x0 = strtold(argv[1], NULL);
    tolerance = strtold(argv[2], NULL);
    if (argc == 4) debug = (int)strtol(argv[3], NULL, 10);

    dual inv_x = d_dual(0.0L);
    dual sqr_x = d_dual(0.0L);
    dual sqr_sin_x = d_dual(0.0L);
    dual sqr_cos_x = d_dual(0.0L);
    dual sqrt_x = d_dual(0.0L);
    dual exp_x = d_dual(0.0L);
    dual neg_exp_x = d_dual(0.0L);
    dual ln_x = d_dual(0.0L);
    dual sin = d_dual(0.0L);
    dual sin_2x = d_dual(0.0L);
    dual cos = d_dual(0.0L);
    dual cos_2x = d_dual(0.0L);
    dual tan = d_dual(0.0L);

    dual c1 = d_dual(1.0L);
    dual x = d_var(x0);
    dual xpx = d_scale(x, 2.0L);

    int x_positive = x.val > 0.0L;
    int x_non_zero = x.val != 0.0L;
    int x_lt_pi_2 = x.val < PI_2;

    printf("\n");
    sqr_x = d_sqr(x);
    if (x_non_zero) inv_x = d_inv(x);
    if (x_positive) sqrt_x = d_sqrt(x);
    char* name = "x * x == sqr(x)";
    compare(name, d_mul(x, x), sqr_x);

    name = "sqr(x) / x == x";
    x_non_zero ? compare(name, d_div(sqr_x, x), x) : skip(name);

    name = "x * 1 / x == 1";
    x_non_zero ? compare(name, d_mul(x, inv_x), c1) : skip(name);

    name = "sqrt(x) * sqrt(x) == x";
    x_positive ? compare(name, d_mul(sqrt_x, sqrt_x), x) : skip(name);

    name = "x / sqrt(x) == sqrt(x)";
    x_positive ? compare(name, d_div(x, sqrt_x), sqrt_x) : skip(name);

    if (debug != 0) printf("\n");

    name = "x^2 == sqr(x)";
    x_positive ? compare(name, d_pow(x, 2.0L), sqr_x) : skip(name);

    name = "x^1 == x";
    x_positive ? compare(name, d_pow(x, 1.0L), x) : skip(name);

    name = "x^0.5 == sqrt(x)";
    x_positive ? compare(name, d_pow(x, 0.5L), sqrt_x): skip(name);

    name = "x^0 == 1";
    x_positive ? compare(name, d_pow(x, 0.0L), c1) : skip(name);

    name = "x^-0.5 == 1 / sqrt(x)";
    x_positive ? compare(name, d_pow(x, -0.5L), d_inv(sqrt_x)) : skip(name);

    name = "x^-1 == 1 / x";
    x_positive ? compare(name, d_pow(x, -1.0L), inv_x) : skip(name);

    name = "x^-2 == 1 / sqr(x)";
    x_positive ? compare(name, d_pow(x, -2.0L), d_inv(sqr_x)) : skip(name);

    if (debug != 0) printf("\n");

    name = "sqr(x) * x^-3 == 1 / x";
    x_positive ? compare(name, d_mul(sqr_x, d_pow(x, -3.0L)), inv_x) : skip(name);

    name = "sqr(x)^0.5 == |x|";
    x_non_zero ? compare(name, d_pow(sqr_x, 0.5L), d_abs(x)) : skip(name);

    if (debug != 0) printf("\n");
    if (x_positive) ln_x = d_log(x);

    name = "log(e^x) == x";
    compare(name, d_log(d_exp(x)), x);

    name = "log(sqr(x)) == log(x) * 2";
    if (x_positive) d_log(x);
    x_positive ? compare(name, d_log(sqr_x), d_scale(ln_x, 2.0L)) : skip(name);

    name = "log(sqrt(x)) == log(x) / 2";
    x_positive ? compare(name, d_log(sqrt_x), d_scale(ln_x, 0.5L)) : skip(name);

    name = "log(1 / x) == - log(x)";
    x_positive ? compare(name, d_log(inv_x), d_scale(ln_x, -1.0L)) : skip(name);

    name = "log(x^-3) == - 3 * log(x)";
    x_positive ? compare(name, d_log(d_pow(x, -3.0L)), d_scale(ln_x, -3.0L)) : skip(name);

    if (debug != 0) printf("\n");
    sin = d_sinh(x);
    cos = d_cosh(x);
    tan = d_tanh(x);
    sqr_sin_x = d_sqr(sin);
    sqr_cos_x = d_sqr(cos);
    sin_2x = d_sinh(xpx);
    cos_2x = d_cosh(xpx);

    name = "cosh^2(x) - sinh^2(x) == 1";
    compare(name, d_sub(sqr_cos_x, sqr_sin_x), c1);

    name = "tanh(x) == sinh(x) / cosh(x)";
    compare(name, tan, d_div(sin, cos));

    name = "sinh(2x) == 2 * sinh(x) * cosh(x)";
    compare(name, sin_2x, d_scale(d_mul(sin, cos), 2.0L));

    name = "cosh(2x) == cosh^2(x) + sinh^2(x)";
    compare(name, cos_2x, d_add(sqr_cos_x, sqr_sin_x));

    if (debug != 0) printf("\n");
    exp_x = d_exp(x);
    neg_exp_x = d_exp(d_scale(x, -1.0L));

    name = "cosh(x) == (e^x + e^-x) / 2";
    compare(name, cos, d_scale(d_add(exp_x, neg_exp_x), 0.5L));

    name = "sinh(x) == (e^x - e^-x) / 2";
    compare(name, sin, d_scale(d_sub(exp_x, neg_exp_x), 0.5L));

    if (debug != 0) printf("\n");
    sin = d_sin(x);
    cos = d_cos(x);
    tan = d_tan(x);
    sqr_sin_x = d_sqr(sin);
    sqr_cos_x = d_sqr(cos);
    sin_2x = d_sin(xpx);
    cos_2x = d_cos(xpx);

    name = "cos^2(x) + sin^2(x) == 1";
    compare(name, d_add(sqr_cos_x, sqr_sin_x), c1);

    name = "tan(x) == sin(x) / cos(x)";
    x_lt_pi_2 ? compare(name, tan, d_div(sin, cos)) : skip(name);

    name = "sin(2x) == 2 * sin(x) * cos(x)";
    compare(name, sin_2x, d_scale(d_mul(sin, cos), 2.0L));

    name = "cos(2x) == cos^2(x) - sin^2(x)";
    compare(name, cos_2x, d_sub(sqr_cos_x, sqr_sin_x));

    if (debug != 0) printf("\n");
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
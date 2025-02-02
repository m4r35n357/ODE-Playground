/*
 * Dual numbers, validation checks
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dual.h"

static int dp, debug = 0, total = 0, passed = 0, skipped = 0;

static real tolerance;

static void skip (char* name) {
    total++;
    skipped++;
    if (debug) fprintf(stderr, "%s SKIP%s %s%s%s\n", YLW, NRM, GRY, name, NRM);
}

static void compare (char* name, dual a, dual b) {
    real delta_max = 0.0L;
    total++;
    real delta_val = fabsl(a.val - b.val);
    if (delta_val > delta_max) {
        delta_max = delta_val;
    }
    if (debug == 2) {
        if (delta_val > tolerance) {
            fprintf(stderr, "  %sval  %s% .*Le % .*Le  %s%.1Le%s\n", RED, NRM, dp, a.val, dp, b.val, RED, delta_val, NRM);
        } else {
            fprintf(stderr, "  val  %s% .*Le % .*Le%s  %.1Le\n", GRY, dp, a.val, dp, b.val, NRM, delta_val);
        }
    }
    real delta_dot = fabsl(a.dot - b.dot);
    if (delta_dot > delta_max) {
        delta_max = delta_dot;
    }
    if (debug == 2) {
        if (delta_dot > tolerance) {
            fprintf(stderr, "  %sdot  %s% .*Le % .*Le%s  %.1Le%s\n", RED, NRM, dp, a.dot, dp, b.dot, RED, delta_dot, NRM);
        } else {
            fprintf(stderr, "  dot  %s% .*Le % .*Le%s  %.1Le\n", GRY, dp, a.dot, dp, b.dot, NRM, delta_dot);
        }
    }
    bool failed = delta_max > tolerance;
    if (debug) {
        if (failed) {
            fprintf(stderr, "%s FAIL%s %s%s%s\n", RED, NRM, WHT, name, NRM);
        } else {
            fprintf(stderr, "%s PASS%s %s\n", GRN, NRM, name);
        }
    }
    if (debug == 2) fprintf(stderr, "\n");
    if (!failed) passed++;
}

int main (int argc, char **argv) {
    PRINT_ARGS(argc, argv);
    CHECK(argc == 4 || argc == 5);

    dp = (int)strtol(argv[1], NULL, BASE);
    dual x = (dual){.val = strtold(argv[2], NULL), .dot = 0.5L};
    tolerance = strtold(argv[3], NULL); CHECK(tolerance > 0.0L);
    if (argc == 5) {
        debug = (int)strtol(argv[4], NULL, BASE); CHECK(debug == 0 || debug == 1 || debug == 2);
    }

    fprintf(stderr, "Dual Numbers %sx = %s%.1Lf%s\n", GRY, WHT, x.val, NRM);
    bool positive = x.val > 0.0L, non_zero = x.val != 0.0L, lt_pi_2 = fabsl(x.val) < 0.5L * acosl(-1.0L);
    dual D1 = d_dual(1.0L), xpx = d_scale(x, 2.0L);
    dual abs_x, inv_x, sqrt_x, ln_x;
    if (non_zero) abs_x = d_abs(x);
    if (positive || non_zero) inv_x = d_rec(x);
    if (positive) sqrt_x = d_sqrt(x);
    if (positive) ln_x = d_ln(x);
    dual sin_x = d_sin(x), cos_x = d_cos(x), tan_x = d_tan(x);
    dual sin2_x = d_sqr(sin_x), cos2_x = d_sqr(cos_x), sin_2x = d_sin(xpx), cos_2x = d_cos(xpx);
    dual sinh_x = d_sinh(x), cosh_x = d_cosh(x), tanh_x = d_tanh(x);
    dual sinh2_x = d_sqr(sinh_x), cosh2_x = d_sqr(cosh_x), sinh_2x = d_sinh(xpx), cosh_2x = d_cosh(xpx);
    dual sqr_x = d_sqr(x), exp_x = d_exp(x), neg_exp_x = d_exp(d_scale(x, -1.0L));
    dual gd_1 = d_ln(d_abs(d_div(d_add(sin_x, D1), cos_x)));

    char* name = "x * x == sqr(x)"; compare(name, d_mul(x, x), sqr_x);
    name = "sqr(x) / x == x"; non_zero ? compare(name, d_div(sqr_x, x), x) : skip(name);
    name = "x * 1 / x == 1"; non_zero ? compare(name, d_mul(x, inv_x), D1) : skip(name);
    name = "sqrt(x) * sqrt(x) == x"; positive ? compare(name, d_mul(sqrt_x, sqrt_x), x) : skip(name);
    name = "x / sqrt(x) == sqrt(x)"; positive ? compare(name, d_div(x, sqrt_x), sqrt_x) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "x^2.0 == sqr(x)"; positive ? compare(name, d_pow(x, 2.0L), sqr_x) : skip(name);
    name = "x^1.0 == x"; positive ? compare(name, d_pow(x, 1.0L), x) : skip(name);
    name = "x^0.5 == sqrt(x)"; positive ? compare(name, d_pow(x, 0.5L), sqrt_x): skip(name);
    name = "x^0.0 == 1"; positive ? compare(name, d_pow(x, 0.0L), D1) : skip(name);
    name = "x^-0.5 == 1 / sqrt(x)"; positive ? compare(name, d_pow(x, -0.5L), d_rec(sqrt_x)) : skip(name);
    name = "x^-1.0 == 1 / x"; positive ? compare(name, d_pow(x, -1.0L), inv_x) : skip(name);
    name = "x^-2.0 == 1 / sqr(x)"; positive ? compare(name, d_pow(x, -2.0L), d_rec(sqr_x)) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "sqr(x) * x^-3 == 1 / x"; positive ? compare(name, d_mul(sqr_x, d_pow(x, -3.0L)), inv_x) : skip(name);
    name = "sqr(x)^0.5 == |x|"; non_zero ? compare(name, d_pow(sqr_x, 0.5L), abs_x) : skip(name);
    name = "sqrt(sqr(x) == |x|"; non_zero ? compare(name, d_sqrt(sqr_x), abs_x) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "ln(e^x) == x"; compare(name, d_ln(d_exp(x)), x);
    name = "ln(sqr(x)) == ln(x) * 2"; positive ? compare(name, d_ln(sqr_x), d_scale(ln_x, 2.0L)) : skip(name);
    name = "ln(sqrt(x)) == ln(x) / 2"; positive ? compare(name, d_ln(sqrt_x), d_scale(ln_x, 0.5L)) : skip(name);
    name = "ln(1 / x) == - ln(x)"; positive ? compare(name, d_ln(inv_x), d_scale(ln_x, -1.0L)) : skip(name);
    name = "ln(x^-3) == -3*ln(x)"; positive ? compare(name, d_ln(d_pow(x, -3.0L)), d_scale(ln_x, -3.0L)) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "cosh^2(x) == 1 + sinh^2(x)"; compare(name, cosh2_x, d_add(D1, sinh2_x));
    name = "tanh(x) == sinh(x) / cosh(x)"; compare(name, tanh_x, d_div(sinh_x, cosh_x));
    name = "sinh(2x) == 2 * sinh(x) * cosh(x)"; compare(name, sinh_2x, d_scale(d_mul(sinh_x, cosh_x), 2.0L));
    name = "cosh(2x) == cosh^2(x) + sinh^2(x)"; compare(name, cosh_2x, d_add(cosh2_x, sinh2_x));

    if (debug) fprintf(stderr, "\n");

    name = "cosh(x) == (e^x + e^-x) / 2"; compare(name, cosh_x, d_scale(d_add(exp_x, neg_exp_x), 0.5L));
    name = "sinh(x) == (e^x - e^-x) / 2"; compare(name, sinh_x, d_scale(d_sub(exp_x, neg_exp_x), 0.5L));

    if (debug) fprintf(stderr, "\n");

    name = "arsinh(sinh(x)) == x"; compare(name, d_asinh(sinh_x), x);
    name = "arcosh(cosh(x)) == |x|"; non_zero ? compare(name, d_acosh(cosh_x), abs_x) : skip(name);
    name = "artanh(tanh(x)) == x"; compare(name, d_atanh(tanh_x), x);

    if (debug) fprintf(stderr, "\n");

    name = "cos^2(x) == 1 - sin^2(x)"; compare(name, cos2_x, d_sub(D1, sin2_x));
    name = "tan(x) == sin(x) / cos(x)"; lt_pi_2 ? compare(name, tan_x, d_div(sin_x, cos_x)) : skip(name);
    name = "sin(2x) == 2 * sin(x) * cos(x)"; compare(name, sin_2x, d_scale(d_mul(sin_x, cos_x), 2.0L));
    name = "cos(2x) == cos^2(x) - sin^2(x)"; compare(name, cos_2x, d_sub(cos2_x, sin2_x));

    if (debug) fprintf(stderr, "\n");

    name = "arcsin(sin(x)) == x"; if (lt_pi_2) {compare(name, d_asin(sin_x), x);} else skip(name);
    name = "arccos(cos(x)) == |x|"; non_zero ? compare(name, d_acos(cos_x), abs_x) : skip(name);
    name = "arctan(tan(x)) == x"; if (lt_pi_2) {compare(name, d_atan(tan_x), x);} else skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "arsinh(tan(x)) == gd^-1 x"; if (lt_pi_2) {compare(name, gd_1, d_asinh(tan_x));} else skip(name);
    name = "artanh(sin(x)) == gd^-1 x"; compare(name, gd_1, d_atanh(sin_x));
    name = "arcsin(tanh(gd^-1 x)) == x"; if (lt_pi_2) {compare(name, d_asin(d_tanh(gd_1)), x);} else skip(name);
    name = "arctan(sinh(gd^-1 x)) == x"; if (lt_pi_2) {compare(name, d_atan(d_sinh(gd_1)), x);} else skip(name);

    if (debug) fprintf(stderr, "\n");

    fprintf(stderr, "%sTotal%s %d  %sPASSED%s %d", WHT, NRM, total, GRN, NRM, passed);
    if (skipped) fprintf(stderr, "  %sSKIPPED%s %d", YLW, NRM, skipped);
    if (passed < total - skipped) {
        fprintf(stderr, "  %sFAILED%s %d\n", RED, NRM, total - passed - skipped);
        return 1;
    }
    fprintf(stderr, "\n");
    return 0;
}

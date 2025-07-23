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
    dual u = (dual){.val = strtold(argv[2], NULL), .dot = 0.5L};
    tolerance = strtold(argv[3], NULL); CHECK(tolerance > 0.0L);
    if (argc == 5) {
        debug = (int)strtol(argv[4], NULL, BASE); CHECK(debug == 0 || debug == 1 || debug == 2);
    }

    fprintf(stderr, "Dual Numbers %su = %s%.1Lf%s\n", GRY, WHT, u.val, NRM);
    bool positive = u.val > 0.0L, non_zero = u.val != 0.0L, lt_pi_2 = fabsl(u.val) < 0.5L * acosl(-1.0L);
    dual D1 = d_dual(1.0L), upu = d_scale(u, 2.0L);
    dual abs_u, inv_u, sqrt_u, ln_u;
    if (non_zero) abs_u = d_abs(u);
    if (positive || non_zero) inv_u = d_rec(u);
    if (positive) sqrt_u = d_sqrt(u);
    if (positive) ln_u = d_ln(u);
    dual sin_u = d_sin(u), cos_u = d_cos(u), tan_u = d_tan(u);
    dual sin2_u = d_sqr(sin_u), cos2_u = d_sqr(cos_u), sin_2u = d_sin(upu), cos_2u = d_cos(upu);
    dual sinh_u = d_sinh(u), cosh_u = d_cosh(u), tanh_u = d_tanh(u);
    dual sinh2_u = d_sqr(sinh_u), cosh2_u = d_sqr(cosh_u), sinh_2u = d_sinh(upu), cosh_2u = d_cosh(upu);
    dual sqr_u = d_sqr(u), exp_u = d_exp(u), neg_exp_u = d_exp(d_scale(u, -1.0L));
    dual gd_1 = d_ln(d_abs(d_div(d_add(sin_u, D1), cos_u)));

    char* name = "u * u == sqr(u)"; compare(name, d_mul(u, u), sqr_u);
    name = "sqr(u) / u == u"; non_zero ? compare(name, d_div(sqr_u, u), u) : skip(name);
    name = "u * 1 / u == 1"; non_zero ? compare(name, d_mul(u, inv_u), D1) : skip(name);
    name = "sqrt(u) * sqrt(u) == u"; positive ? compare(name, d_mul(sqrt_u, sqrt_u), u) : skip(name);
    name = "u / sqrt(u) == sqrt(u)"; positive ? compare(name, d_div(u, sqrt_u), sqrt_u) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "u^2.0 == sqr(u)"; positive ? compare(name, d_pow(u, 2.0L), sqr_u) : skip(name);
    name = "u^1.0 == u"; positive ? compare(name, d_pow(u, 1.0L), u) : skip(name);
    name = "u^0.5 == sqrt(u)"; positive ? compare(name, d_pow(u, 0.5L), sqrt_u): skip(name);
    name = "u^0.0 == 1"; positive ? compare(name, d_pow(u, 0.0L), D1) : skip(name);
    name = "u^-0.5 == 1 / sqrt(u)"; positive ? compare(name, d_pow(u, -0.5L), d_rec(sqrt_u)) : skip(name);
    name = "u^-1.0 == 1 / u"; positive ? compare(name, d_pow(u, -1.0L), inv_u) : skip(name);
    name = "u^-2.0 == 1 / sqr(u)"; positive ? compare(name, d_pow(u, -2.0L), d_rec(sqr_u)) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "sqr(u) * u^-3 == 1 / u"; positive ? compare(name, d_mul(sqr_u, d_pow(u, -3.0L)), inv_u) : skip(name);
    name = "sqr(u)^0.5 == |u|"; non_zero ? compare(name, d_pow(sqr_u, 0.5L), abs_u) : skip(name);
    name = "sqrt(sqr(u) == |u|"; non_zero ? compare(name, d_sqrt(sqr_u), abs_u) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "ln(e^u) == u"; compare(name, d_ln(d_exp(u)), u);
    name = "ln(sqr(u)) == ln(u) * 2"; positive ? compare(name, d_ln(sqr_u), d_scale(ln_u, 2.0L)) : skip(name);
    name = "ln(sqrt(u)) == ln(u) / 2"; positive ? compare(name, d_ln(sqrt_u), d_scale(ln_u, 0.5L)) : skip(name);
    name = "ln(1 / u) == - ln(u)"; positive ? compare(name, d_ln(inv_u), d_scale(ln_u, -1.0L)) : skip(name);
    name = "ln(u^-3) == -3*ln(u)"; positive ? compare(name, d_ln(d_pow(u, -3.0L)), d_scale(ln_u, -3.0L)) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "cosh^2(u) == 1 + sinh^2(u)"; compare(name, cosh2_u, d_add(D1, sinh2_u));
    name = "tanh(u) == sinh(u) / cosh(u)"; compare(name, tanh_u, d_div(sinh_u, cosh_u));
    name = "sinh(2u) == 2 * sinh(u) * cosh(u)"; compare(name, sinh_2u, d_scale(d_mul(sinh_u, cosh_u), 2.0L));
    name = "cosh(2u) == cosh^2(u) + sinh^2(u)"; compare(name, cosh_2u, d_add(cosh2_u, sinh2_u));

    if (debug) fprintf(stderr, "\n");

    name = "cosh(u) == (e^u + e^-u) / 2"; compare(name, cosh_u, d_scale(d_add(exp_u, neg_exp_u), 0.5L));
    name = "sinh(u) == (e^u - e^-u) / 2"; compare(name, sinh_u, d_scale(d_sub(exp_u, neg_exp_u), 0.5L));
    name = "tanh(u) == (e^u - e^-u) / (e^u + e^-u)"; compare(name, tanh_u, d_div(d_sub(exp_u, neg_exp_u), d_add(exp_u, neg_exp_u)));

    if (debug) fprintf(stderr, "\n");

    name = "arsinh(sinh(u)) == u"; compare(name, d_asinh(sinh_u), u);
    name = "arcosh(cosh(u)) == |u|"; non_zero ? compare(name, d_acosh(cosh_u), abs_u) : skip(name);
    name = "artanh(tanh(u)) == u"; compare(name, d_atanh(tanh_u), u);

    if (debug) fprintf(stderr, "\n");

    name = "cos^2(u) == 1 - sin^2(u)"; compare(name, cos2_u, d_sub(D1, sin2_u));
    name = "tan(u) == sin(u) / cos(u)"; lt_pi_2 ? compare(name, tan_u, d_div(sin_u, cos_u)) : skip(name);
    name = "sin(2u) == 2 * sin(u) * cos(u)"; compare(name, sin_2u, d_scale(d_mul(sin_u, cos_u), 2.0L));
    name = "cos(2u) == cos^2(u) - sin^2(u)"; compare(name, cos_2u, d_sub(cos2_u, sin2_u));

    if (debug) fprintf(stderr, "\n");

    name = "arcsin(sin(u)) == u"; if (lt_pi_2) {compare(name, d_asin(sin_u), u);} else skip(name);
    name = "arccos(cos(u)) == |u|"; non_zero ? compare(name, d_acos(cos_u), abs_u) : skip(name);
    name = "arctan(tan(u)) == u"; if (lt_pi_2) {compare(name, d_atan(tan_u), u);} else skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "arsinh(tan(u)) == gd^-1 u"; if (lt_pi_2) {compare(name, gd_1, d_asinh(tan_u));} else skip(name);
    name = "artanh(sin(u)) == gd^-1 u"; compare(name, gd_1, d_atanh(sin_u));
    name = "arcsin(tanh(gd^-1 u)) == u"; if (lt_pi_2) {compare(name, d_asin(d_tanh(gd_1)), u);} else skip(name);
    name = "arctan(sinh(gd^-1 u)) == u"; if (lt_pi_2) {compare(name, d_atan(d_sinh(gd_1)), u);} else skip(name);

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

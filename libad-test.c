/*
 * Tests for core Taylor Series Method, Horner's Method, and recurrence relations
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "taylor-ode.h"
#include "dual.h"

#define NRM "\x1B[0;37m"
#define WHT "\x1B[1;37m"
#define GRY "\x1B[1;30m"
#define GRN "\x1B[1;32m"
#define YLW "\x1B[1;33m"
#define RED "\x1B[1;31m"
#define CYN "\x1B[0;36m"
#define MGT "\x1B[0;35m"

static int k_max = 0, dp, n, debug = 0, total = 0, passed = 0, skipped = 0;

static real delta, delta_val, delta_dot, delta_max = 0.0L, tolerance;

static char *name_max = "N/A", *field_max = "N/A";

struct Parameters { real a, b, c; };

static series ad_scale (series s, series u, real a) { for (int k = 0; k < n; k++) s[k] = u[k] * a; return s; }

static series ad_add (series p, series u, series v) { for (int k = 0; k < n; k++) p[k] = u[k] + v[k]; return p; }

static series ad_sub (series m, series u, series v) { for (int k = 0; k < n; k++) m[k] = u[k] - v[k]; return m; }

static series ad_abs (series a, series u) { for (int k = 0; k < n; k++) a[k] = t_abs(u, k); return a; }

static series ad_mul (series p, series u, series v) { for (int k = 0; k < n; k++) p[k] = t_mul(u, v, k); return p; }

static series ad_div (series q, series u, series v) { for (int k = 0; k < n; k++) t_div(q, u, v, k); return q; }

static series ad_sqr (series s, series u) { for (int k = 0; k < n; k++) s[k] = t_sqr(u, k); return s; }

static series ad_sqrt (series r, series u) { for (int k = 0; k < n; k++) t_sqrt(r, u, k); return r; }

static series ad_exp (series e, series u) { for (int k = 0; k < n; k++) t_exp(e, u, k); return e; }

static void ad_sin_cos (series s, series c, series u, bool trig) { for (int k = 0; k < n; k++) t_sin_cos(s, c, u, k, trig); }

static void ad_tan_sec2 (series t, series s2, series u, bool trig) { for (int k = 0; k < n; k++) t_tan_sec2(t, s2, u, k, trig); }

static series ad_pwr (series p, series u, real a) { for (int k = 0; k < n; k++) t_pwr(p, u, a, k); return p; }

static series ad_ln (series l, series u) { for (int k = 0; k < n; k++) t_ln(l, u, k); return l; }

static void ad_asin (series as, series du_df, series u, bool trig) { for (int k = 0; k < n; k++) t_asin(as, du_df, u, k, trig); }

static void ad_acos (series ac, series du_df, series u, bool trig) { for (int k = 0; k < n; k++) t_acos(ac, du_df, u, k, trig); }

static void ad_atan (series at, series du_df, series u, bool trig) { for (int k = 0; k < n; k++) t_atan(at, du_df, u, k, trig); }

triplet ode (series x, series y, series z, parameters *p, int k) {
    return (triplet) {
        .x = p->a * x[k],
        .y = p->b * y[k],
        .z = p->c * z[k]
    };
}

static void skip (char* name) {
    total++;
    skipped++;
    if (debug) fprintf(stderr, "%s SKIP%s %s%s%s\n", YLW, NRM, GRY, name, NRM);
}

static void compare (char* name, series a, series b) {
    total++;
    for (int k = 0; k < n; k++) {
        delta = fabsl(a[k] - b[k]);
        if (delta > delta_max) {
            delta_max = delta;
            name_max = name;
            k_max = k;
        }
        if (!isfinite(delta) || delta > tolerance) {
            fprintf(stderr, "%s FAIL%s %s%s%s\n  k=%d  LHS: %+.*Le  RHS: %+.*Le  (%.1Le)\n",
                    RED, NRM, WHT, name, NRM, k, dp, a[k], dp, b[k], delta);
            return;
        }
        if (debug >= 2) {
            if (!k) fprintf(stderr, "\n");
            fprintf(stderr, "%s  DEBUG%s  k: %2d  %+.*Le %+.*Le  (%.1Le)\n",
                    NRM, NRM, k, dp, a[k], dp, b[k], delta);
        }
    }
    if (debug) fprintf(stderr, "%s PASS%s %s%s%s\n", GRN, NRM, WHT, name, NRM);
    passed++;
}

static void compare_s_d (char* name, series a, dual b) {
    total++;
    delta_val = fabsl(a[0] - b.val);
    if (delta_val > delta_max) {
        delta_max = delta_val;
        name_max = name;
        field_max = "VAL";
    }
    if (fabsl(delta_val) > tolerance) {
        fprintf(stderr, "%s FAIL%s %s%s%s\n  [0]VAL  LHS: %+.*Le  RHS: %+.*Le  (%+.3Le)\n", RED, NRM, WHT, name, NRM, dp, a[0], dp, b.val, delta_val);
        return;
    }
    delta_dot = fabsl(a[1] - b.dot);
    if (delta_dot > delta_max) {
        delta_max = delta_dot;
        name_max = name;
        field_max = "DOT";
    }
    if (fabsl(delta_dot) > tolerance) {
        fprintf(stderr, "%s FAIL%s %s%s%s\n  [1]DOT  LHS: %+.*Le  RHS: %+.*Le  (%+.3Le)\n", RED, NRM, WHT, name, NRM, dp, a[1], dp, b.dot, delta_dot);
        return;
    }
    if (debug >= 2) fprintf(stderr, "\n");
    if (debug >= 2) fprintf(stderr, "%s  DEBUG%s  %+.*Le %+.*Le  diff %+.3Le\n", NRM, NRM, dp, a[0], dp, b.val, delta_val);
    if (debug >= 2) fprintf(stderr, "%s  DEBUG%s  %+.*Le %+.*Le  diff %+.3Le\n", NRM, NRM, dp, a[1], dp, b.dot, delta_dot);
    if (debug) fprintf(stderr, "%s PASS%s %s%s%s\n", GRN, NRM, WHT, name, NRM);
    passed++;
}

int main (int argc, char **argv) {
    PRINT_ARGS(argc, argv);
    CHECK(argc == 5 || argc == 6);

    dp = (int)strtol(argv[1], NULL, BASE);
    n = (int)strtol(argv[2], NULL, BASE); CHECK(n > 8);
    series x = tsm_jet(n + 1);
    for (int k = 0; k <= n; k++) {
        x[k] = !k ? strtold(argv[3], NULL) : 0.5L / (k * k);
    }
    series xs = tsm_jet(2);
    xs[0] = x[0];
    xs[1] = x[1];
    dual xd = (dual){.val = x[0], .dot = 0.5L};
    dual D1 = d_dual(1.0L);
    tolerance = strtold(argv[4], NULL); CHECK(tolerance > 0.0L);
    if (argc == 6) {
        debug = (int)strtol(argv[5], NULL, BASE); CHECK(debug == 0 || debug == 1 || debug == 2);
    }

    fprintf(stderr, "Horner Summation ");
    series s = tsm_jet(n);
    s[0] = 1.0L; s[1] = 3.0L; s[2] = 0.0L; s[3] = 2.0L;
    CHECK(horner(s, 3, 2.0L) == 23.0L); fprintf(stderr, ".");
    s[0] = 3.0L; s[1] = -1.0L; s[2] = 2.0L; s[3] = -4.0L; s[4] = 0.0L; s[5] = 1.0L;
    CHECK(horner(s, 5, 3.0L) == 153.0L); fprintf(stderr, ".");
    s[0] = 1.0L; s[1] = -4.0L; s[2] = 0.0L; s[3] = 0.0L; s[4] = 2.0L; s[5] = 3.0L; s[6] = 0.0L; s[7] = -2.0L;
    CHECK(horner(s, 7, -2.0L) == 201.0L); fprintf(stderr, ". %sOK%s", GRN, NRM);

    fprintf(stderr, ", Taylor Series Method ");
    controls c = {.order=n, .step=0, .steps=10, .step_size=0.1L};
    parameters p = {.a=1.0L, .b=0.0L, .c=-1.0L};
    series3 *j = malloc(sizeof (series3)); CHECK(j);
    j->x = tsm_const(n + 1, 1.0L);
    j->y = tsm_const(n + 1, 1.0L);
    j->z = tsm_const(n + 1, 1.0L);
    while (tsm_gen(&c, j, &p)) fprintf(stderr, ".");
    CHECK(fabsl(j->x[0] - expl(p.a)) < tolerance);
    CHECK(fabsl(j->y[0] - expl(p.b)) < tolerance);
    CHECK(fabsl(j->z[0] - expl(p.c)) < tolerance);
    fprintf(stderr, " %sOK%s\n", GRN, NRM);

    fprintf(stderr, "Recurrence Relations %sx = %s%s%.1Lf%s\n", MGT, NRM, WHT, x[0], NRM);
    bool positive = x[0] > 0.0L, non_zero = x[0] != 0.0L, lt_1 = fabsl(x[0]) < 1.0L, gt_1 = x[0] > 1.0L, lt_pi_2 = fabsl(x[0]) < 0.5L * acosl(-1.0L);
    series r1 = tsm_jet(n), r2 = tsm_jet(n), r3 = tsm_jet(n), S1 = tsm_const(n, 1.0L);
    series abs_x = tsm_jet(n), inv_x = tsm_jet(n), sqrt_x = tsm_jet(n), ln_x = tsm_jet(n);
    if (non_zero) ad_abs(abs_x, x);
    if (non_zero) ad_div(inv_x, S1, x);
    if (positive) ad_sqrt(sqrt_x, x);
    if (positive) ad_ln(ln_x, x);
    series sin_x = tsm_jet(n), cos_x = tsm_jet(n), tan_x = tsm_jet(n), sec2_x = tsm_jet(n);
    series sin2_x = tsm_jet(n), cos2_x = tsm_jet(n), tan2_x = tsm_jet(n), sin_2x = tsm_jet(n), cos_2x = tsm_jet(n);
    ad_sin_cos(sin_x, cos_x, x, true);
    ad_tan_sec2(tan_x, sec2_x, x, true);
    ad_sqr(sin2_x, sin_x);
    ad_sqr(cos2_x, cos_x);
    ad_sin_cos(sin_2x, cos_2x, ad_scale(r1, x, 2.0L), true);
    series sinh_x = tsm_jet(n), cosh_x = tsm_jet(n), tanh_x = tsm_jet(n), sech2_x = tsm_jet(n);
    series sinh2_x = tsm_jet(n), cosh2_x = tsm_jet(n), tanh2_x = tsm_jet(n), sinh_2x = tsm_jet(n), cosh_2x = tsm_jet(n);
    ad_sin_cos(sinh_x, cosh_x, x, false);
    ad_tan_sec2(tanh_x, sech2_x, x, false);
    ad_sqr(sinh2_x, sinh_x);
    ad_sqr(cosh2_x, cosh_x);
    ad_sin_cos(sinh_2x, cosh_2x, ad_scale(r1, x, 2.0L), false);
    series sqr_x = tsm_jet(n), exp_x = tsm_jet(n), neg_exp_x = tsm_jet(n), gd_1 = tsm_jet(n);
    ad_sqr(sqr_x, x);
    ad_exp(exp_x, x);
    ad_exp(neg_exp_x, ad_scale(r1, x, -1.0L));
    ad_ln(gd_1, ad_abs(r3, ad_div(r2, ad_add(r1, sin_x, S1), cos_x)));

    char* name = "inv(x)"; non_zero ? compare_s_d(name, ad_div(r1, S1, xs), d_div(D1, xd)) : skip(name);
    name = "sqr(x)"; compare_s_d(name, ad_sqr(r1, xs), d_sqr(xd));
    name = "sqrt(x)"; positive ? compare_s_d(name, ad_sqrt(r1, xs), d_sqrt(xd)) : skip(name);

    name = "exp(x)"; compare_s_d(name, ad_exp(r1, xs), d_exp(xd));

    ad_sin_cos(r1, r2, xs, true);
    name = "sin(x)"; compare_s_d(name, r1, d_sin(xd));
    name = "cos(x)"; compare_s_d(name, r2, d_cos(xd));
    name = "tan(x)"; ad_tan_sec2(r1, r2, xs, true); compare_s_d(name, r1, d_tan(xd));

    ad_sin_cos(r1, r2, xs, false);
    name = "sinh(x)"; compare_s_d(name, r1, d_sinh(xd));
    name = "cosh(x)"; compare_s_d(name, r2, d_cosh(xd));
    name = "tanh(x)"; ad_tan_sec2(r1, r2, xs, false); compare_s_d(name, r1, d_tanh(xd));

    name = "ln(x)"; positive ? compare_s_d(name, ad_ln(r1, xs), d_ln(xd)) : skip(name);

    name = "asin(x)"; if (lt_1) {ad_asin(r1, r2, xs, true); compare_s_d(name, r1, d_asin(xd));} else skip(name);
    name = "acos(x)"; if (lt_1) {ad_acos(r1, r2, xs, true); compare_s_d(name, r1, d_acos(xd));} else skip(name);
    name = "atan(x)"; ad_atan(r1, r2, xs, true); compare_s_d(name, r1, d_atan(xd));

    ad_asin(r1, r2, xs, false);
    name = "asinh(x)"; compare_s_d(name, r1, d_asinh(xd));
    name = "acosh(x)"; if (gt_1) {ad_acos(r1, r2, xs, false); compare_s_d(name, r1, d_acosh(xd));} else skip(name);
    name = "atanh(x)"; if (lt_1) {ad_atan(r1, r2, xs, false); compare_s_d(name, r1, d_atanh(xd));}  else skip(name);

    name = "x^1.5"; positive ? compare_s_d(name, ad_pwr(r1, xs, 1.5L), d_pow(xd, 1.5L)) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "x * x == sqr(x)"; compare(name, ad_mul(r1, x, x), sqr_x);
    name = "sqr(x) / x == x"; non_zero ? compare(name, ad_div(r1, sqr_x, x), x) : skip(name);
    name = "x * (1 / x) == 1"; non_zero ? compare(name, ad_mul(r1, x, inv_x), S1) : skip(name);
    name = "sqrt(x) * sqrt(x) == x"; positive ? compare(name, ad_mul(r1, sqrt_x, sqrt_x), x) : skip(name);
    name = "x / sqrt(x) == sqrt(x)"; positive ? compare(name, ad_div(r1, x, sqrt_x), sqrt_x) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "x^2.0 == sqr(x)"; positive ? compare(name, ad_pwr(r1, x, 2.0L), sqr_x) : skip(name);
    name = "x^1.0 == x"; positive ? compare(name, ad_pwr(r1, x, 1.0L), x) : skip(name);
    name = "x^0.5 == sqrt(x)"; positive ? compare(name, ad_pwr(r1, x, 0.5L), sqrt_x): skip(name);
    name = "x^0.0 == 1"; positive ? compare(name, ad_pwr(r1, x, 0.0L), S1) : skip(name);
    name = "x^-0.5 == 1 / sqrt(x)"; positive ? compare(name, ad_pwr(r1, x, -0.5L), ad_div(r2, S1, sqrt_x)) : skip(name);
    name = "x^-1.0 == 1 / x"; positive ? compare(name, ad_pwr(r1, x, -1.0L), inv_x) : skip(name);
    name = "x^-2.0 == 1 / sqr(x)"; positive ? compare(name, ad_pwr(r1, x, -2.0L), ad_div(r2, S1, sqr_x)) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "sqr(x) * x^-3 == 1 / x"; positive ? compare(name, ad_mul(r1, sqr_x, ad_pwr(r2, x, -3.0L)), inv_x) : skip(name);
    name = "sqr(x)^0.5 == |x|"; non_zero ? compare(name, ad_pwr(r1, sqr_x, 0.5L), abs_x) : skip(name);
    name = "sqrt(sqr(x) == |x|"; non_zero ? compare(name, ad_sqrt(r1, sqr_x), abs_x) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "ln(e^x) == x"; compare(name, ad_ln(r1, exp_x), x);
    name = "ln(sqr(x)) == ln(x) * 2"; positive ? compare(name, ad_ln(r1, sqr_x), ad_scale(r2, ln_x, 2.0L)) : skip(name);
    name = "ln(sqrt(x)) == ln(x) / 2"; positive ? compare(name, ad_ln(r1, sqrt_x), ad_scale(r2, ln_x, 0.5L)) : skip(name);
    name = "ln(1 / x) == -ln(x)"; positive ? compare(name, ad_ln(r1, inv_x), ad_scale(r2, ln_x, -1.0L)) : skip(name);
    name = "ln(x^-3) == -3ln(x)"; positive ? compare(name, ad_ln(r1, ad_pwr(r2, x, -3.0L)), ad_scale(r3, ln_x, -3.0L)) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "cosh^2(x) - sinh^2(x) == 1"; compare(name, ad_sub(r1, cosh2_x, sinh2_x), S1);
    name = "sech^2(x) + tanh^2(x) == 1"; compare(name, ad_add(r1, sech2_x, ad_sqr(tanh2_x, tanh_x)), S1);
    name = "tanh(x) == sinh(x) / cosh(x)"; compare(name, tanh_x, ad_div(r1, sinh_x, cosh_x));
    name = "sech^2(x) == 1 / cosh^2(x)"; compare(name, sech2_x, ad_div(r1, S1, cosh2_x));
    name = "sinh(2x) == 2 * sinh(x) * cosh(x)"; compare(name, sinh_2x, ad_scale(r1, ad_mul(r2, sinh_x, cosh_x), 2.0L));
    name = "cosh(2x) == cosh^2(x) + sinh^2(x)"; compare(name, cosh_2x, ad_add(r1, cosh2_x, sinh2_x));

    if (debug) fprintf(stderr, "\n");

    name = "cosh(x) == (e^x + e^-x) / 2"; compare(name, cosh_x, ad_scale(r1, ad_add(r2, exp_x, neg_exp_x), 0.5L));
    name = "sinh(x) == (e^x - e^-x) / 2"; compare(name, sinh_x, ad_scale(r1, ad_sub(r2, exp_x, neg_exp_x), 0.5L));

    if (debug) fprintf(stderr, "\n");

    name = "arsinh(sinh(x)) == x"; ad_asin(r1, r2, sinh_x, false); compare(name, r1, x);
    name = "arcosh(cosh(x)) == |x|"; if (non_zero) {ad_acos(r1, r2, cosh_x, false); compare(name, r1, abs_x);} else skip(name);
    name = "artanh(tanh(x)) == x"; ad_atan(r1, r2, tanh_x, false); compare(name, r1, x);

    if (debug) fprintf(stderr, "\n");

    name = "cos^2(x) + sin^2(x) == 1"; compare(name, ad_add(r1, cos2_x, sin2_x), S1);
    name = "sec^2(x) - tan^2(x) == 1"; lt_pi_2 ? compare(name, ad_sub(r1, sec2_x, ad_sqr(tan2_x, tan_x)), S1) : skip(name);
    name = "tan(x) == sin(x) / cos(x)"; lt_pi_2 ? compare(name, tan_x, ad_div(r1, sin_x, cos_x)) : skip(name);
    name = "sec^2(x) == 1 / cos^2(x)"; lt_pi_2 ? compare(name, sec2_x, ad_div(r1, S1, cos2_x)) : skip(name);
    name = "sin(2x) == 2 * sin(x) * cos(x)"; compare(name, sin_2x, ad_scale(r1, ad_mul(r2, sin_x, cos_x), 2.0L));
    name = "cos(2x) == cos^2(x) - sin^2(x)"; compare(name, cos_2x, ad_sub(r1, cos2_x, sin2_x));

    if (debug) fprintf(stderr, "\n");

    name = "arcsin(sin(x)) == x"; if (lt_pi_2) {ad_asin(r1, r2, sin_x, true); compare(name, r1, x);} else skip(name);
    name = "arccos(cos(x)) == |x|"; if (non_zero) {ad_acos(r1, r2, cos_x, true); compare(name, r1, abs_x);} else skip(name);
    name = "arctan(tan(x)) == x"; if (lt_pi_2) {ad_atan(r1, r2, tan_x, true); compare(name, r1, x);} else skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "arsinh(tan(x)) == gd^-1 x"; if (lt_pi_2) {ad_asin(r1, r2, tan_x, false); compare(name, gd_1, r1);} else skip(name);
    name = "artanh(sin(x)) == gd^-1 x"; ad_atan(r1, r2, sin_x, false); compare(name, gd_1, r1);
    name = "arcsin(tanh(gd^-1 x)) == x"; if (lt_pi_2) {ad_tan_sec2(r3, r2, gd_1, false); ad_asin(r1, r2, r3, true); compare(name, r1, x);} else skip(name);
    name = "arctan(sinh(gd^-1 x)) == x"; if (lt_pi_2) {ad_sin_cos(r3, r2, gd_1, false); ad_atan(r1, r2, r3, true); compare(name, r1, x);} else skip(name);

    if (debug) fprintf(stderr, "\n");
    fprintf(stderr, "%sTotal%s: %d, %sPASSED%s %d", WHT, NRM, total, GRN, NRM, passed);
    if (skipped) fprintf(stderr, ", %sSKIPPED%s %d", YLW, NRM, skipped);
    if (passed == total - skipped) {
        fprintf(stderr, "\nDelta %s%.1Le%s %s%s%s k == %d\n", WHT, delta_max, NRM, MGT, name_max, NRM, k_max);
        return 0;
    } else {
        fprintf(stderr, ", %sFAILED%s %d\n\n", RED, NRM, total - passed - skipped);
        return 3;
    }
}

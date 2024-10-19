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

static int dp, n, debug = 0, total = 0, passed = 0, skipped = 0;

static real tolerance;

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

static series ad_ln (series u, series e) { for (int k = 0; k < n; k++) t_ln(u, e, k); return u; }

static void ad_asin_cos (series u, series c, series s, bool trig) { for (int k = 0; k < n; k++) t_asin_cos(u, c, s, k, trig); }

static void ad_acos_sin (series u, series s, series c, bool trig) { for (int k = 0; k < n; k++) t_acos_sin(u, s, c, k, trig); }

static void ad_atan_sec2 (series u, series s, series t, bool trig) { for (int k = 0; k < n; k++) t_atan_sec2(u, s, t, k, trig); }

static series ad_pwr (series p, series u, real a) { for (int k = 0; k < n; k++) t_pwr(p, u, a, k); return p; }

triplet ode (series x, series y, series z, const model *p, int k) {
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
    real delta_max = 0.0L;
    total++;
    for (int k = 0; k < n; k++) {
        real delta = fabsl(a[k] - b[k]);
        if (delta > delta_max) {
            delta_max = delta;
        }
        if (debug == 2) {
            if (!k) fprintf(stderr, "\n");
            if (delta > tolerance) {
                fprintf(stderr, "  %s%2d  %s% .*Le % .*Le%s  %.1Le%s\n", RED, k, NRM, dp, a[k], dp, b[k], RED, delta, NRM);
            } else {
                fprintf(stderr, "  %2d  %s% .*Le % .*Le%s  %.1Le\n", k, GRY, dp, a[k], dp, b[k], NRM, delta);
            }
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
    if (!failed) passed++;
}

static void compare_s_d (char* name, series a, dual b) {
    int n_tmp = n;
    n = 2;
    compare(name, a, (real []){b.val, b.dot});
    n = n_tmp;
}

int main (int argc, char **argv) {
    PRINT_ARGS(argc, argv);
    CHECK(argc == 5 || argc == 6);

    dp = (int)strtol(argv[1], NULL, BASE);
    n = (int)strtol(argv[2], NULL, BASE); CHECK(n > 8);
    series x = tsm_jet(n + 1);
    for (int k = 0, s = 1; k <= n; k++, s *= -1) {
        x[k] = !k ? strtold(argv[3], NULL) : 0.5L * s / SQR(k);
    }
    dual xd = {.val = x[0], .dot = x[1]};
    tolerance = strtold(argv[4], NULL); CHECK(tolerance > 0.0L);
    if (argc == 6) {
        debug = (int)strtol(argv[5], NULL, BASE); CHECK(debug == 0 || debug == 1 || debug == 2);
    }

    fprintf(stderr, "%sHorner Summation ", GRY);
    series s = tsm_jet(n);
    s[0] = 1.0L; s[1] = 3.0L; s[2] = 0.0L; s[3] = 2.0L;
    CHECK(horner(s, 3, 2.0L) == 23.0L); fprintf(stderr, ".");
    s[0] = 3.0L; s[1] = -1.0L; s[2] = 2.0L; s[3] = -4.0L; s[4] = 0.0L; s[5] = 1.0L;
    CHECK(horner(s, 5, 3.0L) == 153.0L); fprintf(stderr, ".");
    s[0] = 1.0L; s[1] = -4.0L; s[2] = 0.0L; s[3] = 0.0L; s[4] = 2.0L; s[5] = 3.0L; s[6] = 0.0L; s[7] = -2.0L;
    CHECK(horner(s, 7, -2.0L) == 201.0L); fprintf(stderr, ".%s OK%s", NRM, GRY);

    fprintf(stderr, ", Taylor Series Method ");
    controls c = {.order=n, .step=0, .steps=10, .h=0.1L};
    model p = {.a=1.0L, .b=0.0L, .c=-1.0L};
    xyz *_ = malloc(sizeof (xyz)); CHECK(_);
    _->x = tsm_jet(n + 1); _->x[0] = 1.0L;
    _->y = tsm_jet(n + 1); _->y[0] = 1.0L;
    _->z = tsm_jet(n + 1); _->z[0] = 1.0L;
    while (tsm_gen(&c, _, &p)) fprintf(stderr, ".");
    CHECK(fabsl(_->x[0] - expl(p.a)) < tolerance);
    CHECK(fabsl(_->y[0] - expl(p.b)) < tolerance);
    CHECK(fabsl(_->z[0] - expl(p.c)) < tolerance);
    fprintf(stderr, "%s OK\n", NRM);

    fprintf(stderr, "Taylor Arithmetic %sx = %s%.1Lf%s\n", GRY, WHT, x[0], NRM);
    bool positive = x[0] > 0.0L, non_zero = x[0] != 0.0L, lt_1 = fabsl(x[0]) < 1.0L, gt_1 = x[0] > 1.0L,
         lt_pi_2 = fabsl(x[0]) < 0.5L * acosl(-1.0L);
    series r1 = tsm_jet(n), r2 = tsm_jet(n), r3 = tsm_jet(n), S1 = tsm_jet(n); S1[0] = 1.0L;
    series abs_x = tsm_jet(n), inv_x = tsm_jet(n), sqrt_x = tsm_jet(n), ln_x = tsm_jet(n);
    if (non_zero) ad_abs(abs_x, x);
    if (non_zero) ad_div(inv_x, S1, x);
    if (positive) ad_sqrt(sqrt_x, x);
    if (positive) ad_ln(ln_x, x);
    series sin_x = tsm_jet(n), cos_x = tsm_jet(n), tan_x = tsm_jet(n), sec2_x = tsm_jet(n);
    series asin_x = tsm_jet(n), acos_x = tsm_jet(n), atan_x = tsm_jet(n);
    series sin2_x = tsm_jet(n), cos2_x = tsm_jet(n), tan2_x = tsm_jet(n), sin_2x = tsm_jet(n), cos_2x = tsm_jet(n);
    ad_sin_cos(sin_x, cos_x, x, true);
    if (lt_pi_2) ad_tan_sec2(tan_x, sec2_x, x, true);
    ad_sqr(sin2_x, sin_x);
    ad_sqr(cos2_x, cos_x);
    ad_sin_cos(sin_2x, cos_2x, ad_scale(r1, x, 2.0L), true);
    series sinh_x = tsm_jet(n), cosh_x = tsm_jet(n), tanh_x = tsm_jet(n), sech2_x = tsm_jet(n);
    series asinh_x = tsm_jet(n), acosh_x = tsm_jet(n), atanh_x = tsm_jet(n);
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

    char* name = "inv(x)"; non_zero ? compare_s_d(name, inv_x, d_inv(xd)) : skip(name);
    name = "sqr(x)"; compare_s_d(name, sqr_x, d_sqr(xd));
    name = "sqrt(x)"; positive ? compare_s_d(name, sqrt_x, d_sqrt(xd)) : skip(name);

    name = "exp(x)"; compare_s_d(name, exp_x, d_exp(xd));

    name = "sin(x)"; compare_s_d(name, sin_x, d_sin(xd));
    name = "cos(x)"; compare_s_d(name, cos_x, d_cos(xd));
    name = "tan(x)"; lt_pi_2 ? compare_s_d(name, tan_x, d_tan(xd)) : skip(name);

    name = "sinh(x)"; compare_s_d(name, sinh_x, d_sinh(xd));
    name = "cosh(x)"; compare_s_d(name, cosh_x, d_cosh(xd));
    name = "tanh(x)"; compare_s_d(name, tanh_x, d_tanh(xd));

    name = "ln(x)"; positive ? compare_s_d(name, ln_x, d_ln(xd)) : skip(name);

    name = "asin(x)"; if (lt_1) {ad_asin_cos(asin_x, r2, x, true); compare_s_d(name, asin_x, d_asin(xd));} else skip(name);
    name = "acos(x)"; if (lt_1) {ad_acos_sin(acos_x, r2, x, true); compare_s_d(name, acos_x, d_acos(xd));} else skip(name);
    name = "atan(x)"; ad_atan_sec2(atan_x, r2, x, true); compare_s_d(name, atan_x, d_atan(xd));

    name = "asinh(x)"; ad_asin_cos(asinh_x, r2, x, false); compare_s_d(name, asinh_x, d_asinh(xd));
    name = "acosh(x)"; if (gt_1) {ad_acos_sin(acosh_x, r2, x, false); compare_s_d(name, acosh_x, d_acosh(xd));} else skip(name);
    name = "atanh(x)"; if (lt_1) {ad_atan_sec2(atanh_x, r2, x, false); compare_s_d(name, atanh_x, d_atanh(xd));}  else skip(name);

    name = "x^1.5"; positive ? compare_s_d(name, ad_pwr(r1, x, 1.5L), d_pow(xd, 1.5L)) : skip(name);

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

    name = "cosh^2(x) == 1 + sinh^2(x)"; compare(name, cosh2_x, ad_add(r1, S1, sinh2_x));
    name = "sech^2(x) == 1 - tanh^2(x)"; compare(name, sech2_x, ad_sub(r1, S1, ad_sqr(tanh2_x, tanh_x)));
    name = "tanh(x) == sinh(x) / cosh(x)"; compare(name, tanh_x, ad_div(r1, sinh_x, cosh_x));
    name = "sech^2(x) == 1 / cosh^2(x)"; compare(name, sech2_x, ad_div(r1, S1, cosh2_x));
    name = "sinh(2x) == 2 * sinh(x) * cosh(x)"; compare(name, sinh_2x, ad_scale(r1, ad_mul(r2, sinh_x, cosh_x), 2.0L));
    name = "cosh(2x) == cosh^2(x) + sinh^2(x)"; compare(name, cosh_2x, ad_add(r1, cosh2_x, sinh2_x));

    if (debug) fprintf(stderr, "\n");

    name = "cosh(x) == (e^x + e^-x) / 2"; compare(name, cosh_x, ad_scale(r1, ad_add(r2, exp_x, neg_exp_x), 0.5L));
    name = "sinh(x) == (e^x - e^-x) / 2"; compare(name, sinh_x, ad_scale(r1, ad_sub(r2, exp_x, neg_exp_x), 0.5L));

    if (debug) fprintf(stderr, "\n");

    name = "arsinh(sinh(x)) == x"; ad_asin_cos(r1, r2, sinh_x, false); compare(name, r1, x);
    name = "arcosh(cosh(x)) == |x|"; if (non_zero) {ad_acos_sin(r1, r2, cosh_x, false); compare(name, r1, abs_x);} else skip(name);
    name = "artanh(tanh(x)) == x"; ad_atan_sec2(r1, r2, tanh_x, false); compare(name, r1, x);

    if (debug) fprintf(stderr, "\n");

    name = "cos^2(x) == 1 - sin^2(x)"; compare(name, cos2_x, ad_sub(r1, S1, sin2_x));
    name = "sec^2(x) == 1 + tan^2(x)"; lt_pi_2 ? compare(name, sec2_x, ad_add(r1, S1, ad_sqr(tan2_x, tan_x))) : skip(name);
    name = "tan(x) == sin(x) / cos(x)"; lt_pi_2 ? compare(name, tan_x, ad_div(r1, sin_x, cos_x)) : skip(name);
    name = "sec^2(x) == 1 / cos^2(x)"; lt_pi_2 ? compare(name, sec2_x, ad_div(r1, S1, cos2_x)) : skip(name);
    name = "sin(2x) == 2 * sin(x) * cos(x)"; compare(name, sin_2x, ad_scale(r1, ad_mul(r2, sin_x, cos_x), 2.0L));
    name = "cos(2x) == cos^2(x) - sin^2(x)"; compare(name, cos_2x, ad_sub(r1, cos2_x, sin2_x));

    if (debug) fprintf(stderr, "\n");

    name = "arcsin(sin(x)) == x"; if (lt_pi_2) {ad_asin_cos(r1, r2, sin_x, true); compare(name, r1, x);} else skip(name);
    name = "arccos(cos(x)) == |x|"; if (non_zero) {ad_acos_sin(r1, r2, cos_x, true); compare(name, r1, abs_x);} else skip(name);
    name = "arctan(tan(x)) == x"; if (lt_pi_2) {ad_atan_sec2(r1, r2, tan_x, true); compare(name, r1, x);} else skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "arsinh(tan(x)) == gd^-1 x"; if (lt_pi_2) {ad_asin_cos(r1, r2, tan_x, false); compare(name, gd_1, r1);} else skip(name);
    name = "artanh(sin(x)) == gd^-1 x"; ad_atan_sec2(r1, r2, sin_x, false); compare(name, gd_1, r1);
    name = "arcsin(tanh(gd^-1 x)) == x"; if (lt_pi_2) {ad_tan_sec2(r3, r2, gd_1, false); ad_asin_cos(r1, r2, r3, true); compare(name, r1, x);} else skip(name);
    name = "arctan(sinh(gd^-1 x)) == x"; if (lt_pi_2) {ad_sin_cos(r3, r2, gd_1, false); ad_atan_sec2(r1, r2, r3, true); compare(name, r1, x);} else skip(name);

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

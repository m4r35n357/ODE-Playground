/*
 * Tests for core Taylor Series Method, Horner's Method, and recurrence relations
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "taylor-ode.h"
#include "dual.h"

static int k_max = 0, dp, n, debug = 0, total = 0, passed = 0, skipped = 0;

static real delta, delta_max, tolerance;

real D0 = 0.0L;
real D01 = 0.1L;
real D05 = 0.5L;
real D_05 = -0.5L;
real D1 = 1.0L;
real D2 = 2.0L;
real D3 = 3.0L;
real D_1 = -1.0L;
real D_2 = -2.0L;
real D_3 = -3.0L;

static char *name_max = "N/A";

struct Parameters { real a, b, c; };

static series ad_scale (series s, series u, real a) { for (int k = 0; k < n; k++) s[k] = u[k] * a; return s; }

static series ad_add (series p, series u, series v) { for (int k = 0; k < n; k++) p[k] = u[k] + v[k]; return p; }

static series ad_sub (series m, series u, series v) { for (int k = 0; k < n; k++) m[k] = u[k] - v[k]; return m; }

static series ad_abs (series a, series u) { for (int k = 0; k < n; k++) a[k] = t_abs(u, k); return a; }

static series ad_mul (series p, series u, series v) { for (int k = 0; k < n; k++) p[k] = t_mul(u, v, k); return p; }

static series ad_div (series q, series u, series v) { for (int k = 0; k < n; k++) t_div(q, u, v, k); return q; }

static series ad_rec (series r, series v) { for (int k = 0; k < n; k++) t_div(r, NULL, v, k); return r; }

static series ad_sqr (series s, series u) { for (int k = 0; k < n; k++) s[k] = t_sqr(u, k); return s; }

static series ad_sqrt (series r, series u) { for (int k = 0; k < n; k++) t_sqrt(r, u, k); return r; }

static series ad_pwr (series p, series u, real a) { for (int k = 0; k < n; k++) t_pwr(p, u, a, k); return p; }

static series ad_exp (series e, series u) { for (int k = 0; k < n; k++) t_exp(e, u, k); return e; }

static void ad_sin_cos (series s, series c, series u, bool trig) { for (int k = 0; k < n; k++) t_sin_cos(s, c, u, k, trig); }

static void ad_tan_sec2 (series t, series s2, series u, bool trig) { for (int k = 0; k < n; k++) t_tan_sec2(t, s2, u, k, trig); }

static series ad_ln (series u, series e) { for (int k = 0; k < n; k++) t_ln(u, e, k); return u; }

static void ad_asin_cos (series u, series c, series s, bool trig) { for (int k = 0; k < n; k++) t_asin_cos(u, c, s, k, trig); }

static void ad_acos_sin (series u, series s, series c, bool trig) { for (int k = 0; k < n; k++) t_acos_sin(u, s, c, k, trig); }

static void ad_atan_sec2 (series u, series s, series t, bool trig) { for (int k = 0; k < n; k++) t_atan_sec2(u, s, t, k, trig); }

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
    total++;
    for (int k = 0; k < n; k++) {
        delta = fabsl(a[k] - b[k]);
        if (delta > delta_max) {
            delta_max = delta;
            name_max = name;
            k_max = k;
        }
        if (! isfinite(delta) || delta_max > tolerance) {
            fprintf(stderr, "%s FAIL%s %s%s%s\n", RED, NRM, WHT, name, NRM);
            return;
        }
        if (debug == 2) {
            if (!k) fprintf(stderr, "\n");
            fprintf(stderr, "  %2d  %s% .*Le % .*Le%s  %.1Le\n", k, GRY, dp, a[k], dp, b[k], NRM, delta);
        }
    }
    if (debug) fprintf(stderr, "%s PASS%s %s\n", GRN, NRM, name);
    passed++;
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
    n = (int)strtol(argv[2], NULL, BASE); CHECK(n > 0 && n <=64);
    series u = tsm_jet(n + 1);
    for (int k = 0, s = 1; k <= n; k++, s *= -1) {
        u[k] = !k ? strtold(argv[3], NULL) : 0.5L * s / SQR(k);
    }
    dual ud = {.val = u[0], .dot = u[1]};
    tolerance = strtold(argv[4], NULL); CHECK(tolerance > 0.0L);
    if (argc == 6) {
        debug = (int)strtol(argv[5], NULL, BASE); CHECK(debug == 0 || debug == 1 || debug == 2);
    }

    fprintf(stderr, "%sHorner Summation ", GRY);
    series s = tsm_jet(8);
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

    fprintf(stderr, "Taylor Arithmetic %su = %s%.1Lf%s\n", GRY, WHT, u[0], NRM);
    bool positive = u[0] > 0.0L, non_zero = u[0] != 0.0L, lt_1 = fabsl(u[0]) < 1.0L, gt_1 = u[0] > 1.0L,
         lt_pi_2 = fabsl(u[0]) < 0.5L * acosl(-1.0L);
    series r1 = tsm_jet(n), r2 = tsm_jet(n), r3 = tsm_jet(n), S1 = tsm_jet(n); S1[0] = 1.0L;
    series abs_u = tsm_jet(n), rec_u = tsm_jet(n), sqrt_u = tsm_jet(n), ln_u = tsm_jet(n);
    if (non_zero) ad_abs(abs_u, u);
    if (non_zero) ad_rec(rec_u, u);
    if (positive) ad_sqrt(sqrt_u, u);
    if (positive) ad_ln(ln_u, u);
    series sin_u = tsm_jet(n), cos_u = tsm_jet(n), tan_u = tsm_jet(n), sec2_u = tsm_jet(n);
    series asin_u = tsm_jet(n), acos_u = tsm_jet(n), atan_u = tsm_jet(n);
    series sin2_u = tsm_jet(n), cos2_u = tsm_jet(n), tan2_u = tsm_jet(n), sin_2u = tsm_jet(n), cos_2u = tsm_jet(n);
    ad_sin_cos(sin_u, cos_u, u, true);
    if (lt_pi_2) ad_tan_sec2(tan_u, sec2_u, u, true);
    ad_sqr(sin2_u, sin_u);
    ad_sqr(cos2_u, cos_u);
    ad_sin_cos(sin_2u, cos_2u, ad_scale(r1, u, D2), true);
    series sinh_u = tsm_jet(n), cosh_u = tsm_jet(n), tanh_u = tsm_jet(n), sech2_u = tsm_jet(n);
    series asinh_u = tsm_jet(n), acosh_u = tsm_jet(n), atanh_u = tsm_jet(n);
    series sinh2_u = tsm_jet(n), cosh2_u = tsm_jet(n), tanh2_u = tsm_jet(n), sinh_2u = tsm_jet(n), cosh_2u = tsm_jet(n);
    ad_sin_cos(sinh_u, cosh_u, u, false);
    ad_tan_sec2(tanh_u, sech2_u, u, false);
    ad_sqr(sinh2_u, sinh_u);
    ad_sqr(cosh2_u, cosh_u);
    ad_sin_cos(sinh_2u, cosh_2u, ad_scale(r1, u, D2), false);
    series sqr_u = tsm_jet(n), exp_u = tsm_jet(n), neg_exp_u = tsm_jet(n), gd_1 = tsm_jet(n);
    ad_sqr(sqr_u, u);
    ad_exp(exp_u, u);
    ad_exp(neg_exp_u, ad_scale(r1, u, D_1));
    ad_ln(gd_1, ad_abs(r3, ad_div(r2, ad_add(r1, sin_u, S1), cos_u)));

    char* name = "inv(u)"; non_zero ? compare_s_d(name, rec_u, d_rec(ud)) : skip(name);
    name = "sqr(u)"; compare_s_d(name, sqr_u, d_sqr(ud));
    name = "sqrt(u)"; positive ? compare_s_d(name, sqrt_u, d_sqrt(ud)) : skip(name);

    name = "exp(u)"; compare_s_d(name, exp_u, d_exp(ud));

    name = "sin(u)"; compare_s_d(name, sin_u, d_sin(ud));
    name = "cos(u)"; compare_s_d(name, cos_u, d_cos(ud));
    name = "tan(u)"; lt_pi_2 ? compare_s_d(name, tan_u, d_tan(ud)) : skip(name);

    name = "sinh(u)"; compare_s_d(name, sinh_u, d_sinh(ud));
    name = "cosh(u)"; compare_s_d(name, cosh_u, d_cosh(ud));
    name = "tanh(u)"; compare_s_d(name, tanh_u, d_tanh(ud));

    name = "ln(u)"; positive ? compare_s_d(name, ln_u, d_ln(ud)) : skip(name);

    name = "asin(u)"; if (lt_1) {ad_asin_cos(asin_u, r2, u, true); compare_s_d(name, asin_u, d_asin(ud));} else skip(name);
    name = "acos(u)"; if (lt_1) {ad_acos_sin(acos_u, r2, u, true); compare_s_d(name, acos_u, d_acos(ud));} else skip(name);
    name = "atan(u)"; ad_atan_sec2(atan_u, r2, u, true); compare_s_d(name, atan_u, d_atan(ud));

    name = "asinh(u)"; ad_asin_cos(asinh_u, r2, u, false); compare_s_d(name, asinh_u, d_asinh(ud));
    name = "acosh(u)"; if (gt_1) {ad_acos_sin(acosh_u, r2, u, false); compare_s_d(name, acosh_u, d_acosh(ud));} else skip(name);
    name = "atanh(u)"; if (lt_1) {ad_atan_sec2(atanh_u, r2, u, false); compare_s_d(name, atanh_u, d_atanh(ud));}  else skip(name);

    name = "u^1.5"; positive ? compare_s_d(name, ad_pwr(r1, u, 1.5L), d_pow(ud, 1.5L)) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "u * u == sqr(u)"; compare(name, ad_mul(r1, u, u), sqr_u);
    name = "u / u == 1"; non_zero ? compare(name, ad_div(r1, u, u), S1) : skip(name);
    name = "u * (1 / u) == 1"; non_zero ? compare(name, ad_mul(r1, u, rec_u), S1) : skip(name);
    name = "sqr(u) / u == u"; non_zero ? compare(name, ad_div(r1, sqr_u, u), u) : skip(name);
    name = "sqr(u) * (1 / u) == u"; non_zero ? compare(name, ad_mul(r1, sqr_u, rec_u), u) : skip(name);

    name = "sqrt(u) * sqrt(u) == u"; positive ? compare(name, ad_mul(r1, sqrt_u, sqrt_u), u) : skip(name);
    name = "u / sqrt(u) == sqrt(u)"; positive ? compare(name, ad_div(r1, u, sqrt_u), sqrt_u) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "u^2.0 == sqr(u)"; positive ? compare(name, ad_pwr(r1, u, D2), sqr_u) : skip(name);
    name = "u^1.0 == u"; positive ? compare(name, ad_pwr(r1, u, D1), u) : skip(name);
    name = "u^0.5 == sqrt(u)"; positive ? compare(name, ad_pwr(r1, u, D05), sqrt_u): skip(name);
    name = "u^0.0 == 1"; positive ? compare(name, ad_pwr(r1, u, D0), S1) : skip(name);
    name = "u^-0.5 == 1 / sqrt(u)"; positive ? compare(name, ad_pwr(r1, u, D_05), ad_rec(r2, sqrt_u)) : skip(name);
    name = "u^-1.0 == 1 / u"; positive ? compare(name, ad_pwr(r1, u, D_1), rec_u) : skip(name);
    name = "u^-2.0 == 1 / sqr(u)"; positive ? compare(name, ad_pwr(r1, u, D_2), ad_rec(r2, sqr_u)) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "sqr(u) * u^-3 == 1 / u"; positive ? compare(name, ad_mul(r1, sqr_u, ad_pwr(r2, u, D_3)), rec_u) : skip(name);
    name = "sqr(u)^0.5 == |u|"; non_zero ? compare(name, ad_pwr(r1, sqr_u, D05), abs_u) : skip(name);
    name = "sqrt(sqr(u) == |u|"; non_zero ? compare(name, ad_sqrt(r1, sqr_u), abs_u) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "ln(e^u) == u"; compare(name, ad_ln(r1, exp_u), u);
    name = "ln(sqr(u)) == ln(u) * 2"; positive ? compare(name, ad_ln(r1, sqr_u), ad_scale(r2, ln_u, D2)) : skip(name);
    name = "ln(sqrt(u)) == ln(u) / 2"; positive ? compare(name, ad_ln(r1, sqrt_u), ad_scale(r2, ln_u, D05)) : skip(name);
    name = "ln(1 / u) == -ln(u)"; positive ? compare(name, ad_ln(r1, rec_u), ad_scale(r2, ln_u, D_1)) : skip(name);
    name = "ln(u^-3) == -3ln(u)"; positive ? compare(name, ad_ln(r1, ad_pwr(r2, u, D_3)), ad_scale(r3, ln_u, D_3)) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "cosh^2(u) - sinh^2(u) == 1"; compare(name, ad_sub(r1, cosh2_u, sinh2_u), S1);
    name = "sech^2(u) + tanh^2(u) == 1"; compare(name, ad_add(r1, sech2_u, ad_sqr(tanh2_u, tanh_u)), S1);
    name = "tanh(u) == sinh(u) / cosh(u)"; compare(name, tanh_u, ad_div(r1, sinh_u, cosh_u));
    name = "sech^2(u) == 1 / cosh^2(u)"; compare(name, sech2_u, ad_rec(r1, cosh2_u));
    name = "sinh(2u) == 2 * sinh(u) * cosh(u)"; compare(name, sinh_2u, ad_scale(r1, ad_mul(r2, sinh_u, cosh_u), D2));
    name = "cosh(2u) == cosh^2(u) + sinh^2(u)"; compare(name, cosh_2u, ad_add(r1, cosh2_u, sinh2_u));

    if (debug) fprintf(stderr, "\n");

    name = "cosh(u) == (e^u + e^-u) / 2"; compare(name, cosh_u, ad_scale(r1, ad_add(r2, exp_u, neg_exp_u), D05));
    name = "sinh(u) == (e^u - e^-u) / 2"; compare(name, sinh_u, ad_scale(r1, ad_sub(r2, exp_u, neg_exp_u), D05));

    if (debug) fprintf(stderr, "\n");

    name = "arsinh(sinh(u)) == u"; ad_asin_cos(r1, r2, sinh_u, false); compare(name, r1, u);
    name = "arcosh(cosh(u)) == |u|"; if (non_zero) {ad_acos_sin(r1, r2, cosh_u, false); compare(name, r1, abs_u);} else skip(name);
    name = "artanh(tanh(u)) == u"; ad_atan_sec2(r1, r2, tanh_u, false); compare(name, r1, u);

    if (debug) fprintf(stderr, "\n");

    name = "cos^2(u) + sin^2(u) == 1"; compare(name, ad_add(r1, cos2_u, sin2_u), S1);
    name = "sec^2(u) - tan^2(u) == 1"; lt_pi_2 ? compare(name, ad_sub(r1, sec2_u, ad_sqr(tan2_u, tan_u)), S1) : skip(name);
    name = "tan(u) == sin(u) / cos(u)"; lt_pi_2 ? compare(name, tan_u, ad_div(r1, sin_u, cos_u)) : skip(name);
    name = "sec^2(u) == 1 / cos^2(u)"; lt_pi_2 ? compare(name, sec2_u, ad_rec(r1, cos2_u)) : skip(name);
    name = "sin(2u) == 2 * sin(u) * cos(u)"; compare(name, sin_2u, ad_scale(r1, ad_mul(r2, sin_u, cos_u), D2));
    name = "cos(2u) == cos^2(u) - sin^2(u)"; compare(name, cos_2u, ad_sub(r1, cos2_u, sin2_u));

    if (debug) fprintf(stderr, "\n");

    name = "arcsin(sin(u)) == u"; if (lt_pi_2) {ad_asin_cos(r1, r2, sin_u, true); compare(name, r1, u);} else skip(name);
    name = "arccos(cos(u)) == |u|"; if (non_zero) {ad_acos_sin(r1, r2, cos_u, true); compare(name, r1, abs_u);} else skip(name);
    name = "arctan(tan(u)) == u"; if (lt_pi_2) {ad_atan_sec2(r1, r2, tan_u, true); compare(name, r1, u);} else skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "arsinh(tan(u)) == gd^-1 u"; if (lt_pi_2) {ad_asin_cos(r1, r2, tan_u, false); compare(name, gd_1, r1);} else skip(name);
    name = "artanh(sin(u)) == gd^-1 u"; ad_atan_sec2(r1, r2, sin_u, false); compare(name, gd_1, r1);
    name = "arcsin(tanh(gd^-1 u)) == u"; if (lt_pi_2) {ad_tan_sec2(r3, r2, gd_1, false); ad_asin_cos(r1, r2, r3, true); compare(name, r1, u);} else skip(name);
    name = "arctan(sinh(gd^-1 u)) == u"; if (lt_pi_2) {ad_sin_cos(r3, r2, gd_1, false); ad_atan_sec2(r1, r2, r3, true); compare(name, r1, u);} else skip(name);

    if (debug) fprintf(stderr, "\n");
    fprintf(stderr, "%sTotal%s %d  %sPASSED%s %d", WHT, NRM, total, GRN, NRM, passed);
    if (skipped) fprintf(stderr, "  %sSKIPPED%s %d", YLW, NRM, skipped);
    if (passed == total - skipped) {
        fprintf(stderr, "\n%sDelta%s %.1Le %s%s%s %sk == %d%s\n", GRY, NRM, delta_max, BLU, name_max, NRM, GRY, k_max, NRM);
        return 0;
    } else {
        fprintf(stderr, "  %sFAILED%s %d\n\n", RED, NRM, total - passed - skipped);
        return 3;
    }
}

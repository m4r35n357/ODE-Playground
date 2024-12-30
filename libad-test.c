/*
 * Tests for core Taylor Series Method, Horner's Method, and recurrence relations
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
#include "taylor-ode.h"
#include "ad.h"

static int k_max = 0, dp, n, debug = 0, total = 0, passed = 0, skipped = 0;

static real delta, delta_max, tolerance, D0, D01, D05, D_05, D1, D_1, D2, D_2, D3, D_3;

static char *name_max = "N/A";

static void libad_test_init (void) {
    ad_init(n);
    mpfr_inits(delta, delta_max, NULL);
    mpfr_set_zero(delta_max, 1);
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

struct Parameters { real a, b, c; };

model *tsm_init_p (int argc, char **argv, int order) { (void)argc; (void)argv; (void)order;
    model *p = malloc(sizeof (model));
    mpfr_init_set_si(p->a, 1, RND);
    mpfr_init_set_si(p->b, 0, RND);
    mpfr_init_set_si(p->c, -1, RND);
    return p;
}

void ode (triplet *v, series x, series y, series z, model *p, int k) {
    mpfr_mul(v->x, p->a, x[k], RND);
    mpfr_mul(v->y, p->b, y[k], RND);
    mpfr_mul(v->z, p->c, z[k], RND);
}

static void skip (char* name) {
    total++;
    skipped++;
    if (debug) fprintf(stderr, "%s SKIP%s %s%s%s\n", YLW, NRM, GRY, name, NRM);
}

static void compare (char* name, series a, series b) {
    total++;
    for (int k = 0; k < n; k++) {
        mpfr_sub(delta, a[k], b[k], RND);
        mpfr_abs(delta, delta, RND);
        if (mpfr_cmp(delta, delta_max) > 0) {
            mpfr_set(delta_max, delta, RND);
            name_max = name;
            k_max = k;
        }
        if (mpfr_number_p(delta) == 0 || mpfr_cmpabs(delta, tolerance) > 0) {
            fprintf(stderr, "%s FAIL%s %s%s%s\n  k=%d  LHS: %+.*Le  RHS: %+.*Le  (%.1Le)\n",
                    RED, NRM, MGT, name, NRM, k, dp, mpfr_get_ld(a[k], RND), dp, mpfr_get_ld(b[k], RND), mpfr_get_ld(delta, RND));
            return;
        }
        if (debug >= 2) {
            if (!k) fprintf(stderr, "\n");
            fprintf(stderr, "%s  DEBUG%s  k: %2d  %+.*Le %+.*Le  (%.1Le)\n",
                    NRM, NRM, k, dp, mpfr_get_ld(a[k], RND), dp, mpfr_get_ld(b[k], RND), mpfr_get_ld(delta, RND));
        }
    }
    if (debug) fprintf(stderr, "%s PASS%s %s\n", GRN, NRM, name);
    passed++;
}

int main (int argc, char **argv) {
    PRINT_ARGS(argc, argv);
    CHECK(argc == 6 || argc == 7);

    real PI_2;
    dp = (int)strtol(argv[1], NULL, BASE);
    mpfr_set_default_prec((int)strtol(argv[2], NULL, BASE));
    n = (int)strtol(argv[3], NULL, BASE); CHECK(n > 8);
    libad_test_init();
    series x = tsm_jet(n + 1);
    mpfr_init_set_str(x[0], argv[4], BASE, RND);
    for (int k = 1; k <= n; k++) {
        mpfr_div_si(x[k], D05, k * k, RND);
    }
    mpfr_init_set_str(tolerance, argv[5], BASE, RND); CHECK(mpfr_sgn(tolerance) > 0);
    if (argc == 7) {
        debug = (int)strtol(argv[6], NULL, BASE); CHECK(debug == 0 || debug == 1 || debug == 2);
    }

    mpfr_init(PI_2);
    mpfr_const_pi(PI_2, RND);
    mpfr_div_2si(PI_2, PI_2, 1, RND);

    fprintf(stderr, "%sTaylor Series Method", GRY);
    int steps = 10;
    triplet *v = malloc(sizeof (triplet)); CHECK(v);
    mpfr_inits(v->x, v->y, v->z, NULL);
    xyz *_ = malloc(sizeof (xyz)); CHECK(_);
    _->x = tsm_jet(n + 1); mpfr_set(_->x[0], D1, RND);
    _->y = tsm_jet(n + 1); mpfr_set(_->y[0], D1, RND);
    _->z = tsm_jet(n + 1); mpfr_set(_->z[0], D1, RND);
    tsm_init(dp);
    tsm(n, D01, steps, v, _, tsm_init_p(argc, argv, n), clock());
    fprintf(stdout, "%sCheck: e^1  e^0  e^-1%s\n", WHT, NRM);
    real e1, e0, e_1;
    mpfr_inits(e1, e0, e_1, NULL);
    mpfr_exp(e1, D1, RND);
    mpfr_exp(e0, D0, RND);
    mpfr_exp(e_1, D_1, RND);
    _out_(e1, e0, e_1, D01, steps, 0.0F);

    fprintf(stderr, "%s, Horner Summation ", GRY);
    series p = tsm_jet(n);
    mpfr_set_si(p[0], 1, RND);
    mpfr_set_si(p[1], 3, RND);
    mpfr_set_si(p[2], 0, RND);
    mpfr_set_si(p[3], 2, RND);
    CHECK(!mpfr_cmp_ld(*horner(p, 3, D2), 23.0L)); fprintf(stderr, ".");
    mpfr_set_si(p[0], 3, RND);
    mpfr_set_si(p[1], -1, RND);
    mpfr_set_si(p[2], 2, RND);
    mpfr_set_si(p[3], -4, RND);
    mpfr_set_si(p[4], 0, RND);
    mpfr_set_si(p[5], 1, RND);
    CHECK(!mpfr_cmp_ld(*horner(p, 5, D3), 153.0L)); fprintf(stderr, ".");
    mpfr_set_si(p[0], 1, RND);
    mpfr_set_si(p[1], -4, RND);
    mpfr_set_si(p[2], 0, RND);
    mpfr_set_si(p[3], 0, RND);
    mpfr_set_si(p[4], 2, RND);
    mpfr_set_si(p[5], 3, RND);
    mpfr_set_si(p[6], 0, RND);
    mpfr_set_si(p[7], -2, RND);
    CHECK(!mpfr_cmp_ld(*horner(p, 7, D_2), 201.0L)); fprintf(stderr, ". %sOK%s\n", GRN, NRM);

    fprintf(stderr, "Taylor Arithmetic %sx = %s%.1Lf%s\n", GRY, WHT, mpfr_get_ld(x[0], RND), NRM);
    bool positive = mpfr_sgn(x[0]) > 0, non_zero = mpfr_zero_p(x[0]) == 0, lt_pi_2 = mpfr_cmpabs(x[0], PI_2) < 0;
    series r1 = tsm_jet(n), r2 = tsm_jet(n), r3 = tsm_jet(n), S1 = tsm_jet(n); mpfr_set(S1[0], D1, RND);
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
    ad_sin_cos(sin_2x, cos_2x, ad_scale(r1, x, D2), true);
    series sinh_x = tsm_jet(n), cosh_x = tsm_jet(n), tanh_x = tsm_jet(n), sech2_x = tsm_jet(n);
    series sinh2_x = tsm_jet(n), cosh2_x = tsm_jet(n), tanh2_x = tsm_jet(n), sinh_2x = tsm_jet(n), cosh_2x = tsm_jet(n);
    ad_sin_cos(sinh_x, cosh_x, x, false);
    ad_tan_sec2(tanh_x, sech2_x, x, false);
    ad_sqr(sinh2_x, sinh_x);
    ad_sqr(cosh2_x, cosh_x);
    ad_sin_cos(sinh_2x, cosh_2x, ad_scale(r1, x, D2), false);
    series sqr_x = tsm_jet(n), exp_x = tsm_jet(n), neg_exp_x = tsm_jet(n), gd_1 = tsm_jet(n);
    ad_sqr(sqr_x, x);
    ad_exp(exp_x, x);
    ad_exp(neg_exp_x, ad_scale(r1, x, D_1));
    ad_ln(gd_1, ad_abs(r3, ad_div(r2, ad_add(r1, sin_x, S1), cos_x)));

    char* name = "x * x == sqr(x)"; compare(name, ad_mul(r1, x, x), sqr_x);
    name = "sqr(x) / x == x"; non_zero ? compare(name, ad_div(r1, sqr_x, x), x) : skip(name);
    name = "x * (1 / x) == 1"; non_zero ? compare(name, ad_mul(r1, x, inv_x), S1) : skip(name);
    name = "sqrt(x) * sqrt(x) == x"; positive ? compare(name, ad_mul(r1, sqrt_x, sqrt_x), x) : skip(name);
    name = "x / sqrt(x) == sqrt(x)"; positive ? compare(name, ad_div(r1, x, sqrt_x), sqrt_x) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "x^2.0 == sqr(x)"; positive ? compare(name, ad_pwr(r1, x, D2), sqr_x) : skip(name);
    name = "x^1.0 == x"; positive ? compare(name, ad_pwr(r1, x, D1), x) : skip(name);
    name = "x^0.5 == sqrt(x)"; positive ? compare(name, ad_pwr(r1, x, D05), sqrt_x): skip(name);
    name = "x^0.0 == 1"; positive ? compare(name, ad_pwr(r1, x, D0), S1) : skip(name);
    name = "x^-0.5 == 1 / sqrt(x)"; positive ? compare(name, ad_pwr(r1, x, D_05), ad_div(r2, S1, sqrt_x)) : skip(name);
    name = "x^-1.0 == 1 / x"; positive ? compare(name, ad_pwr(r1, x, D_1), inv_x) : skip(name);
    name = "x^-2.0 == 1 / sqr(x)"; positive ? compare(name, ad_pwr(r1, x, D_2), ad_div(r2, S1, sqr_x)) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "sqr(x) * x^-3 == 1 / x"; positive ? compare(name, ad_mul(r1, sqr_x, ad_pwr(r2, x, D_3)), inv_x) : skip(name);
    name = "sqr(x)^0.5 == |x|"; positive ? compare(name, ad_pwr(r1, sqr_x, D05), abs_x) : skip(name);
    name = "sqrt(sqr(x) == |x|"; positive ? compare(name, ad_sqrt(r1, sqr_x), abs_x) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "ln(e^x) == x"; compare(name, ad_ln(r1, exp_x), x);
    name = "ln(sqr(x)) == ln(x) * 2"; positive ? compare(name, ad_ln(r1, sqr_x), ad_scale(r2, ln_x, D2)) : skip(name);
    name = "ln(sqrt(x)) == ln(x) / 2"; positive ? compare(name, ad_ln(r1, sqrt_x), ad_scale(r2, ln_x, D05)) : skip(name);
    name = "ln(1 / x) == -ln(x)"; positive ? compare(name, ad_ln(r1, inv_x), ad_scale(r2, ln_x, D_1)) : skip(name);
    name = "ln(x^-3) == -3ln(x)"; positive ? compare(name, ad_ln(r1, ad_pwr(r2, x, D_3)), ad_scale(r3, ln_x, D_3)) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "cosh^2(x) - sinh^2(x) == 1"; compare(name, ad_sub(r1, cosh2_x, sinh2_x), S1);
    name = "sech^2(x) + tanh^2(x) == 1"; compare(name, ad_add(r1, sech2_x, ad_sqr(tanh2_x, tanh_x)), S1);
    name = "tanh(x) == sinh(x) / cosh(x)"; compare(name, tanh_x, ad_div(r1, sinh_x, cosh_x));
    name = "sech^2(x) == 1 / cosh^2(x)"; compare(name, sech2_x, ad_div(r1, S1, cosh2_x));
    name = "sinh(2x) == 2 * sinh(x) * cosh(x)"; compare(name, sinh_2x, ad_scale(r1, ad_mul(r2, sinh_x, cosh_x), D2));
    name = "cosh(2x) == cosh^2(x) + sinh^2(x)"; compare(name, cosh_2x, ad_add(r1, cosh2_x, sinh2_x));

    if (debug) fprintf(stderr, "\n");

    name = "cosh(x) == (e^x + e^-x) / 2"; compare(name, cosh_x, ad_scale(r1, ad_add(r2, exp_x, neg_exp_x), D05));
    name = "sinh(x) == (e^x - e^-x) / 2"; compare(name, sinh_x, ad_scale(r1, ad_sub(r2, exp_x, neg_exp_x), D05));

    if (debug) fprintf(stderr, "\n");

    name = "arsinh(sinh(x)) == x"; ad_asin(r1, r2, sinh_x, false); compare(name, r1, x);
    name = "arcosh(cosh(x)) == |x|"; if (non_zero) {ad_acos(r1, r2, cosh_x, false); compare(name, r1, abs_x);} else skip(name);
    name = "artanh(tanh(x)) == x"; ad_atan(r1, r2, tanh_x, false); compare(name, r1, x);

    if (debug) fprintf(stderr, "\n");

    name = "cos^2(x) + sin^2(x) == 1"; compare(name, ad_add(r1, cos2_x, sin2_x), S1);
    name = "sec^2(x) - tan^2(x) == 1"; lt_pi_2 ? compare(name, ad_sub(r1, sec2_x, ad_sqr(tan2_x, tan_x)), S1) : skip(name);
    name = "tan(x) == sin(x) / cos(x)"; lt_pi_2 ? compare(name, tan_x, ad_div(r1, sin_x, cos_x)) : skip(name);
    name = "sec^2(x) == 1 / cos^2(x)"; lt_pi_2 ? compare(name, sec2_x, ad_div(r1, S1, cos2_x)) : skip(name);
    name = "sin(2x) == 2 * sin(x) * cos(x)"; compare(name, sin_2x, ad_scale(r1, ad_mul(r2, sin_x, cos_x), D2));
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
    fprintf(stderr, "%sTotal%s %d  %sPASSED%s %d", WHT, NRM, total, GRN, NRM, passed);
    if (skipped) fprintf(stderr, "  %sSKIPPED%s %d", YLW, NRM, skipped);
    if (passed == total - skipped) {
        fprintf(stderr, "\n%sDelta%s %.1Le %s%s%s %sk == %s%s\n", GRY, NRM, mpfr_get_ld(delta_max, RND), BLU, name_max, NRM, GRY, name_max, NRM);
        return 0;
    } else {
        fprintf(stderr, "  %sFAILED%s %d\n\n", RED, NRM, total - passed - skipped);
        return 3;
    }
}

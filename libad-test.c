/*
 * Tests for core Taylor Series Method, Horner's Method, and recurrence relations
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

static int k_max = 0, dp, n, debug = 0, total = 0, passed = 0, skipped = 0;

static real delta, delta_max, tolerance, D0, D01, D05, D_05, D1, D_1, D2, D_2, D3, D_3, PI_2;

static void libad_test_init (void) {
    mpfr_inits(delta, delta_max, PI_2, NULL);
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
    mpfr_const_pi(PI_2, RND);
    mpfr_div_2si(PI_2, PI_2, 1, RND);
}

static char *name_max = "N/A";

struct Parameters { real a, b, c; };

static series ad_scale (series s, series u, real a) { for (int k = 0; k < n; k++) mpfr_mul(s[k], u[k], a, RND); return s; }

static series ad_add (series p, series u, series v) { for (int k = 0; k < n; k++) mpfr_add(p[k], u[k], v[k], RND); return p; }

static series ad_sub (series m, series u, series v) { for (int k = 0; k < n; k++) mpfr_sub(m[k], u[k], v[k], RND); return m; }

static series ad_abs (series a, series u) { for (int k = 0; k < n; k++) mpfr_set(a[k], *t_abs(u, k), RND); return a; }

static series ad_mul (series p, series u, series v) { for (int k = 0; k < n; k++) mpfr_set(p[k], *t_mul(u, v, k), RND); return p; }

static series ad_div (series q, series u, series v) { for (int k = 0; k < n; k++) t_div(q, u, v, k); return q; }

static series ad_rec (series r, series v) { for (int k = 0; k < n; k++) t_div(r, NULL, v, k); return r; }

static series ad_sqr (series s, series u) { for (int k = 0; k < n; k++) mpfr_set(s[k], *t_sqr(u, k), RND); return s; }

static series ad_sqrt (series r, series u) { for (int k = 0; k < n; k++) t_sqrt(r, u, k); return r; }

static series ad_pwr (series p, series u, real a) { for (int k = 0; k < n; k++) t_pwr(p, u, a, k); return p; }

static series ad_exp (series e, series u) { for (int k = 0; k < n; k++) t_exp(e, u, k); return e; }

static void ad_sin_cos (series s, series c, series u, bool trig) { for (int k = 0; k < n; k++) t_sin_cos(s, c, u, k, trig); }

static void ad_tan_sec2 (series t, series s2, series u, bool trig) { for (int k = 0; k < n; k++) t_tan_sec2(t, s2, u, k, trig); }

static series ad_ln (series u, series e) { for (int k = 0; k < n; k++) t_ln(u, e, k); return u; }

static void ad_asin_cos (series u, series c, series s, bool trig) { for (int k = 0; k < n; k++) t_asin_cos(u, c, s, k, trig); }

static void ad_acos_sin (series u, series s, series c, bool trig) { for (int k = 0; k < n; k++) t_acos_sin(u, s, c, k, trig); }

static void ad_atan_sec2 (series u, series s, series t, bool trig) { for (int k = 0; k < n; k++) t_atan_sec2(u, s, t, k, trig); }

model *tsm_init_p (int argc, char **argv, int o) { (void)argc; (void)argv; (void)o;
    model *p = malloc(sizeof (model));
    mpfr_init_set_si(p->a, 1, RND);
    mpfr_init_set_si(p->b, 0, RND);
    mpfr_init_set_si(p->c, -1, RND);
    return p;
}

void ode (triplet *v, series x, series y, series z, const model *p, int k) {
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

    dp = (int)strtol(argv[1], NULL, BASE);
    mpfr_set_default_prec((int)strtol(argv[2], NULL, BASE));
    n = (int)strtol(argv[3], NULL, BASE); CHECK(n > 0 && n <=64);
    libad_test_init();
    series u = tsm_jet(n + 1);
    mpfr_init_set_str(u[0], argv[4], BASE, RND);
    for (int k = 1; k <= n; k++) {
        mpfr_mul_si(u[k], u[k - 1], -1, RND);
    }
    mpfr_init_set_str(tolerance, argv[5], BASE, RND); CHECK(mpfr_sgn(tolerance) > 0);
    if (argc == 7) {
        debug = (int)strtol(argv[6], NULL, BASE); CHECK(debug == 0 || debug == 1 || debug == 2);
    }

    fprintf(stderr, "%sHorner Summation ", GRY);
    tsm_init(dp);
    series p = tsm_jet(8);
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

    fprintf(stderr, "%sTaylor Series Method . . . %s\n", GRY, NRM);
    int steps = 10;
    triplet *v = malloc(sizeof (triplet)); CHECK(v);
    mpfr_inits(v->x, v->y, v->z, NULL);
    series x = tsm_jet(n + 1); mpfr_set(x[0], D1, RND);
    series y = tsm_jet(n + 1); mpfr_set(y[0], D1, RND);
    series z = tsm_jet(n + 1); mpfr_set(z[0], D1, RND);
    tsm(n, D01, steps, v, x, y, z, tsm_init_p(argc, argv, n), clock());
    fprintf(stdout, "%sCheck: e^1  e^0  e^-1%s\n", WHT, NRM);
    real e1, e0, e_1;
    mpfr_inits(e1, e0, e_1, NULL);
    mpfr_exp(e1, D1, RND);
    mpfr_exp(e0, D0, RND);
    mpfr_exp(e_1, D_1, RND);
    _out_(e1, e0, e_1, D01, steps, 0.0F);

    fprintf(stderr, "Taylor Arithmetic %su = %s%.1Lf%s\n", GRY, WHT, mpfr_get_ld(u[0], RND), NRM);
    bool positive = mpfr_sgn(u[0]) > 0, non_zero = mpfr_zero_p(u[0]) == 0, lt_pi_2 = mpfr_cmpabs(u[0], PI_2) < 0;
    series r1 = tsm_jet(n), r2 = tsm_jet(n), r3 = tsm_jet(n), S1 = tsm_jet(n); mpfr_set(S1[0], D1, RND);
    series abs_u = tsm_jet(n), rec_u = tsm_jet(n), sqrt_u = tsm_jet(n), ln_u = tsm_jet(n);
    if (non_zero) ad_abs(abs_u, u);
    if (non_zero) ad_rec(rec_u, u);
    if (positive) ad_sqrt(sqrt_u, u);
    if (positive) ad_ln(ln_u, u);
    series sin_u = tsm_jet(n), cos_u = tsm_jet(n), tan_u = tsm_jet(n), sec2_u = tsm_jet(n);
    series sin2_u = tsm_jet(n), cos2_u = tsm_jet(n), tan2_u = tsm_jet(n), sin_2u = tsm_jet(n), cos_2u = tsm_jet(n);
    ad_sin_cos(sin_u, cos_u, u, true);
    if (lt_pi_2) ad_tan_sec2(tan_u, sec2_u, u, true);
    ad_sqr(sin2_u, sin_u);
    ad_sqr(cos2_u, cos_u);
    ad_sin_cos(sin_2u, cos_2u, ad_scale(r1, u, D2), true);
    series sinh_u = tsm_jet(n), cosh_u = tsm_jet(n), tanh_u = tsm_jet(n), sech2_u = tsm_jet(n);
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

    char* name = "u * u == sqr(u)"; compare(name, ad_mul(r1, u, u), sqr_u);
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

    name = "cosh(u) == (e^u + e^-u) / 2"; compare(name, cosh_u, ad_scale(r1, ad_add(r1, exp_u, neg_exp_u), D05));
    name = "sinh(u) == (e^u - e^-u) / 2"; compare(name, sinh_u, ad_scale(r1, ad_sub(r2, exp_u, neg_exp_u), D05));
    name = "tanh(u) == (e^u - e^-u) / (e^u + e^-u)"; compare(name, tanh_u, ad_div(r3, ad_sub(r2, exp_u, neg_exp_u), ad_add(r1, exp_u, neg_exp_u)));

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
        fprintf(stderr, "\n%sDelta%s %.1Le %s%s%s %sk == %d%s\n", GRY, NRM, mpfr_get_ld(delta_max, RND), BLU, name_max, NRM, GRY, k_max, NRM);
        return 0;
    } else {
        fprintf(stderr, "  %sFAILED%s %d\n\n", RED, NRM, total - passed - skipped);
        return 3;
    }
}

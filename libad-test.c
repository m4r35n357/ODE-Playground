/*
 * Automatic Differentiation of Taylor Series, recurrence ralations validation checks
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "taylor-ode.h"

#define NRM "\x1B[0;37m"
#define WHT "\x1B[1;37m"
#define GRN "\x1B[1;32m"
#define YLW "\x1B[1;33m"
#define RED "\x1B[1;31m"
#define CYN "\x1B[0;36m"

static int n, debug = 0, total = 0, passed = 0, skipped = 0;

static real delta, tolerance;

static const real PLUS1 = 1.0L, ZERO = 0.0L, MINUS1 = -1.0L;

typedef struct Parameters { real a, b, c; } parameters;

static series ad_scale (series s, series u, real a) {
    for (int k = 0; k < n; k++) {
        s[k] = u[k] * a;
    }
    return s;
}

static series ad_add (series p, series u, series v) {
    for (int k = 0; k < n; k++) {
        p[k] = u[k] + v[k];
    }
    return p;
}

static series ad_sub (series m, series u, series v) {
    for (int k = 0; k < n; k++) {
        m[k] = u[k] - v[k];
    }
    return m;
}

static series ad_abs (series a, series u) {
    for (int k = 0; k < n; k++) {
        a[k] = t_abs(u, k);
    }
    return a;
}

static series ad_mul (series p, series u, series v) {
    for (int k = 0; k < n; k++) {
        p[k] = t_mul(u, v, k);
    }
    return p;
}

static series ad_div (series q, series u, series v) {
    for (int k = 0; k < n; k++) {
        t_div(q, u, v, k);
    }
    return q;
}

static series ad_inv (series i, series v) {
    for (int k = 0; k < n; k++) {
        t_div(i, NULL, v, k);
    }
    return i;
}

static series ad_sqr (series s, series u) {
    for (int k = 0; k < n; k++) {
        s[k] = t_sqr(u, k);
    }
    return s;
}

static series ad_sqrt (series r, series u) {
    for (int k = 0; k < n; k++) {
        t_sqrt(r, u, k);
    }
    return r;
}

static series ad_exp (series e, series u) {
    for (int k = 0; k < n; k++) {
        t_exp(e, u, k);
    }
    return e;
}

static void ad_sin_cos (series s, series c, series u, bool trig) {
    for (int k = 0; k < n; k++) {
        t_sin_cos(s, c, u, k, trig);
    }
}

static void ad_tan_sec2 (series t, series s2, series u, bool trig) {
    for (int k = 0; k < n; k++) {
        t_tan_sec2(t, s2, u, k, trig);
    }
}

static series ad_pwr (series p, series u, real a) {
    for (int k = 0; k < n; k++) {
        t_pwr(p, u, a, k);
    }
    return p;
}

static series ad_ln (series l, series u) {
    for (int k = 0; k < n; k++) {
        t_ln(l, u, k);
    }
    return l;
}

static void ad_asin (series as, series du_df, series u, bool trig) {
    for (int k = 0; k < n; k++) {
        t_asin(as, du_df, u, k, trig);
    }
}

static void ad_acos (series ac, series du_df, series u, bool trig) {
    for (int k = 0; k < n; k++) {
        t_acos(ac, du_df, u, k, trig);
    }
}

static void ad_atan (series at, series du_df, series u, bool trig) {
    for (int k = 0; k < n; k++) {
        t_atan(at, du_df, u, k, trig);
    }
}

void *get_p (int argc, char **argv, int order) { (void)argc; (void)argv; (void)order;
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
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
    if (debug) fprintf(stderr, "%s SKIP%s %s\n", YLW, NRM, name);
}

static void compare (char* name, series a, series b) {
    total++;
    for (int k = 0; k < n; k++) {
        delta = a[k] - b[k];
        if (!isfinite(delta) || fabsl(delta) > tolerance) {
            fprintf(stderr, "%s FAIL%s %s\n  k=%d  LHS: %+.6Le  RHS: %+.6Le  (%+.3Le)\n", RED, NRM, name, k, a[k], b[k], delta);
            return;
        }
        if (debug >= 2) {
            if (!k) fprintf(stderr, "\n");
            fprintf(stderr, "%s  DEBUG%s  k: %2d  %+.6Le %+.6Le  (%+.3Le)\n", NRM, NRM, k, a[k], b[k], delta);
        }
    }
    if (debug) fprintf(stderr, "%s PASS%s %s\n", GRN, NRM, name);
    passed++;
}

int main (int argc, char **argv) {
    PRINT_ARGS(argc, argv);
    CHECK(argc == 4 || argc == 5);

    n = (int)strtol(argv[1], NULL, BASE); CHECK(n > 0);
    series x = t_jet(n + 1);
    x[0] = strtold(argv[2], NULL);
    for (int k = 1; k <= n; k++) {
        x[k] = x[0] / (k * k);
    }
    tolerance = strtold(argv[3], NULL); CHECK(tolerance > 0.0L);
    if (argc == 5) {
        debug = (int)strtol(argv[4], NULL, BASE); CHECK(debug == 0 || debug == 1 || debug == 2);
     }

    bool positive = x[0] > 0.0L, non_zero = x[0] != 0.0L, lt_pi_2 = fabsl(x[0]) < 0.5L * acosl(-1.0L);

    series abs_x = t_jet(n), inv_x = t_jet(n), sqr_x = t_jet(n), sqrt_x = t_jet(n);
    series sin2_x = t_jet(n), cos2_x = t_jet(n), tan2_x = t_jet(n);
    series sinh2_x = t_jet(n), cosh2_x = t_jet(n), tanh2_x = t_jet(n);
    series exp_x = t_jet(n), neg_exp_x = t_jet(n), ln_x = t_jet(n);
    series sin_x = t_jet(n), sin_2x = t_jet(n), sinh_x = t_jet(n), sinh_2x = t_jet(n);
    series cos_x = t_jet(n), cos_2x = t_jet(n), cosh_x = t_jet(n), cosh_2x = t_jet(n);
    series tan_x = t_jet(n), sec2_x = t_jet(n), tanh_x = t_jet(n), sech2_x = t_jet(n);
    series gd_1 = t_jet(n), r1 = t_jet(n), r2 = t_jet(n), r3 = t_jet(n);

    series S1 = t_const(n, 1.0L);

    fprintf(stdout, "%sHorner Summation%s\n", WHT, NRM);
    series s = t_jet(n >= 8 ? n : 8);
    s[0] = 1.0L; s[1] = 3.0L; s[2] = 0.0L; s[3] = 2.0L;
    fprintf(stdout, " 23 %8.3Lf\n", t_horner(s, 3, 2.0L));
    s[0] = 3; s[1] = -1.0L; s[2] = 2.0L; s[3] = -4.0L; s[4] = 0.0L; s[5] = 1.0L;
    fprintf(stdout, "153 %8.3Lf\n", t_horner(s, 5, 3.0L));
    s[0] = 1.0L; s[1] = -4.0L; s[2] = 0.0L; s[3] = 0.0L; s[4] = 2.0L; s[5] = 3.0L; s[6] = 0.0L; s[7] = -2.0L;
    fprintf(stdout, "201 %8.3Lf\n", t_horner(s, 7, -2.0L));

    int dp = 12;
    controls c = {.order=n, .step=0, .steps=10, .step_size=0.1L};
    void *p = get_p(argc, argv, n);
    series3 *j = malloc(sizeof (series3)); CHECK(j);

    fprintf(stdout, "%sTaylor Series Method (stdout): x'=1  y'=0  z'=-1%s\n", WHT, NRM);
    j->x = t_const(n + 1, 1.0L);
    j->y = t_const(n + 1, 1.0L);
    j->z = t_const(n + 1, 1.0L);
    tsm_stdout(dp, &c, j, p, clock());

    fprintf(stdout, "%sTaylor Series Method (generator): x'=1  y'=0  z'=-1%s\n", WHT, NRM);
    j->x = t_const(n + 1, 1.0L);
    j->y = t_const(n + 1, 1.0L);
    j->z = t_const(n + 1, 1.0L);
    fprintf(stdout, "%+.*Le %+.*Le %+.*Le %.6Le\n", dp, j->x[0], dp, j->y[0], dp, j->z[0], 0.0L);
    while (tsm_gen(&c, j, p)) {
        fprintf(stdout, "%+.*Le %+.*Le %+.*Le %.6Le\n", dp, j->x[0], dp, j->y[0], dp, j->z[0], (c.step + 1) * c.step_size);
    }

    fprintf(stdout, "%sTaylor Series Method Check: e^1  e^0  e^-1%s\n", WHT, NRM);
    t_out(dp, expl(PLUS1), expl(ZERO), expl(MINUS1), c.step_size * c.steps, "_", "_", "_", clock());

    fprintf(stderr, "%sRecurrence Relations: %s%sx = %.1Lf%s\n", WHT, NRM, CYN, x[0], NRM);

    ad_sqr(sqr_x, x);
    if (non_zero) ad_inv(inv_x, x);
    if (positive) ad_sqrt(sqrt_x, x);

    char* name = "x * x == sqr(x)"; compare(name, ad_mul(r1, x, x), sqr_x);

    name = "sqr(x) / x == x"; non_zero ? compare(name, ad_div(r1, sqr_x, x), x) : skip(name);

    name = "x * 1 / x == 1"; non_zero ? compare(name, ad_mul(r1, x, inv_x), S1) : skip(name);

    name = "sqrt(x) * sqrt(x) == x"; positive ? compare(name, ad_mul(r1, sqrt_x, sqrt_x), x) : skip(name);

    name = "x / sqrt(x) == sqrt(x)"; positive ? compare(name, ad_div(r1, x, sqrt_x), sqrt_x) : skip(name);

    if (debug) fprintf(stderr, "\n");

    name = "x^2.0 == sqr(x)"; positive ? compare(name, ad_pwr(r1, x, 2.0L), sqr_x) : skip(name);

    name = "x^1.0 == x"; positive ? compare(name, ad_pwr(r1, x, 1.0L), x) : skip(name);

    name = "x^0.5 == sqrt(x)"; positive ? compare(name, ad_pwr(r1, x, 0.5L), sqrt_x): skip(name);

    name = "x^0.0 == 1"; positive ? compare(name, ad_pwr(r1, x, 0.0L), S1) : skip(name);

    name = "x^-0.5 == 1 / sqrt(x)"; positive ? compare(name, ad_pwr(r1, x, -0.5L), ad_inv(r2, sqrt_x)) : skip(name);

    name = "x^-1.0 == 1 / x"; positive ? compare(name, ad_pwr(r1, x, -1.0L), inv_x) : skip(name);

    name = "x^-2.0 == 1 / sqr(x)"; positive ? compare(name, ad_pwr(r1, x, -2.0L), ad_inv(r2, sqr_x)) : skip(name);

    if (debug) fprintf(stderr, "\n");
    ad_abs(abs_x, x);

    name = "sqr(x) * x^-3 == 1 / x"; positive ? compare(name, ad_mul(r1, sqr_x, ad_pwr(r2, x, -3.0L)), inv_x) : skip(name);

    name = "sqr(x)^0.5 == |x|"; non_zero ? compare(name, ad_pwr(r1, sqr_x, 0.5L), abs_x) : skip(name);

    name = "sqrt(sqr(x) == |x|"; non_zero ? compare(name, ad_sqrt(r1, sqr_x), abs_x) : skip(name);

    if (debug) fprintf(stderr, "\n");
    ad_exp(exp_x, x);
    if (positive) ad_ln(ln_x, x);

    name = "ln(e^x) == x"; compare(name, ad_ln(r1, exp_x), x);

    name = "ln(sqr(x)) == ln(x) * 2"; positive ? compare(name, ad_ln(r1, sqr_x), ad_scale(r2, ln_x, 2.0L)) : skip(name);

    name = "ln(sqrt(x)) == ln(x) / 2"; positive ? compare(name, ad_ln(r1, sqrt_x), ad_scale(r2, ln_x, 0.5L)) : skip(name);

    name = "ln(1 / x) == -ln(x)"; positive ? compare(name, ad_ln(r1, inv_x), ad_scale(r2, ln_x, -1.0L)) : skip(name);

    name = "ln(x^-3) == -3ln(x)"; positive ? compare(name, ad_ln(r1, ad_pwr(r2, x, -3.0L)), ad_scale(r3, ln_x, -3.0L)) : skip(name);

    if (debug) fprintf(stderr, "\n");
    ad_sin_cos(sinh_x, cosh_x, x, false);
    ad_tan_sec2(tanh_x, sech2_x, x, false);
    ad_sqr(sinh2_x, sinh_x);
    ad_sqr(cosh2_x, cosh_x);
    ad_sin_cos(sinh_2x, cosh_2x, ad_scale(r1, x, 2.0L), false);

    name = "cosh^2(x) - sinh^2(x) == 1"; compare(name, ad_sub(r1, cosh2_x, sinh2_x), S1);

    name = "sech^2(x) + tanh^2(x) == 1"; compare(name, ad_add(r1, sech2_x, ad_sqr(tanh2_x, tanh_x)), S1);

    name = "tanh(x) == sinh(x) / cosh(x)"; compare(name, tanh_x, ad_div(r1, sinh_x, cosh_x));

    name = "sech^2(x) == 1 / cosh^2(x)"; compare(name, sech2_x, ad_inv(r1, cosh2_x));

    name = "sinh(2x) == 2 * sinh(x) * cosh(x)"; compare(name, sinh_2x, ad_scale(r1, ad_mul(r2, sinh_x, cosh_x), 2.0L));

    name = "cosh(2x) == cosh^2(x) + sinh^2(x)"; compare(name, cosh_2x, ad_add(r1, cosh2_x, sinh2_x));

    if (debug) fprintf(stderr, "\n");
    ad_exp(neg_exp_x, ad_scale(r1, x, -1.0L));

    name = "cosh(x) == (e^x + e^-x) / 2"; compare(name, cosh_x, ad_scale(r1, ad_add(r2, exp_x, neg_exp_x), 0.5L));

    name = "sinh(x) == (e^x - e^-x) / 2"; compare(name, sinh_x, ad_scale(r1, ad_sub(r2, exp_x, neg_exp_x), 0.5L));

    if (debug) fprintf(stderr, "\n");

    name = "arcsinh(sinh(x)) == x"; ad_asin(r1, r2, sinh_x, false); compare(name, r1, x);

    name = "arccosh(cosh(x)) == |x|"; ad_acos(r1, r2, cosh_x, false); non_zero ? compare(name, r1, abs_x) : skip(name);

    name = "arctanh(tanh(x)) == x"; ad_atan(r1, r2, tanh_x, false); compare(name, r1, x);

    if (debug) fprintf(stderr, "\n");
    ad_sin_cos(sin_x, cos_x, x, true);
    ad_tan_sec2(tan_x, sec2_x, x, true);
    ad_sqr(sin2_x, sin_x);
    ad_sqr(cos2_x, cos_x);
    ad_sin_cos(sin_2x, cos_2x, ad_scale(r1, x, 2.0L), true);

    name = "cos^2(x) + sin^2(x) == 1"; compare(name, ad_add(r1, cos2_x, sin2_x), S1);

    name = "sec^2(x) - tan^2(x) == 1"; lt_pi_2 ? compare(name, ad_sub(r1, sec2_x, ad_sqr(tan2_x, tan_x)), S1) : skip(name);

    name = "tan(x) == sin(x) / cos(x)"; lt_pi_2 ? compare(name, tan_x, ad_div(r1, sin_x, cos_x)) : skip(name);

    name = "sec^2(x) == 1 / cos^2(x)"; lt_pi_2 ? compare(name, sec2_x, ad_inv(r1, cos2_x)) : skip(name);

    name = "sin(2x) == 2 * sin(x) * cos(x)"; compare(name, sin_2x, ad_scale(r1, ad_mul(r2, sin_x, cos_x), 2.0L));

    name = "cos(2x) == cos^2(x) - sin^2(x)"; compare(name, cos_2x, ad_sub(r1, cos2_x, sin2_x));

    if (debug) fprintf(stderr, "\n");

    name = "arcsin(sin(x)) == x"; ad_asin(r1, r2, sin_x, true); compare(name, r1, x);

    name = "arccos(cos(x)) == |x|"; ad_acos(r1, r2, cos_x, true); non_zero ? compare(name, r1, abs_x) : skip(name);

    name = "arctan(tan(x)) == x"; ad_atan(r1, r2, tan_x, true); compare(name, r1, x);

    if (debug) fprintf(stderr, "\n");
    ad_ln(gd_1, ad_abs(r3, ad_div(r2, ad_add(r1, sin_x, S1), cos_x)));

    name = "arsin(tan(x)) == gd^-1 x"; ad_asin(r1, r2, tan_x, false); compare(name, gd_1, r1);

    name = "artan(sin(x)) == gd^-1 x"; ad_atan(r1, r2, sin_x, false); compare(name, gd_1, r1);

    name = "arcsin(tanh(gd^-1 x)) == x"; ad_tan_sec2(r3, r2, gd_1, false); ad_asin(r1, r2, r3, true); compare(name, r1, x);

    name = "arctan(sinh(gd^-1 x)) == x"; ad_sin_cos(r3, r2, gd_1, false); ad_atan(r1, r2, r3, true); compare(name, r1, x);

    if (debug) fprintf(stderr, "\n");
    fprintf(stderr, "%sTotal%s: %d, %sPASSED%s %d", WHT, NRM, total, GRN, NRM, passed);
    if (skipped) fprintf(stderr, ", %sSKIPPED%s %d", YLW, NRM, skipped);
    if (passed == total - skipped) {
        fprintf(stderr, "\n");
        return 0;
    } else {
        fprintf(stderr, ", %sFAILED%s %d\n\n", RED, NRM, total - passed - skipped);
        return 3;
    }
}

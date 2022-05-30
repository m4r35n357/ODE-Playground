/*
 * Automatic Differentiation of Taylor Series, newest validation checks
 *
 * Example: ./libad-test-dbg 237 64 .5 1e-64 [ 0 | 1 | 2 ]
 *
./libad-test-dbg $(yad --title="Taylor Arithmetic Tests" --form --separator=" " --align=right \
    --field="Precision in Bits":NUM --field="Order":NUM --field="Value":NUM --field="Deviation":NUM --field="Detail":CB \
    -- '237!11..999!2' '64!4..256!1' '0.5!-1.0..1.0!0.1!1' '18!3..36!3' '0!1!2')
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"
#include "ad.h"

#define NRM "\x1B[0;37m"
#define WHT "\x1B[1;37m"
#define GRN "\x1B[1;32m"
#define YLW "\x1B[1;33m"
#define RED "\x1B[1;31m"
#define CYN "\x1B[0;36m"

static int n, debug = 0, total = 0, passed = 0, skipped = 0;

static mpfr_t delta, tolerance, D0, D01, D05, D_05, D1, D_1, D2, D_2, D3, D_3;

static void libad_test_init (void) {
    ad_init(n);
    t_init(12);
    mpfr_init(delta);
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

typedef struct { mpfr_t a, b, c; } parameters;

void *get_p (int argc, char **argv, int order) {
    (void)argc; (void)argv; (void)order;
    parameters *p = malloc(sizeof (parameters));
    mpfr_init_set(p->a, D1, RND);
    mpfr_init_set(p->b, D0, RND);
    mpfr_init_set(p->c, D_1, RND);
    return p;
}

void ode (components *v, series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    mpfr_mul(v->x, p->a, x[k], RND);
    mpfr_mul(v->y, p->b, y[k], RND);
    mpfr_mul(v->z, p->c, z[k], RND);
}

static void skip (char* name) {
    total++;
    skipped++;
    if (debug >= 1) fprintf(stderr, "%s SKIP%s %s\n", YLW, NRM, name);
}

static void compare (char* name, series a, series b) {
    total++;
    for (int k = 0; k < n; k++) {
        mpfr_sub(delta, a[k], b[k], RND);
        if (mpfr_number_p(delta) == 0 || mpfr_cmpabs(delta, tolerance) > 0) {
            fprintf(stderr, "%s FAIL%s %s  k=%d  LHS: %.6e  RHS: %.6e  (%.3e)\n",
                    RED, NRM, name, k, mpfr_get_d(a[k], RND), mpfr_get_d(b[k], RND), mpfr_get_d(delta, RND));
            return;
        }
        if (debug >= 2) {
            if (k == 0) fprintf(stderr, "\n");
            fprintf(stderr, "%s  DEBUG%s  k: %2d  %+.6e %+.6e  (%+.3e)\n",
                    NRM, NRM, k, mpfr_get_d(a[k], RND), mpfr_get_d(b[k], RND), mpfr_get_d(delta, RND));
        }
    }
    if (debug >= 1) fprintf(stderr, "%s PASS%s %s\n", GRN, NRM, name);
    passed++;
}

int main (int argc, char **argv) {
    mpfr_t PI_2;

    assert(argc == 5 || argc == 6);
    mpfr_set_default_prec((int)strtol(argv[1], NULL, BASE));
    n = (int)strtol(argv[2], NULL, BASE);
    assert(n > 1);
    libad_test_init();
    series x = t_jet(n + 1);
    mpfr_init_set_str(x[0], argv[3], BASE, RND);
    for (int k = 1; k <= n; k++) {
        mpfr_div_si(x[k], x[0], k * k, RND);
    }
    mpfr_init_set_str(tolerance, argv[4], BASE, RND);
    if (argc == 6) debug = (int)strtol(argv[5], NULL, BASE);

    mpfr_init(PI_2);
    mpfr_const_pi(PI_2, RND);
    mpfr_div_2si(PI_2, PI_2, 1, RND);

    int x_positive = mpfr_sgn(x[0]) > 0;
    int x_non_zero = mpfr_zero_p(x[0]) == 0;
    int x_lt_pi_2 = mpfr_cmpabs(x[0], PI_2) < 0;

    series r1 = t_jet(n);
    series r2 = t_jet(n);
    series r3 = t_jet(n);
    series abs_x = t_jet(n);
    series inv_x = t_jet(n);
    series sqr_x = t_jet(n);
    series sqr_sin_x = t_jet(n);
    series sqr_cos_x = t_jet(n);
    series sqr_tan_x = t_jet(n);
    series sqr_sinh_x = t_jet(n);
    series sqr_cosh_x = t_jet(n);
    series sqr_tanh_x = t_jet(n);
    series sqrt_x = t_jet(n);
    series exp_x = t_jet(n);
    series neg_exp_x = t_jet(n);
    series ln_x = t_jet(n);
    series sin = t_jet(n);
    series sin_2x = t_jet(n);
    series cos = t_jet(n);
    series cos_2x = t_jet(n);
    series tan = t_jet(n);
    series sec2 = t_jet(n);
    series sinh = t_jet(n);
    series sinh_2x = t_jet(n);
    series cosh = t_jet(n);
    series cosh_2x = t_jet(n);
    series tanh = t_jet(n);
    series sech2 = t_jet(n);
    series gd_1 = t_jet(n);

    series S1 = ad_const(t_jet(n), D1);

    fprintf(stdout, "\n");
    fprintf(stdout, "%sHorner%s\n", WHT, NRM);
    series p = t_jet(n >= 7 ? n : 7);
    mpfr_set_si(p[0], 1, RND);
    mpfr_set_si(p[1], 3, RND);
    mpfr_set_si(p[2], 0, RND);
    mpfr_set_si(p[3], 2, RND);
    mpfr_fprintf(stdout, " 23 %8.3RNf\n", *t_horner(p, 3, D2));
    mpfr_set_si(p[0], 3, RND);
    mpfr_set_si(p[1], -1, RND);
    mpfr_set_si(p[2], 2, RND);
    mpfr_set_si(p[3], -4, RND);
    mpfr_set_si(p[4], 0, RND);
    mpfr_set_si(p[5], 1, RND);
    mpfr_fprintf(stdout, "153 %8.3RNf\n", *t_horner(p, 5, D3));
    mpfr_set_si(p[0], 1, RND);
    mpfr_set_si(p[1], -4, RND);
    mpfr_set_si(p[2], 0, RND);
    mpfr_set_si(p[3], 0, RND);
    mpfr_set_si(p[4], 2, RND);
    mpfr_set_si(p[5], 3, RND);
    mpfr_set_si(p[6], 0, RND);
    mpfr_set_si(p[7], -2, RND);
    mpfr_fprintf(stdout, "201 %8.3RNf\n", *t_horner(p, 7, D_2));

    fprintf(stdout, "\n");
    fprintf(stdout, "%sTaylor Series Method: x'=1  y'=0  z'=-1%s\n", WHT, NRM);
    int steps = 10;
    tsm(n, D01, steps, D1, D1, D1, get_p(argc, argv, n));
    fprintf(stdout, "%sCheck: e^1  e^0  e^-1%s\n", WHT, NRM);
    mpfr_t e1, e0, e_1;
    mpfr_inits(e1, e0, e_1, NULL);
    mpfr_exp(e1, D1, RND);
    mpfr_exp(e0, D0, RND);
    mpfr_exp(e_1, D_1, RND);
    t_output(e1, e0, e_1, D01, steps, 0.0);

    fprintf(stderr, "\n");
    fprintf(stderr, "%sRecurrence Relations: %s%sx = %.1Lf%s\n", WHT, NRM, CYN, mpfr_get_ld(x[0], RND), NRM);

    ad_sqr(sqr_x, x);
    if (x_non_zero) ad_inv(inv_x, x);
    if (x_positive) ad_sqrt(sqrt_x, x);

    char* name = "x * x == sqr(x)";
    compare(name, ad_mul(r1, x, x), sqr_x);

    name = "sqr(x) / x == x";
    x_non_zero ? compare(name, ad_div(r1, sqr_x, x), x) : skip(name);

    name = "x * 1 / x == 1";
    x_non_zero ? compare(name, ad_mul(r1, x, inv_x), S1) : skip(name);

    name = "sqrt(x) * sqrt(x) == x";
    x_positive ? compare(name, ad_mul(r1, sqrt_x, sqrt_x), x) : skip(name);

    name = "x / sqrt(x) == sqrt(x)";
    x_positive ? compare(name, ad_div(r1, x, sqrt_x), sqrt_x) : skip(name);

    if (debug != 0) fprintf(stderr, "\n");

    name = "x^2.0 == sqr(x)";
    x_positive ? compare(name, ad_pwr(r1, x, D2), sqr_x) : skip(name);

    name = "x^1.0 == x";
    x_positive ? compare(name, ad_pwr(r1, x, D1), x) : skip(name);

    name = "x^0.5 == sqrt(x)";
    x_positive ? compare(name, ad_pwr(r1, x, D05), sqrt_x): skip(name);

    name = "x^0.0 == 1";
    x_positive ? compare(name, ad_pwr(r1, x, D0), S1) : skip(name);

    name = "x^-0.5 == 1 / sqrt(x)";
    x_positive ? compare(name, ad_pwr(r1, x, D_05), ad_inv(r2, sqrt_x)) : skip(name);

    name = "x^-1.0 == 1 / x";
    x_positive ? compare(name, ad_pwr(r1, x, D_1), inv_x) : skip(name);

    name = "x^-2.0 == 1 / sqr(x)";
    x_positive ? compare(name, ad_pwr(r1, x, D_2), ad_inv(r2, sqr_x)) : skip(name);

    if (debug != 0) fprintf(stderr, "\n");

    name = "x^2 == sqr(x)";
    x_positive ? compare(name, ad_ipwr(r1, x, 2), sqr_x) : skip(name);

    name = "x^1 == x";
    x_positive ? compare(name, ad_ipwr(r1, x, 1), x) : skip(name);

    name = "x^0 == 1";
    x_positive ? compare(name, ad_ipwr(r1, x, 0), S1) : skip(name);

    name = "x^-1 == 1 / x";
    x_positive ? compare(name, ad_ipwr(r1, x, -1), inv_x) : skip(name);

    name = "x^-2 == 1 / sqr(x)";
    x_positive ? compare(name, ad_ipwr(r1, x, -2), ad_inv(r2, sqr_x)) : skip(name);

    if (debug != 0) fprintf(stderr, "\n");
    ad_abs(abs_x, x);

    name = "sqr(x) * x^-3 == 1 / x";
    x_positive ? compare(name, ad_mul(r1, sqr_x, ad_ipwr(r2, x, -3)), inv_x) : skip(name);

    name = "sqr(x)^0.5 == |x|";
    x_non_zero ? compare(name, ad_pwr(r1, sqr_x, D05), abs_x) : skip(name);

    name = "sqrt(sqr(x) == |x|";
    x_non_zero ? compare(name, ad_sqrt(r1, sqr_x), abs_x) : skip(name);

    if (debug != 0) fprintf(stderr, "\n");
    ad_exp(exp_x, x);
    if (x_positive) ad_ln(ln_x, x);

    name = "log(e^x) == x";
    compare(name, ad_ln(r1, exp_x), x);

    name = "log(sqr(x)) == log(x) * 2";
    x_positive ? compare(name, ad_ln(r1, sqr_x), ad_scale(r2, ln_x, D2)) : skip(name);

    name = "log(sqrt(x)) == log(x) / 2";
    x_positive ? compare(name, ad_ln(r1, sqrt_x), ad_scale(r2, ln_x, D05)) : skip(name);

    name = "log(1 / x) == - log(x)";
    x_positive ? compare(name, ad_ln(r1, inv_x), ad_scale(r2, ln_x, D_1)) : skip(name);

    name = "log(x^-2) == - 2 * log(x)";
    x_positive ? compare(name, ad_ln(r1, ad_ipwr(r2, x, -2)), ad_scale(r3, ln_x, D_2)) : skip(name);

    if (debug != 0) fprintf(stderr, "\n");
    ad_sin_cos(sinh, cosh, x, HYP);
    ad_tan_sec2(tanh, sech2, x, HYP);
    ad_sqr(sqr_sinh_x, sinh);
    ad_sqr(sqr_cosh_x, cosh);
    ad_sin_cos(sinh_2x, cosh_2x, ad_scale(r1, x, D2), HYP);

    name = "cosh^2(x) - sinh^2(x) == 1";
    compare(name, ad_sub(r1, sqr_cosh_x, sqr_sinh_x), S1);

    name = "sech^2(x) + tanh^2(x) == 1";
    compare(name, ad_add(r1, sech2, ad_sqr(sqr_tanh_x, tanh)), S1);

    name = "tanh(x) == sinh(x) / cosh(x)";
    compare(name, tanh, ad_div(r1, sinh, cosh));

    name = "sech^2(x) == 1 / cosh^2(x)";
    compare(name, sech2, ad_inv(r1, sqr_cosh_x));

    name = "sinh(2x) == 2 * sinh(x) * cosh(x)";
    compare(name, sinh_2x, ad_scale(r1, ad_mul(r2, sinh, cosh), D2));

    name = "cosh(2x) == cosh^2(x) + sinh^2(x)";
    compare(name, cosh_2x, ad_add(r1, sqr_cosh_x, sqr_sinh_x));

    if (debug != 0) fprintf(stderr, "\n");
    ad_exp(neg_exp_x, ad_scale(r1, x, D_1));

    name = "cosh(x) == (e^x + e^-x) / 2";
    compare(name, cosh, ad_scale(r1, ad_add(r2, exp_x, neg_exp_x), D05));

    name = "sinh(x) == (e^x - e^-x) / 2";
    compare(name, sinh, ad_scale(r1, ad_sub(r2, exp_x, neg_exp_x), D05));

    if (debug != 0) fprintf(stderr, "\n");

    ad_asin(r1, r2, sinh, HYP);
    name = "arcsinh(sinh(x)) == x";
    compare(name, r1, x);

    ad_acos(r1, r2, cosh, HYP);
    name = "arccosh(cosh(x)) == |x|";
    x_non_zero ? compare(name, r1, abs_x) : skip(name);

    ad_atan(r1, r2, tanh, HYP);
    name = "arctanh(tanh(x)) == x";
    compare(name, r1, x);

    if (debug != 0) fprintf(stderr, "\n");
    ad_sin_cos(sin, cos, x, TRIG);
    ad_tan_sec2(tan, sec2, x, TRIG);
    ad_sqr(sqr_sin_x, sin);
    ad_sqr(sqr_cos_x, cos);
    ad_sin_cos(sin_2x, cos_2x, ad_scale(r1, x, D2), TRIG);

    name = "cos^2(x) + sin^2(x) == 1";
    compare(name, ad_add(r1, sqr_cos_x, sqr_sin_x), S1);

    name = "sec^2(x) - tan^2(x) == 1";
    x_lt_pi_2 ? compare(name, ad_sub(r1, sec2, ad_sqr(sqr_tan_x, tan)), S1) : skip(name);

    name = "tan(x) == sin(x) / cos(x)";
    x_lt_pi_2 ? compare(name, tan, ad_div(r1, sin, cos)) : skip(name);

    name = "sec^2(x) == 1 / cos^2(x)";
    x_lt_pi_2 ? compare(name, sec2, ad_inv(r1, sqr_cos_x)) : skip(name);

    name = "sin(2x) == 2 * sin(x) * cos(x)";
    compare(name, sin_2x, ad_scale(r1, ad_mul(r2, sin, cos), D2));

    name = "cos(2x) == cos^2(x) - sin^2(x)";
    compare(name, cos_2x, ad_sub(r1, sqr_cos_x, sqr_sin_x));

    if (debug != 0) fprintf(stderr, "\n");

    ad_asin(r1, r2, sin, TRIG);
    name = "arcsin(sin(x)) == x";
    compare(name, r1, x);

    ad_acos(r1, r2, cos, TRIG);
    name = "arccos(cos(x)) == |x|";
    x_non_zero ? compare(name, r1, abs_x) : skip(name);

    ad_atan(r1, r2, tan, TRIG);
    name = "arctan(tan(x)) == x";
    compare(name, r1, x);

    if (debug != 0) fprintf(stderr, "\n");
    ad_ln(gd_1, ad_abs(r3, ad_div(r2, ad_add(r1, sin, S1), cos)));

    name = "arsin(tan(x)) == gd^-1 x";
    ad_asin(r1, r2, tan, HYP);
    compare(name, gd_1, r1);

    name = "artan(sin(x)) == gd^-1 x";
    ad_atan(r1, r2, sin, HYP);
    compare(name, gd_1, r1);

    ad_tan_sec2(r3, r2, gd_1, HYP);
    ad_asin(r1, r2, r3, TRIG);
    name = "arcsin(tanh(gd^-1 x)) == x";
    compare(name, r1, x);

    ad_sin_cos(r3, r2, gd_1, HYP);
    ad_atan(r1, r2, r3, TRIG);
    name = "arctan(sinh(gd^-1 x)) == x";
    compare(name, r1, x);

    if (debug != 0) fprintf(stderr, "\n");
    fprintf(stderr, "%sTotal%s: %d, %sPASSED%s %d", WHT, NRM, total, GRN, NRM, passed);
    if (skipped > 0) {
        fprintf(stderr, ", %sSKIPPED%s %d", YLW, NRM, skipped);
    }
    if (passed == total - skipped) {
        fprintf(stderr, "\n\n");
        return 0;
    } else {
        fprintf(stderr, ", %sFAILED%s %d\n\n", RED, NRM, total - passed - skipped);
        return 3;
    }
}

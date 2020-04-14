/*
 * Automatic Differentiation of Taylor Series, newest validation checks
 *
 * Example: ./libad-test-dbg 32 20 2 1e-18
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"
#include "ad.h"

#define KNRM "\x1B[0;37m"
#define KGRN "\x1B[1;32m"
#define KYLW "\x1B[1;33m"
#define KRED "\x1B[1;31m"

typedef enum {PASS, FAIL} result;

static int total = 0, passed = 0, skipped = 0;

static mpfr_t delta;

static result compare (char* name, mpfr_t *a, mpfr_t *b, int size, mpfr_t threshold) {
    total++;
    for (int k = 0; k < size; k++) {
        mpfr_sub(delta, a[k], b[k], RND);
        if (mpfr_cmp_abs(delta, threshold) > 0) {
            printf("%sFAILED%s %s  k: %d  diff: %.3e  a: %.3e  b: %.3e\n",
                    KRED, KNRM, name, k, mpfr_get_d(delta, RND), mpfr_get_d(a[k], RND), mpfr_get_d(b[k], RND));
            return FAIL;
        }
    }
    printf("%sPASSED%s %s\n", KGRN, KNRM, name);
    passed++;
    return PASS;
}

int main (int argc, char **argv) {
    mpfr_t x0, DA, D05, D1, tol;

    assert(argc == 5);

    mpfr_init(delta);
    mpfr_set_default_prec(strtod(argv[1], NULL) * 3.322);
    ad_tempvars();
    long n = strtol(argv[2], NULL, BASE) + 1;
    assert(n > 1);
    mpfr_init_set_str(x0, argv[3], BASE, RND);
    mpfr_init_set_str(tol, argv[4], BASE, RND);

    double da = -3.7;
    mpfr_init_set_d(DA, da, RND);
    mpfr_init_set_str(D05, "0.5", BASE, RND);
    mpfr_init_set_ui(D1, 1, RND);

    mpfr_t *a = t_jet(n);
    mpfr_t *p = t_jet(n);
    mpfr_t *q = t_jet(n);
    mpfr_t *m = t_jet(n);
    mpfr_t *s1 = t_jet(n);
    mpfr_t *s2 = t_jet(n);
    mpfr_t *p1 = t_jet(n);
    mpfr_t *p2 = t_jet(n);
    mpfr_t *p3 = t_jet(n);
    mpfr_t *e1 = t_jet(n);
    mpfr_t *e2 = t_jet(n);
    mpfr_t *l1 = t_jet(n);
    mpfr_t *l2 = t_jet(n);
    mpfr_t *s = t_jet(n);
    mpfr_t *c = t_jet(n);
    mpfr_t *t = t_jet(n);

    printf("\n");

    mpfr_t *c1 = t_jet_c(n, D1);
    mpfr_t *x = t_jet_c(n, x0);
    set_ad_status(x, VARIABLE);

    compare("x * x == sqr(x)", ad_prod(p, x, x, n), ad_sqr(s1, x, n), n, tol);
    if (mpfr_sgn(x[0]) > 0) {
        compare("x * x == x^2", p, ad_pwr(p1, x, 2.0, n), n, tol);
    } else skipped++;

    if (mpfr_sgn(x[0]) > 0) {
        compare("x^1 == x", ad_pwr(p1, x, 1.0, n), x, n, tol);
    } else skipped++;
    if (mpfr_sgn(x[0]) > 0) {
        compare("x^0.5 == sqrt(x)", ad_pwr(p1, x, 0.5, n), ad_sqrt(s, x, n), n, tol);
    } else skipped++;
    if (mpfr_sgn(x[0]) > 0) {
        compare("x^0 == x / x", ad_pwr(p1, x, 0.0, n), ad_quot(q, x, x, n), n, tol);
    } else skipped++;
    if (mpfr_sgn(x[0]) > 0) {
        compare("x^-0.5 == 1 / sqrt(x)", ad_pwr(p1, x, -0.5, n), ad_quot(q, c1, ad_sqrt(s, x, n), n), n, tol);
    } else skipped++;
    if (mpfr_sgn(x[0]) > 0) {
        compare("x^-1 == 1 / x", ad_pwr(p1, x, -1.0, n), ad_quot(q, c1, x, n), n, tol);
    } else skipped++;
    if (mpfr_sgn(x[0]) > 0) {
        compare("x^2 * x^-5 == x^-3", ad_prod(p, ad_pwr(p1, x, 2.0, n), ad_pwr(p2, x, -5.0, n), n), ad_pwr(p3, x, -3.0, n), n, tol);
    } else skipped++;

    if (mpfr_zero_p(x[0]) == 0) {
        compare("sqrt(x * x) == |x|", ad_sqrt(s, ad_prod(p, x, x, n), n), ad_abs(a, x, n), n, tol);
    } else skipped++;

    if (mpfr_sgn(x[0]) > 0) {
        compare("log(x^a) == a * log(x)", ad_ln(l1, ad_pwr(p1, x, da, n), n), ad_scale(s, ad_ln(l2, x, n), DA, n), n, tol);
    } else skipped++;
    compare("log(e^x) == x", ad_ln(l1, ad_exp(e1, x, n), n), x, n, tol);

    ad_sin_cos(s, c, x, n, TRIG);
    compare("cos^2(x) + sin^2(x) == 1", ad_plus(p, ad_sqr(s1, c, n), ad_sqr(s2, s, n), n), c1, n, tol);

    ad_sin_cos(s, c, x, n, HYP);
    compare("cosh^2(x) - sinh^2(x) == 1", ad_minus(m, ad_sqr(s1, c, n), ad_sqr(s2, s, n), n), c1, n, tol);

    ad_exp(e1, x, n);
    ad_exp(e2, ad_neg(m, x, n), n);
    compare("sinh(x) == 0.5 * (e^x - e^-x)", s, ad_scale(s, ad_minus(m, e1, e2, n), D05, n), n, tol);
    compare("cosh(x) == 0.5 * (e^x + e^-x)", c, ad_scale(s, ad_plus(p, e1, e2, n), D05, n), n, tol);

    ad_tan_sec2(t, s2, x, n, TRIG);
    compare("sec^2(x) - tan^2(x) == 1", ad_minus(m, s2, ad_sqr(s1, t, n), n), c1, n, tol);

    ad_sin_cos(s, c, x, n, TRIG);
    compare("tan(x) == sin(x) / cos(x)", t, ad_quot(q, s, c, n), n, tol);
    compare("sec^2(x) == 1 / cos^2(x)", s2, ad_quot(q, c1, ad_sqr(s1, c, n), n), n, tol);

    ad_tan_sec2(t, s2, x, n, HYP);
    compare("sech^2(x) + tanh^2(x) == 1", ad_plus(p, s2, ad_sqr(s1, t, n), n), c1, n, tol);

    ad_sin_cos(s, c, x, n, HYP);
    compare("tanh(x) == sinh(x) / cosh(x)", t, ad_quot(q, s, c, n), n, tol);
    compare("sech^2(x) == 1 / cosh^2(x)", s2, ad_quot(q, c1, ad_sqr(s1, c, n), n), n, tol);

    printf("Total: %d, Passed: %d", total, passed);
    if (skipped > 0) {
        printf(", %sSKIPPED%s %d", KYLW, KNRM, skipped);
    }
    if (passed == total) {
        printf("\n\n");
        return 0;
    } else {
        printf(", %sFAILED%s %d\n\n", KRED, KNRM, total - passed);
        return 1;
    }
}

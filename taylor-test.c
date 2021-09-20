/*
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 *
 * Example:  grep -n '[^a-zA-Z]t_[a-z0-9_]*' taylor-test.c
 *
 * Example:  f_list='t_const t_abs t_prod t_sqr t_quot t_inv t_sqrt t_exp t_sin_cos t_tan_sec2 t_pwr t_ln'
 *
 * Example:  for f in $f_list; do echo $f; grep -n $f taylor-test.c tsm-*.c; done
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "taylor-ode.h"

#define BBLK "\x1B[0;40m"
#define KNRM "\x1B[0;37m"
#define KCYN "\x1B[36m"
#define KGRY "\x1B[2;37m"
#define KBLD "\x1B[1;37m"

static const real P_A = 1.0L, P_B = 0.0L, P_C = -1.0L, X0 = 1.0L, Y0 = 1.0L, Z0 = 1.0L;

typedef struct { real a, b, c; } parameters;

void *get_p (int argc, char **argv, long order) {
    (void)argc;
    (void)argv;
    (void)order;
    parameters *p = malloc(sizeof (parameters));
    p->a = P_A;
    p->b = P_B;
    p->c = P_C;
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

static void output (long dp, long n, series x, series y, series z) {
    char fs[128];
    sprintf(fs, " %s%s[%%2ld]%s%s%s %%+.%ldLe %%+.%ldLe %%+.%ldLe%s\n", KBLD, KGRY, KNRM, BBLK, KGRY, dp, dp, dp, KNRM);
    for (int k = 0; k <= n; k++) {
        printf(fs, k, x[k], y[k], z[k]);
    }
}

int main (int argc, char **argv) {
    assert(argc == 5);
    long dp = strtol(argv[1], NULL, 10);
    long n = strtol(argv[2], NULL, 10);
    real x = strtold(argv[3], NULL); assert(x > 0.0L);
    real y = strtold(argv[4], NULL); assert(y != 0.0L);
    series in_X = t_jet(n + 1);
    in_X[0] = x;
    in_X[1] = 1.0L;
    series in_Y = t_jet(n + 1);
    in_Y[0] = y;
    in_Y[1] = 1.0L;
    series in_HYP = t_jet(n + 1);
    in_HYP[0] = MY_PI / x;
    in_HYP[1] = 1.0L;
    series in_TRIG = t_jet(n + 1);
    in_TRIG[0] = MY_PI / y;
    in_TRIG[1] = 1.0L;
    series out1 = t_jet(n + 1);
    series out2 = t_jet(n + 1);
    series out3 = t_jet(n + 1);
    series out4 = t_jet(n + 1);
    series out5 = t_jet(n + 1);
    series out6 = t_jet(n + 1);
    series target = t_jet(n + 1);

    printf("%s%s%s\n", KCYN, "Horner's Method, should = 128.000", KNRM);
    out1[0] = -19.0L; out1[1] = 7.0L; out1[2] = -4.0L; out1[3] = 6.0L;
    t_horner(out1, 3, 3.0L);
    printf("%.3Lf\n",  out1[0]);
    printf("\n");

    printf("%s%s%s\n", KCYN, "Taylor Series Method, should be e^1.0, e^0.0, e^-1.0", KNRM);
    tsm(argc, argv, 12, 10, 0.1L, 10, X0, Y0, Z0);
    printf("%s%s+%.12Le %+.12Le %+.12Le%s\n", KBLD, KGRY, expl(P_A), expl(P_B), expl(P_C), KNRM);
    printf("\n");

    printf("%s%s%s\n", KCYN, "Additive constant", KNRM);
    for (int k = 0; k <= n; k++) {
        out1[k] = t_const(y, k);
        out2[k] = t_const(-y, k);
        target[k] = out1[k] + out2[k];
    }
    output(dp, n, out1, out2, target);
    printf("\n");

    printf("%s%.1Lf * %.1Lf = %.1Lf%s\n", KCYN, in_X[0], in_Y[0], in_X[0] * in_Y[0], KNRM);
    for (int k = 0; k <= n; k++) {
        target[k] = t_prod(in_X, in_Y, k);
    }
    output(dp, n, in_X, in_Y, target);

    printf("%s%.1Lf / %.1Lf = %.1Lf%s\n", KCYN, in_X[0], in_Y[0], in_X[0] / in_Y[0], KNRM);
    for (int k = 0; k <= n; k++) {
        t_quot(target, in_X, in_Y, k);
    }
    output(dp, n, in_X, in_Y, target);

    printf("%s|%.1Lf| * |%.1Lf| = |%.1Lf * %.1Lf|%s\n", KCYN, in_X[0], in_Y[0], in_X[0], in_Y[0], KNRM);
    for (int k = 0; k <= n; k++) {
        out1[k] = t_abs(in_X, k);
        out2[k] = t_abs(in_Y, k);
        target[k] = t_prod(out1, out2, k);
    }
    output(dp, n, out1, out2, target);

    printf("%s|%.1Lf| / |%.1Lf| = |%.1Lf / %.1Lf|%s\n", KCYN, in_X[0], in_Y[0], in_X[0], in_Y[0], KNRM);
    for (int k = 0; k <= n; k++) {
        out1[k] = t_abs(in_X, k);
        out2[k] = t_abs(in_Y, k);
        t_quot(target, out1, out2, k);
    }
    output(dp, n, out1, out2, target);

    printf("%s%.1Lf / %.1Lf = 1%s\n", KCYN, in_Y[0], in_Y[0], KNRM);
    for (int k = 0; k <= n; k++) {
        t_quot(target, in_Y, in_Y, k);
    }
    output(dp, n, in_Y, in_Y, target);

    printf("%s%.1Lf * (1 / %.1Lf) = 1%s\n", KCYN, in_X[0], in_X[0], KNRM);
    for (int k = 0; k <= n; k++) {
        t_inv(out1, in_X, k);
        target[k] = t_prod(in_X, out1, k);
    }
    output(dp, n, in_X, out1, target);

    printf("%s%.1Lf * (%.1Lf / %.1Lf) = %.1Lf%s\n", KCYN, in_X[0], in_Y[0], in_X[0], in_Y[0], KNRM);
    for (int k = 0; k <= n; k++) {
        t_quot(out1, in_Y, in_X, k);
        target[k] = t_prod(in_X, out1, k);
    }
    output(dp, n, in_X, out1, target);

    printf("%s(sqr(%.1Lf))^0.5 = |%.1Lf|%s\n", KCYN, in_Y[0], in_Y[0], KNRM);
    for (int k = 0; k <= n; k++) {
        out1[k] = t_sqr(in_Y, k);
        t_pwr(target, out1, 0.5L, k);
    }
    output(dp, n, in_Y, out1, target);

    printf("%ssqrt(%.1Lf^2) = |%.1Lf|%s\n", KCYN, in_X[0], in_X[0], KNRM);
    for (int k = 0; k <= n; k++) {
        t_pwr(out1, in_X, 2.0L, k);
        t_sqrt(target, out1, k);
    }
    output(dp, n, in_X, out1, target);

    printf("%s1 / (%.1Lf^-1) = %.1Lf%s\n", KCYN, in_X[0], in_X[0], KNRM);
    for (int k = 0; k <= n; k++) {
        t_pwr(out1, in_X, -1.0L, k);
        t_inv(target, out1, k);
    }
    output(dp, n, in_X, out1, target);

    printf("%sexp(ln(%.1Lf) = %.1Lf%s\n", KCYN, in_X[0], in_X[0], KNRM);
    for (int k = 0; k <= n; k++) {
        t_ln(out1, in_X, k);
        t_exp(target, out1, k);
    }
    output(dp, n, in_X, out1, target);

    printf("%sln(exp(%.1Lf) = %.1Lf%s\n", KCYN, in_Y[0], in_Y[0], KNRM);
    for (int k = 0; k <= n; k++) {
        t_exp(out1, in_Y, k);
        t_ln(target, out1, k);
    }
    output(dp, n, in_Y, out1, target);

    printf("%ssin^2(%.3Lf) + cos^2(%.3Lf) = 1.0%s\n", KCYN, in_TRIG[0], in_TRIG[0], KNRM);
    for (int k = 0; k <= n; k++) {
        t_sin_cos(out1, out2, in_TRIG, k, TRIG);
        out3[k] = t_prod(out1, out1, k);
        out4[k] = t_sqr(out2, k);
        target[k] = out3[k] + out4[k];
    }
    output(dp, n, out3, out4, target);

    printf("%scosh^2(%.3Lf) - sinh^2(%.3Lf) = 1.0%s\n", KCYN, in_HYP[0], in_HYP[0], KNRM);
    for (int k = 0; k <= n; k++) {
        t_sin_cos(out1, out2, in_HYP, k, HYP);
        out3[k] = t_prod(out1, out1, k);
        out4[k] = t_sqr(out2, k);
        target[k] = out4[k] - out3[k];
    }
    output(dp, n, out4, out3, target);

    printf("%ssec^2(%.3Lf) - tanh^2(%.3Lf) = 1.0%s\n", KCYN, in_TRIG[0], in_TRIG[0], KNRM);
    for (int k = 0; k <= n; k++) {
        t_tan_sec2(out1, out2, in_TRIG, k, TRIG);
        out3[k] = t_sqr(out1, k);
        target[k] = out2[k] - out3[k];
    }
    output(dp, n, out2, out3, target);

    printf("%stanh^2(%.3Lf) + sech^2(%.3Lf) = 1.0%s\n", KCYN, in_HYP[0], in_HYP[0], KNRM);
    for (int k = 0; k <= n; k++) {
        t_tan_sec2(out1, out2, in_HYP, k, HYP);
        out3[k] = t_sqr(out1, k);
        target[k] = out3[k] + out2[k];
    }
    output(dp, n, out3, out2, target);

    printf("%ssin(%.3Lf) / cos(%.3Lf) - tan(%.3Lf) = 0.0%s\n", KCYN, in_TRIG[0], in_TRIG[0], in_TRIG[0], KNRM);
    for (int k = 0; k <= n; k++) {
        t_sin_cos(out1, out2, in_TRIG, k, TRIG);
        t_tan_sec2(out3, out4, in_TRIG, k, TRIG);
        t_quot(out5, out1, out2, k);
        target[k] = out5[k] - out3[k];
    }
    output(dp, n, out5, out3, target);

    printf("%ssinh(%.3Lf) / cosh(%.3Lf) - tanh(%.3Lf) = 0.0%s\n", KCYN, in_HYP[0], in_HYP[0], in_HYP[0], KNRM);
    for (int k = 0; k <= n; k++) {
        t_sin_cos(out1, out2, in_HYP, k, HYP);
        t_tan_sec2(out3, out4, in_HYP, k, HYP);
        t_quot(out5, out1, out2, k);
        target[k] = out5[k] - out3[k];
    }
    output(dp, n, out5, out3, target);

    printf("%s1.0 / cos(%.3Lf) - sqrt(sec^2(%.3Lf)) = 0.0%s\n", KCYN, in_TRIG[0], in_TRIG[0], KNRM);
    for (int k = 0; k <= n; k++) {
        t_sin_cos(out1, out2, in_TRIG, k, TRIG);
        t_tan_sec2(out3, out4, in_TRIG, k, TRIG);
        t_inv(out5, out2, k);
        t_sqrt(out6, out4, k);
        target[k] = out5[k] - out6[k];
    }
    output(dp, n, out5, out6, target);

    printf("%s1.0 / cosh(%.3Lf) - sqrt(sech^2(%.3Lf)) = 0.0%s\n", KCYN, in_HYP[0], in_HYP[0], KNRM);
    for (int k = 0; k <= n; k++) {
        t_sin_cos(out1, out2, in_HYP, k, HYP);
        t_tan_sec2(out3, out4, in_HYP, k, HYP);
        t_inv(out5, out2, k);
        t_sqrt(out6, out4, k);
        target[k] = out5[k] - out6[k];
    }
    output(dp, n, out5, out6, target);

    return 0;
}

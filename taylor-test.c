/*
 *  c99 -g -O3 -Wall  -o test test.c taylor-ode.c -lm
 *
 *  ./taylor-test 6 6 4 -3
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

#define BBLK "\x1B[0;40m"
#define KNRM "\x1B[0;37m"
#define KCYN "\x1B[36m"
#define KGRY "\x1B[2;37m"
#define KBLD "\x1B[1;37m"

static void to_derivs (series derivs, series jet, long n) {
    for (int i = 0; i <= n; i++) {
        derivs[i] = jet[i];
    }
    for (int i = 2; i <= n; i++) {
        for (int j = i; j <= n; j++)
            derivs[j] *= i;
    }
}

static void output (long dp, long n, series d1, series d2, series d3, series x, series y, series z) {
    to_derivs(d1, x, n);
    to_derivs(d2, y, n);
    to_derivs(d3, z, n);
    char fs[128];
    sprintf(fs, "%s%s[%%2ld]%s%s%s %%+.%ldLe %%+.%ldLe %%+.%ldLe%s\n", KBLD, KGRY, KNRM, BBLK, KGRY, dp, dp, dp, KNRM);
    for (int k = 0; k <= n; k++) {
		printf(fs, k, d1[k], d2[k], d3[k]);
	}
}

int main (int argc, char **argv) {
    assert(argc == 5);
    long dp = strtol(argv[1], NULL, 10);
    long n = strtol(argv[2], NULL, 10);
    real x = strtold(argv[3], NULL); assert(x > 0.0L);
    real y = strtold(argv[4], NULL); assert(y != 0.0L);
    series in = t_jet(n + 1);
    in[0] = x;
    in[1] = 1.0L;
    series in_B = t_jet(n + 1);
    in_B[0] = y;
    in_B[1] = 1.0L;
    series out1 = t_jet(n + 1);
    series out2 = t_jet(n + 1);
    series out3 = t_jet(n + 1);
    series d1 = t_jet(n + 1);
    series d2 = t_jet(n + 1);
    series d3 = t_jet(n + 1);

    printf("%s%.1Lf * %.1Lf = %.1Lf%s\n", KCYN, x, y, x * y, KNRM);
    for (int k = 0; k <= n; k++) {
        out3[k] = t_prod(in, in_B, k);
    }
    output(dp, n, d1, d2, d3, in, in_B, out3);

    printf("%s%.1Lf / %.1Lf = %.1Lf%s\n", KCYN, x, y, x / y, KNRM);
    for (int k = 0; k <= n; k++) {
        t_quot(out3, in, in_B, k);
    }
    output(dp, n, d1, d2, d3, in, in_B, out3);

    printf("%s|%.1Lf| * |%.1Lf| = |%.1Lf * %.1Lf|%s\n", KCYN, x, y, x, y, KNRM);
    for (int k = 0; k <= n; k++) {
        out1[k] = t_abs(in, k);
        out2[k] = t_abs(in_B, k);
        out3[k] = t_prod(out1, out2, k);
    }
    output(dp, n, d1, d2, d3, out1, out2, out3);

    printf("%s|%.1Lf| / |%.1Lf| = |%.1Lf / %.1Lf|%s\n", KCYN, x, y, x, y, KNRM);
    for (int k = 0; k <= n; k++) {
        out1[k] = t_abs(in, k);
        out2[k] = t_abs(in_B, k);
        t_quot(out3, out1, out2, k);
    }
    output(dp, n, d1, d2, d3, out1, out2, out3);

    printf("%s%.1Lf * (1 / %.1Lf) = 1%s\n", KCYN, x, x, KNRM);
    for (int k = 0; k <= n; k++) {
        t_inv(out1, in, k);
        out2[k] = t_prod(in, out1, k);
    }
    output(dp, n, d1, d2, d3, in, out1, out2);

    printf("%s%.1Lf * (%.1Lf / %.1Lf) = %.1Lf%s\n", KCYN, x, y, x, y, KNRM);
    for (int k = 0; k <= n; k++) {
        t_quot(out1, in_B, in, k);
        out2[k] = t_prod(in, out1, k);
    }
    output(dp, n, d1, d2, d3, in, out1, out2);

    printf("%s(SQR(%.1Lf))^0.5 = |%.1Lf|%s\n", KCYN, y, y, KNRM);
    for (int k = 0; k <= n; k++) {
        out1[k] = t_sqr(in_B, k);
        t_pwr(out2, out1, 0.5L, k);
    }
    output(dp, n, d1, d2, d3, in_B, out1, out2);

    printf("%sSQRT(%.1Lf^2) = |%.1Lf|%s\n", KCYN, x, x, KNRM);
    for (int k = 0; k <= n; k++) {
        t_pwr(out1, in, 2.0L, k);
        t_sqrt(out2, out1, k);
    }
    output(dp, n, d1, d2, d3, in, out1, out2);

    printf("%s1 / (%.1Lf^-1) = %.1Lf%s\n", KCYN, x, x, KNRM);
    for (int k = 0; k <= n; k++) {
        t_pwr(out1, in, -1.0L, k);
        t_inv(out2, out1, k);
    }
    output(dp, n, d1, d2, d3, in, out1, out2);

    printf("%sLN(EXP(%.1Lf) = %.1Lf%s\n", KCYN, y, y, KNRM);
    for (int k = 0; k <= n; k++) {
        t_exp(out1, in_B, k);
        t_ln(out2, out1, k);
    }
    output(dp, n, d1, d2, d3, in_B, out1, out2);

    in[0] = MY_PI / y;
    printf("%sSIN_COS(%.3Lf)%s\n", KCYN, in[0], KNRM);
    for (int k = 0; k <= n; k++) {
        t_sin_cos(out1, out2, in, k, TRIG);
        out3[k] = t_sqr(out1, k) + t_sqr(out2, k);
    }
    output(dp, n, d1, d2, d3, out1, out2, out3);

    printf("%sSINH_COSH(%.1Lf)%s\n", KCYN, x, KNRM);
    for (int k = 0; k <= n; k++) {
        t_sin_cos(out1, out2, in, k, HYP);
        out3[k] = - t_sqr(out1, k) + t_sqr(out2, k);
    }
    output(dp, n, d1, d2, d3, out1, out2, out3);

    in[0] = MY_PI / x;
    printf("%sTAN_SEC2(%.3Lf)%s\n", KCYN, in[0], KNRM);
    for (int k = 0; k <= n; k++) {
        t_tan_sec2(out1, out2, in, k, TRIG);
        out3[k] = - t_sqr(out1, k) + out2[k];
    }
    output(dp, n, d1, d2, d3, out1, out2, out3);

    printf("%sTANH_SECH2(%.1Lf)%s\n", KCYN, y, KNRM);
    for (int k = 0; k <= n; k++) {
        t_tan_sec2(out1, out2, in_B, k, HYP);
        out3[k] = t_sqr(out1, k) + out2[k];
    }
    output(dp, n, d1, d2, d3, out1, out2, out3);

    return 0;
}

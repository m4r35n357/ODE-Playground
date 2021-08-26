/*
 *  c99 -g -O3 -Wall  -o test test.c taylor-ode.c -lm
 *
 *  ./taylor-test 6 6 4 -3
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "taylor-ode.h"

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
    sprintf(fs, "[%%2ld] %%+.%ldLe %%+.%ldLe %%+.%ldLe\n", dp, dp, dp);
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
    in_B[1] = fabsl(y) / y;
    series out1 = t_jet(n + 1);
    series out2 = t_jet(n + 1);
    series out3 = t_jet(n + 1);
    series d1 = t_jet(n + 1);
    series d2 = t_jet(n + 1);
    series d3 = t_jet(n + 1);

    printf("x * y\n");
    for (int k = 0; k <= n; k++) {
        out3[k] = t_prod(in, in_B, k);
    }
    output(dp, n, d1, d2, d3, in, in_B, out3);

    printf("x / y\n");
    for (int k = 0; k <= n; k++) {
        t_quot(out3, in, in_B, k);
    }
    output(dp, n, d1, d2, d3, in, in_B, out3);

    printf("|x| * |y| = |x * y|\n");
    for (int k = 0; k <= n; k++) {
        out1[k] = t_abs(in, k);
        out2[k] = t_abs(in_B, k);
        out3[k] = t_prod(out1, out2, k);
    }
    output(dp, n, d1, d2, d3, out1, out2, out3);

    printf("|x| / |y| = |x / y|\n");
    for (int k = 0; k <= n; k++) {
        out1[k] = t_abs(in, k);
        out2[k] = t_abs(in_B, k);
        t_quot(out3, out1, out2, k);
    }
    output(dp, n, d1, d2, d3, out1, out2, out3);

    printf("x * 1/x = 1\n");
    for (int k = 0; k <= n; k++) {
        t_inv(out1, in, k);
        out2[k] = t_prod(in, out1, k);
    }
    output(dp, n, d1, d2, d3, in, out1, out2);

    printf("x * y/x = y\n");
    for (int k = 0; k <= n; k++) {
        t_quot(out1, in_B, in, k);
        out2[k] = t_prod(in, out1, k);
    }
    output(dp, n, d1, d2, d3, in, out1, out2);

    printf("SQRT(SQR(y)) = |y|\n");
    for (int k = 0; k <= n; k++) {
        out1[k] = t_sqr(in_B, k);
        t_sqrt(out2, out1, k);
    }
    output(dp, n, d1, d2, d3, in_B, out1, out2);

    printf("SQRT(x^2) = |x|\n");
    for (int k = 0; k <= n; k++) {
        t_pwr(out1, in, 2.0L, k);
        t_sqrt(out2, out1, k);
    }
    output(dp, n, d1, d2, d3, in, out1, out2);

    printf("1/(x^-1) = x\n");
    for (int k = 0; k <= n; k++) {
        t_pwr(out1, in, -1.0L, k);
        t_inv(out2, out1, k);
    }
    output(dp, n, d1, d2, d3, in, out1, out2);

    printf("LN(EXP(y)) = y\n");
    for (int k = 0; k <= n; k++) {
        t_exp(out1, in_B, k);
        t_ln(out2, out1, k);
    }
    output(dp, n, d1, d2, d3, in_B, out1, out2);

    printf("SIN_COS(pi/y)\n");
    in[0] = MY_PI / y;
    for (int k = 0; k <= n; k++) {
        t_sin_cos(out1, out2, in, k, TRIG);
    }
    output(dp, n, d1, d2, d3, in, out1, out2);

    printf("SINH_COSH(0)\n");
    in[0] = 0.0L;
    for (int k = 0; k <= n; k++) {
        t_sin_cos(out1, out2, in, k, HYP);
    }
    output(dp, n, d1, d2, d3, in, out1, out2);

    printf("TAN_SEC2(pi/x)\n");
    in[0] = MY_PI / x;
    for (int k = 0; k <= n; k++) {
        t_tan_sec2(out1, out2, in, k, TRIG);
    }
    output(dp, n, d1, d2, d3, in, out1, out2);

    printf("TANH_SECH2(0)\n");
    in[0] = 0.0L;
    for (int k = 0; k <= n; k++) {
        t_tan_sec2(out1, out2, in, k, HYP);
    }
    output(dp, n, d1, d2, d3, in, out1, out2);

    return 0;
}

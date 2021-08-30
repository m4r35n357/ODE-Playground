/*
 *  c99 -g -O3 -Wall  -o test test.c taylor-ode.c -lm
 *
 *  ./dual-test 6 _ 4 -3
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "dual.h"
#include "taylor-ode.h"

#define BBLK "\x1B[0;40m"
#define KNRM "\x1B[0;37m"
#define KCYN "\x1B[36m"
#define KGRY "\x1B[2;37m"
#define KBLD "\x1B[1;37m"

static void output (long dp, dual x, dual y, dual z) {
    char fs[128];
    sprintf(fs, "%s%s[%%s]%s%s%s %%+.%ldLe %%+.%ldLe %%+.%ldLe%s\n", KBLD, KGRY, KNRM, BBLK, KGRY, dp, dp, dp, KNRM);
	printf(fs, "val", x.val, y.val, z.val);
	printf(fs, "dot", x.dot, y.dot, z.dot);
}

int main (int argc, char **argv) {
    assert(argc == 5);
    long dp = strtol(argv[1], NULL, 10);
    real x = strtold(argv[3], NULL); assert(x > 0.0L);
    real y = strtold(argv[4], NULL); assert(y != 0.0L);
    dual in = d_var(x);
    dual in_B = d_var(y);
    dual out1, out2, out3;

    printf("%s%.1Lf * %.1Lf = %.1Lf%s\n", KCYN, x, y, x * y, KNRM);
    output(dp, in, in_B, d_mul(in, in_B));

    printf("%s%.1Lf / %.1Lf = %.1Lf%s\n", KCYN, x, y, x / y, KNRM);
    output(dp, in, in_B, d_div(in, in_B));

    printf("%s|%.1Lf| * |%.1Lf| = |%.1Lf * %.1Lf|%s\n", KCYN, x, y, x, y, KNRM);
    out1 = d_abs(in);
    out2 = d_abs(in_B);
    output(dp, out1, out2, d_mul(out1, out2));

    printf("%s|%.1Lf| / |%.1Lf| = |%.1Lf / %.1Lf|%s\n", KCYN, x, y, x, y, KNRM);
    out1 = d_abs(in);
    out2 = d_abs(in_B);
    output(dp, out1, out2, d_div(out1, out2));

    printf("%s%.1Lf * (1 / %.1Lf) = 1%s\n", KCYN, x, x, KNRM);
    out1 = d_inv(in);
    out2 = d_mul(in, out1);
    output(dp, in, out1, out2);

    printf("%s%.1Lf * (%.1Lf / %.1Lf) = %.1Lf%s\n", KCYN, x, y, x, y, KNRM);
    out1 = d_div(in_B, in);
    out2 = d_mul(in, out1);
    output(dp, in, out1, out2);

    printf("%s(SQR(%.1Lf))^0.5 = |%.1Lf|%s\n", KCYN, y, y, KNRM);
    out1 = d_sqr(in_B);
    out2 = d_pow(out1, 0.5L);
    output(dp, in_B, out1, out2);

    printf("%sSQRT(%.1Lf^2) = |%.1Lf|%s\n", KCYN, x, x, KNRM);
    out1 = d_pow(in, 2.0L);
    out2 = d_sqrt(out1);
    output(dp, in, out1, out2);

    printf("%s1 / (%.1Lf^-1) = %.1Lf%s\n", KCYN, x, x, KNRM);
    out1 = d_pow(in, -1.0L);
    out2 = d_inv(out1);
    output(dp, in, out1, out2);

    printf("%sLN(EXP(%.1Lf) = %.1Lf%s\n", KCYN, y, y, KNRM);
    out1 = d_exp(in_B);
    out2 = d_log(out1);
    output(dp, in, out1, out2);

    in = d_var(MY_PI / y);
    printf("%sSIN_COS(%.3Lf)%s\n", KCYN, in.val, KNRM);
    out1 = d_sin(in);
    out2 = d_cos(in);
    out3 = d_add(d_sqr(out1), d_sqr(out2));
    output(dp, out1, out2, out3);

    printf("%sSINH_COSH(%.1Lf)%s\n", KCYN, x, KNRM);
    out1 = d_sinh(in);
    out2 = d_cosh(in);
    out3 = d_sub(d_sqr(out2), d_sqr(out1));
    output(dp, out1, out2, out3);

    in = d_var(MY_PI / x);
    printf("%sTAN_SEC2(%.3Lf)%s\n", KCYN, in.val, KNRM);
    out1 = d_tan(in);
    out2 = d_sqr(d_inv(d_cos(in)));
    out3 = d_sub(out2, d_sqr(out1));
    output(dp, out1, out2, out3);

    printf("%sTANH_SECH2(%.1Lf)%s\n", KCYN, y, KNRM);
    out1 = d_tanh(in_B);
    out2 = d_sqr(d_inv(d_cosh(in_B)));
    out3 = d_add(d_sqr(out1), out2);
    output(dp, out1, out2, out3);

    return 0;
}

/*
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 *
 * Example:  grep -n '[^a-zA-Z]d_[a-z0-9_]*' dual-test.c
 *
 * Example:  c99 --coverage -O0 -o dual-test dual-test.c dual.c -lm
 *
 * Example:  ./dual-test 6 _ 4 -3
 *
 * Example:  gcov -k dual.c
 *
 * Example:  f_list='d_abs d_inv d_sqr d_shift d_scale d_add d_sub d_mul d_div d_exp d_log d_sqrt d_pow d_sin d_cos d_tan d_sinh d_cosh d_tanh'
 *
 * Example:  for f in $f_list; do echo $f; grep -n $f h-*.c; done
 *
 * Example:  for f in $f_list; do echo $f; grep -n $f dual-test.c; done
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "dual.h"

#define BBLK "\x1B[0;40m"
#define KNRM "\x1B[0;37m"
#define KCYN "\x1B[36m"
#define KGRY "\x1B[2;37m"
#define KBLD "\x1B[1;37m"

static void output (long dp, dual x, dual y, dual z) {
    char fs[128];
    sprintf(fs, " %s%s[%%s]%s%s%s %%+.%ldLe %%+.%ldLe %%+.%ldLe%s\n", KBLD, KGRY, KNRM, BBLK, KGRY, dp, dp, dp, KNRM);
    printf(fs, "val", x.val, y.val, z.val);
    printf(fs, "dot", x.dot, y.dot, z.dot);
}

int main (int argc, char **argv) {
    assert(argc == 5);
    long dp = strtol(argv[1], NULL, 10);
    real x = strtold(argv[3], NULL); assert(x > 0.0L);
    real y = strtold(argv[4], NULL); assert(y != 0.0L);
    dual in_X = d_var(x);
    dual in_Y = d_var(y);
    dual out1, out2;

    printf("%s%.1Lf * %.1Lf = %.1Lf%s\n", KCYN, x, y, x * y, KNRM);
    output(dp, in_X, in_Y, d_mul(in_X, in_Y));

    printf("%s%.1Lf / %.1Lf = %.1Lf%s\n", KCYN, x, y, x / y, KNRM);
    output(dp, in_X, in_Y, d_div(in_X, in_Y));

    printf("%s|%.1Lf| * |%.1Lf| = |%.1Lf * %.1Lf|%s\n", KCYN, x, y, x, y, KNRM);
    out1 = d_abs(in_X);
    out2 = d_abs(in_Y);
    output(dp, out1, out2, d_mul(out1, out2));

    printf("%s|%.1Lf| / |%.1Lf| = |%.1Lf / %.1Lf|%s\n", KCYN, x, y, x, y, KNRM);
    out1 = d_abs(in_X);
    out2 = d_abs(in_Y);
    output(dp, out1, out2, d_div(out1, out2));

    printf("%s%.1Lf / %.1Lf = 1%s\n", KCYN, y, y, KNRM);
    output(dp, in_Y, in_Y, d_div(in_Y, in_Y));

    printf("%s%.1Lf * (1 / %.1Lf) = 1%s\n", KCYN, x, x, KNRM);
    out1 = d_inv(in_X);
    out2 = d_mul(in_X, out1);
    output(dp, in_X, out1, out2);

    printf("%s%.1Lf * (%.1Lf / %.1Lf) = %.1Lf%s\n", KCYN, x, y, x, y, KNRM);
    out1 = d_div(in_Y, in_X);
    out2 = d_mul(in_X, out1);
    output(dp, in_X, out1, out2);

    printf("%s(sqr(%.1Lf))^0.5 = |%.1Lf|%s\n", KCYN, y, y, KNRM);
    out1 = d_sqr(in_Y);
    out2 = d_pow(out1, 0.5L);
    output(dp, in_Y, out1, out2);

    printf("%ssqrt(%.1Lf^2) = |%.1Lf|%s\n", KCYN, x, x, KNRM);
    out1 = d_pow(in_X, 2.0L);
    out2 = d_sqrt(out1);
    output(dp, in_X, out1, out2);

    printf("%s1 / (%.1Lf^-1) = %.1Lf%s\n", KCYN, x, x, KNRM);
    out1 = d_pow(in_X, -1.0L);
    out2 = d_inv(out1);
    output(dp, in_X, out1, out2);

    printf("%sexp(ln(%.1Lf) = %.1Lf%s\n", KCYN, x, x, KNRM);
    out1 = d_log(in_X);
    out2 = d_exp(out1);
    output(dp, in_X, out1, out2);

    printf("%sln(exp(%.1Lf) = %.1Lf%s\n", KCYN, y, y, KNRM);
    out1 = d_exp(in_Y);
    out2 = d_log(out1);
    output(dp, in_Y, out1, out2);

    dual in_TRIG = d_var(MY_PI / y);
    printf("%ssin^2(%.3Lf) + cos^2(%.3Lf) = 1.0%s\n", KCYN, in_TRIG.val, in_TRIG.val, KNRM);
    out1 = d_sqr(d_sin(in_TRIG));
    out2 = d_sqr(d_cos(in_TRIG));
    output(dp, out1, out2, d_add(out1, out2));

    printf("%scosh^2(%.3Lf) - sinh^2(%.3Lf) = 1.0%s\n", KCYN, x, x, KNRM);
    out1 = d_sqr(d_sinh(in_X));
    out2 = d_sqr(d_cosh(in_X));
    output(dp, out2, out1, d_sub(out2, out1));

    in_TRIG = d_var(MY_PI / x);
    printf("%ssec^2(%.3Lf) - tanh^2(%.3Lf) = 1.0%s\n", KCYN, in_TRIG.val, in_TRIG.val, KNRM);
    out1 = d_sqr(d_tan(in_TRIG));
    out2 = d_shift(out1, 1.0L);
    output(dp, out2, out1, d_sub(out2, out1));

    printf("%stanh^2(%.1Lf) + sech^2(%.1Lf) = 1.0%s\n", KCYN, y, y, KNRM);
    out1 = d_sqr(d_tanh(in_Y));
    out2 = d_scale(d_shift(out1, -1.0L), -1.0L);
    output(dp, out1, out2, d_add(out1, out2));

    dual out3, out4;
    in_TRIG = d_var(MY_PI / y);
    printf("%ssin(%.3Lf) / cos(%.3Lf) - tan(%.3Lf) = 0.0%s\n", KCYN, in_TRIG.val, in_TRIG.val, in_TRIG.val, KNRM);
    out1 = d_sin(in_TRIG);
    out2 = d_cos(in_TRIG);
    out3 = d_div(out1, out2);
    out4 = d_tan(in_TRIG);
    output(dp, out3, out4, d_sub(out3, out4));

    printf("%ssinh(%.3Lf) / cosh(%.3Lf) - tanh(%.3Lf) = 0.0%s\n", KCYN, x, x, x, KNRM);
    out1 = d_sinh(in_X);
    out2 = d_cosh(in_X);
    out3 = d_div(out1, out2);
    out4 = d_tanh(in_X);
    output(dp, out3, out4, d_sub(out3, out4));

    printf("%s1.0 / cos(%.3Lf) - sqrt(sec^2(%.3Lf)) = 0.0%s\n", KCYN, in_TRIG.val, in_TRIG.val, KNRM);
    out2 = d_cos(in_TRIG);
    out3 = d_inv(out2);
    out4 = d_sqrt(d_shift(d_sqr(d_tan(in_TRIG)), 1.0L));
    output(dp, out3, out4, d_sub(out3, out4));

    printf("%s1.0 / cosh(%.3Lf) - sqrt(sech^2(%.3Lf)) = 0.0%s\n", KCYN, x, x, KNRM);
    out2 = d_cosh(in_X);
    out3 = d_inv(out2);
    out4 = d_sqrt(d_shift(d_scale(d_sqr(d_tanh(in_X)), -1.0L), 1.0L));
    output(dp, out3, out4, d_sub(out3, out4));

    return 0;
}

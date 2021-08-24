/*
 *  c99 -g -O3 -Wall  -o test test.c taylor-ode.c -lm
 *
 *  ./test 6 24
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

real jet_to_derivs (series derivs, series jet, int k) {
	for (int i = 0; i <= k; i++) {
		derivs[i] = jet[i];
	}
    for (int i = 2; i <= k; i++) {
		for (int j = i; j <= k; j++)
			derivs[j] *= i;
    }
    return derivs[k];
}

static void output (long dp, int k, real x, real y, real z) {
    char fs[128];
    sprintf(fs, "[%%2ld] %%+.%ldLe %%+.%ldLe %%+.%ldLe\n", dp, dp, dp);
    printf(fs, k, x, y, z);
}

int main (int argc, char **argv) {
	assert(argc == 3);
    long dp = strtol(argv[1], NULL, 10);
    long n = strtol(argv[2], NULL, 10);
	series in = t_jet(n + 1);
	in[0] = 4.0L;
	in[1] = 1.0L;
	series in_B = t_jet(n + 1);
	in_B[0] = 3.0L;
	in_B[1] = 1.0L;
	series tmp1 = t_jet(n + 1);
	series tmp2 = t_jet(n + 1);
	series out = t_jet(n + 1);
	series derivs = t_jet(n + 1);
	printf("x * 1/x = 1\n");
	for (int k = 0; k <= n; k++) {
		t_inv(tmp1, in, k);
		out[k] = t_prod(in, tmp1, k);
		output(dp, k, in[k], jet_to_derivs(derivs, tmp1, k), jet_to_derivs(derivs, out, k));
	}
	printf("x * y/x = y\n");
	for (int k = 0; k <= n; k++) {
		t_quot(tmp1, in_B, in, k);
		out[k] = t_prod(in, tmp1, k);
		output(dp, k, in[k], jet_to_derivs(derivs, tmp1, k), jet_to_derivs(derivs, out, k));
	}
	printf("SQR(SQRT(x)) = x\n");
	for (int k = 0; k <= n; k++) {
		t_sqrt(tmp1, in, k);
		out[k] = t_sqr(tmp1, k);
		output(dp, k, in[k], jet_to_derivs(derivs, tmp1, k), jet_to_derivs(derivs, out, k));
	}
	printf("(SQRT(x))^2 = x\n");
	for (int k = 0; k <= n; k++) {
		t_sqrt(tmp1, in, k);
		t_pwr(out, tmp1, 2.0L, k);
		output(dp, k, in[k], jet_to_derivs(derivs, tmp1, k), jet_to_derivs(derivs, out, k));
	}
	printf("1/(x^-1) = x\n");
	for (int k = 0; k <= n; k++) {
		tmp1[k] = t_pwr(tmp1, in, -1.0L, k);
		t_inv(out, tmp1, k);
		output(dp, k, in[k], jet_to_derivs(derivs, tmp1, k), jet_to_derivs(derivs, out, k));
	}
	printf("LN(EXP(x)) = x\n");
	in[0] = 0.0L;
	for (int k = 0; k <= n; k++) {
		t_exp(tmp1, in, k);
		t_ln(out, tmp1, k);
		output(dp, k, in[k], jet_to_derivs(derivs, tmp1, k), jet_to_derivs(derivs, out, k));
	}
	printf("SIN_COS\n");
	in[0] = MY_PI / 3.0L;
	for (int k = 0; k <= n; k++) {
		t_sin_cos(tmp1, tmp2, in, k, TRIG);
		output(dp, k, in[k], jet_to_derivs(derivs, tmp1, k), jet_to_derivs(derivs, tmp2, k));
	}
	printf("SINH_COSH\n");
	in[0] = 0.0L;
	for (int k = 0; k <= n; k++) {
		t_sin_cos(tmp1, tmp2, in, k, HYP);
		output(dp, k, in[k], jet_to_derivs(derivs, tmp1, k), jet_to_derivs(derivs, tmp2, k));
	}
	printf("TAN_SEC2\n");
	in[0] = MY_PI / 4.0L;
	for (int k = 0; k <= n; k++) {
		t_tan_sec2(tmp1, tmp2, in, k, TRIG);
		output(dp, k, in[k], jet_to_derivs(derivs, tmp1, k), jet_to_derivs(derivs, tmp2, k));
	}
	printf("TANH_SECH2\n");
	in[0] = 0.0L;
	for (int k = 0; k <= n; k++) {
		t_tan_sec2(tmp1, tmp2, in, k, HYP);
		output(dp, k, in[k], jet_to_derivs(derivs, tmp1, k), jet_to_derivs(derivs, tmp2, k));
	}
	
	return 0;
}

/*
 * (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include <math.h>
#include <mpfr.h>
#include "taylor-ode.h"

const int BASE = 10;

const mpfr_rnd_t RND = GMP_RNDN;

void t_line_output (mpfr_t t, int count, ...) {
    va_list data;
    va_start(data, count);
    for (int i = 0; i < count; i++) {
        mpfr_vprintf("%.9RNe ", data);
    }
    va_end(data);
    mpfr_printf("%.5RNe\n", t);
}

void t_stepper (int argc, char **argv, long *n, mpfr_t *t, mpfr_t *h, long *nsteps) {
    mpfr_set_default_prec(strtod(argv[1], NULL) * log(10.0) / log(2.0));
    *n = strtol(argv[2], NULL, BASE);
    assert(*n > 1);
    mpfr_init_set_ui(*t, 0, RND);
    mpfr_init_set_str(*h, argv[3], BASE, RND);
    *nsteps = strtol(argv[4], NULL, BASE);
}

mpfr_t *t_jet (int n) {
    assert(n > 0);
    mpfr_t *jet = malloc(sizeof (mpfr_t) * n);
    for (int i = 0; i < n; i++) {
        mpfr_init(jet[i]);
    }
    return jet;
}

mpfr_t *t_jet_constant (int n, mpfr_t value) {
    assert(n > 0);
    mpfr_t *jet = t_jet(n);
    mpfr_set(jet[0], value, RND);
    for (int i = 1; i < n; i++) {
         mpfr_set_zero(jet[i], 1);
    }
    return jet;
}

void t_horner (mpfr_t *sum, mpfr_t *jet, long n, mpfr_t h) {
    assert(sizeof *sum == sizeof (mpfr_t));
    assert(n > 0);
    mpfr_set(*sum, jet[n], RND);
    for (int j = n - 1; j > - 1; j--) {
        mpfr_fma(*sum, *sum, h, jet[j], RND);
    }
}

void t_square (mpfr_t *S, mpfr_t *U, int k) {
    assert(sizeof *S == sizeof (mpfr_t));
    assert(k >= 0);
    if (k == 0) {
        mpfr_sqr(*S, U[0], RND);
    } else {
        mpfr_set_zero(*S, 1);
        if (k % 2 == 1) {
            for (int j = 0; j < (k - 1) / 2 + 1; j++) {
                mpfr_fma(*S, U[j], U[k - j], *S, RND);
            }
            mpfr_mul_2ui(*S, *S, 1, RND);
        } else {
            for (int j = 0; j < (k - 2) / 2 + 1; j++) {
                mpfr_fma(*S, U[j], U[k - j], *S, RND);
            }
            mpfr_mul_2ui(*S, *S, 1, RND);
            mpfr_fma(*S, U[k / 2], U[k / 2], *S, RND);
        }
    }
}

void t_product (mpfr_t *P, mpfr_t *U, mpfr_t *V, int k) {
    assert(sizeof *P == sizeof (mpfr_t));
    assert(sizeof *U == sizeof *V);
    assert(k >= 0);
    mpfr_set_zero(*P, 1);
    for (int j = 0; j < k + 1; j++) {
        mpfr_fma(*P, U[j], V[k - j], *P, RND);
    }
}

void t_quotient (mpfr_t *Q, mpfr_t *U, mpfr_t *V, int k) {
    assert(mpfr_sgn(V[0]) != 0);
    assert(Q != U && Q != V && U != V);
    assert(sizeof *U == sizeof *Q && sizeof *V == sizeof *Q);
    assert(k >= 0);
    mpfr_set_zero(Q[k], 1);
    for (int j = 1; j < k + 1; j++) {
        mpfr_fma(Q[k], V[j], Q[k - j], Q[k], RND);
    }
    mpfr_sub(Q[k], U[k], Q[k], RND);
    mpfr_div(Q[k], Q[k], V[0], RND);
}

void t_chain (mpfr_t *dFdX, mpfr_t *dFdU, mpfr_t *U, int k, mpfr_t *dUdX) {
    assert (k > 0);
    mpfr_set_zero(dFdX[k], 1);
    for (int j = 0; j < k; j++) {
        mpfr_mul_ui(*dUdX, U[k - j], k - j, RND);
        mpfr_fma(dFdX[k], dFdU[j], *dUdX, dFdX[k], RND);
    }
    mpfr_div_ui(dFdX[k], dFdX[k], k, RND);
}

void t_exp (mpfr_t *E, mpfr_t *U, int k, mpfr_t *tmp) {
    assert(E != U);
    assert(sizeof *U == sizeof *E);
    assert(sizeof *tmp == sizeof (mpfr_t));
    assert(k >= 0);
    if (k == 0) {
        mpfr_exp(E[0], U[0], RND);
    } else {
        t_chain(E, E, U, k, tmp);
    }
}

void t_sin_cos (mpfr_t *S, mpfr_t *C, mpfr_t *U, int k, mpfr_t *tmp, geometry g) {
    assert(S != C && S != U && C != U);
    assert(sizeof *U == sizeof *S && sizeof *U == sizeof *C);
    assert(sizeof *tmp == sizeof (mpfr_t));
    assert(k >= 0);
    if (k == 0) {
        if (g == TRIG) {
            mpfr_sin_cos(S[0], C[0], U[0], RND);
        } else {
            mpfr_sinh_cosh(S[0], C[0], U[0], RND);
        }
    } else {
        t_chain(S, C, U, k, tmp);
        t_chain(C, S, U, k, tmp);
        if (g == TRIG) {
            mpfr_neg(C[k], C[k], RND);
        }
    }
}

void t_tan_sec2 (mpfr_t *T, mpfr_t *S2, mpfr_t *U, int k, mpfr_t *tmp, geometry g) {
    assert(T != S2 && T != U && S2 != U);
    assert(sizeof *U == sizeof *T && sizeof *U == sizeof *S2);
    assert(sizeof *tmp == sizeof (mpfr_t));
    assert(k >= 0);
    if (k == 0) {
        if (g == TRIG) {
            mpfr_tan(T[0], U[0], RND);
            mpfr_sqr(*tmp, T[0], RND);
            mpfr_add_ui(S2[0], *tmp, 1, RND);
        } else {
            mpfr_tanh(T[0], U[0], RND);
            mpfr_sqr(*tmp, T[0], RND);
            mpfr_ui_sub(S2[0], 1, *tmp, RND);
        }
    } else {
        t_chain(T, S2, U, k, tmp);
        t_chain(S2, T, T, k, tmp);
        mpfr_mul_2ui(S2[k], S2[k], 1, RND);
        if (g == HYP) {
            mpfr_neg(S2[k], S2[k], RND);
        }
    }
}

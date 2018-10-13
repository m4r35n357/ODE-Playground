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
    assert(sum != jet);
    assert(n > 0);
    mpfr_set(*sum, jet[n], RND);
    for (int j = n - 1; j > - 1; j--) {
        mpfr_fma(*sum, *sum, h, jet[j], RND);
    }
}

void t_square (mpfr_t *S, mpfr_t *U, int k) {
    assert(S != U);
    assert(k >= 0);
    if (k == 0) {
        mpfr_sqr(*S, U[0], RND);
    } else {
        mpfr_set_zero(*S, 1);
        if (k % 2 == 1) {
            for (int j = 0; j < (k - 1) / 2 + 1; j++) {
                mpfr_fma(*S, U[j], U[k - j], *S, RND);
            }
            mpfr_mul_si(*S, *S, 2, RND);
        } else {
            for (int j = 0; j < (k - 2) / 2 + 1; j++) {
                mpfr_fma(*S, U[j], U[k - j], *S, RND);
            }
            mpfr_mul_si(*S, *S, 2, RND);
            mpfr_fma(*S, U[k / 2], U[k / 2], *S, RND);
        }
    }
}

void t_product (mpfr_t *P, mpfr_t *V, mpfr_t *U, int k) {
    assert(P != U && P != V);
    assert(k >= 0);
    mpfr_set_zero(*P, 1);
    for (int j = 0; j < k + 1; j++) {
        mpfr_fma(*P, U[j], V[k - j], *P, RND);
    }
}

void t_quotient (mpfr_t *Q, mpfr_t *U, mpfr_t *V, int k) {
    assert(mpfr_sgn(V[0]) != 0);
    assert(Q != U && Q != V && U != V);
    assert(k >= 0);
    mpfr_set_zero(Q[k], 1);
    for (int j = 1; j < k + 1; j++) {
        mpfr_fma(Q[k], V[j], Q[k - j], Q[k], RND);
    }
    mpfr_sub(Q[k], U[k], Q[k], RND);
    mpfr_div(Q[k], Q[k], V[0], RND);
}

void t_sqrt (mpfr_t *R, mpfr_t *U, int k) {
    assert(mpfr_sgn(U[0]) > 0);
    assert(R != U);
    assert(k >= 0);
    if (k == 0) {
        mpfr_sqrt(R[0], U[0], RND);
    } else {
        mpfr_set_zero(R[k], RND);
        if (k % 2 == 1) {
            for (int j = 1; j < (k - 1) / 2 + 1; j++) {
                mpfr_fma(R[k], R[j], R[k - j], R[k], RND);
            }
            mpfr_mul_ui(R[k], R[k], 2, RND);
        } else {
            for (int j = 1; j < (k - 2) / 2 + 1; j++) {
                mpfr_fma(R[k], R[j], R[k - j], R[k], RND);
            }
            mpfr_mul_ui(R[k], R[k], 2, RND);
            mpfr_fma(R[k], R[k / 2], R[k / 2], R[k], RND);
        }
        mpfr_sub(R[k], U[k], R[k], RND);
        mpfr_div_ui(R[k], R[k], 2, RND);
        mpfr_div(R[k], R[k], R[0], RND);
    }
}

void t_power (mpfr_t *P, mpfr_t *U, mpfr_t *a, int k) {
    assert(mpfr_sgn(U[0]) > 0);
    assert(P != U);
    assert(k >= 0);
    mpfr_set_zero(P[k], 1);
    if (k == 0) {
        mpfr_pow(P[0], U[0], *a, RND);
    } else {
        for (int j = 0; j < k; j++) {
            mpfr_mul_si(P[k], *a, k - j, RND);
            mpfr_sub_si(P[k], P[k], j, RND);
            mpfr_mul(P[k], P[k], U[k - j], RND);
            mpfr_mul(P[k], P[k], P[j], RND);
        }
        mpfr_div_si(P[k], P[k], k, RND);
        mpfr_div(P[k], P[k], U[0], RND);
    }
}

void t_exp (mpfr_t *E, mpfr_t *U, int k, mpfr_t *tmp) {
    assert(tmp != E && tmp != U);
    assert(E != U);
    assert(k >= 0);
    if (k == 0) {
        mpfr_exp(E[0], U[0], RND);
    } else {
        mpfr_set_zero(E[k], 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_si(*tmp, U[k - j], k - j, RND);
            mpfr_fma(E[k], *tmp, E[j], E[k], RND);
        }
        mpfr_div_si(E[k], E[k], k, RND);
    }
}

void t_sin_cos (mpfr_t *S, mpfr_t *C, mpfr_t *U, int k, mpfr_t *tmp) {
    assert(tmp != S && tmp != C && tmp != U);
    assert(S != C && S != U && C != U);
    assert(k >= 0);
    if (k == 0) {
        mpfr_sin_cos(S[0], C[0], U[0], RND);
    } else {
        mpfr_set_zero(S[k], 1);
        mpfr_set_zero(C[k], 1);
        for (int j = 0; j < k; j++) {
            mpfr_mul_si(*tmp, U[k - j], k - j, RND);
            mpfr_fma(S[k], *tmp, C[j], S[k], RND);
            mpfr_fma(C[k], *tmp, S[j], C[k], RND);
        }
        mpfr_div_si(S[k], S[k], k, RND);
        mpfr_div_si(C[k], C[k], - k, RND);
    }
}

/*
 * (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"
#include "ad.h"

void jet_output (mpfr_t *jet, long n, char* f_colour, char *fk_colour) {
    mpfr_printf("%s%9.6RNf ", f_colour, jet[0]);
    for (int i = 1; i < n; i++) {
        mpfr_printf("%s%9.6RNf ", fk_colour, jet[i]);
    }
    printf("%s\n", f_colour);
}

void jet_to_derivs (mpfr_t *jet, long n) {
    for (int i = 1, fac = 1; i < n; i++) {
        fac *= i;
        mpfr_mul_si(jet[i], jet[i], fac, RND);
    }
}

void derivative_output (mpfr_t *jet, long n, char* f_colour, char *fk_colour) {
    jet_to_derivs(jet, n);
    jet_output(jet, n, f_colour, fk_colour);
}

int ad_bisect (model m, mpfr_t *xa, mpfr_t *xb, int max_it, mpfr_t f_tol, mpfr_t x_tol) {
    int counter = 0;
    mpfr_t delta, tmp, *fa, *xc, *fc;
    mpfr_inits(tmp, delta, NULL);
    xc = t_jet(1);
    fa = t_jet(1);
    fc = t_jet(1);
    m(fa, xa, 1);
    mpfr_sub(delta, xb[0], xa[0], RND);
    while(mpfr_cmp_abs(fc[0], f_tol) >= 0 || mpfr_cmp_abs(delta, x_tol) >= 0) {
        mpfr_add(xc[0], xa[0], xb[0], RND);
        mpfr_div_ui(xc[0], xc[0], 2, RND);
        m(fc, xc, 1);
        mpfr_mul(tmp, fa[0], fc[0], RND);
        if (mpfr_sgn(tmp) < 0) {
            mpfr_set(xb[0], xc[0], RND);
        } else {
            mpfr_set(xa[0], xc[0], RND);
        }
        mpfr_sub(delta, xb[0], xa[0], RND);
        counter++;
        if (counter > max_it) return 1;
    }
    mpfr_fprintf(stderr, "%3d %20.12RNe %20.12RNe %20.12RNe\n", counter, xc[0], fc[0], delta);
    return 0;
}

int ad_newton (model m, mpfr_t *f, mpfr_t *x, int max_it, mpfr_t f_tol, mpfr_t x_tol) {
    int counter = 0;
    mpfr_t delta;
    mpfr_init_set_ui(delta, 1, RND);
    while(mpfr_cmp_abs(f[0], f_tol) >= 0 || mpfr_cmp_abs(delta, x_tol) >= 0) {
        mpfr_set_si(x[1], 1, RND);
        m(f, x, 2);
        mpfr_div(delta, f[0], f[1], RND);
        mpfr_sub(x[0], x[0], delta, RND);
        counter++;
        if (counter > max_it) return 1;
    }
    mpfr_fprintf(stderr, "%3d %20.12RNe %20.12RNe %20.12RNe\n", counter, x[0], f[0], delta);
    return 0;
}

int ad_householder (model m, mpfr_t *f, mpfr_t *x, long n, int max_it, mpfr_t f_tol, mpfr_t x_tol) {
    int counter = 0;
    mpfr_t delta, D1, *f_reciprocal, *w1;
    mpfr_init_set_ui(delta, 1, RND);
    f_reciprocal = t_jet(n);
    mpfr_init_set_si(D1, 1, RND);
    w1 = t_jet_constant(n, D1);
    while(mpfr_cmp_abs(f[0], f_tol) >= 0 || mpfr_cmp_abs(delta, x_tol) >= 0) {
        mpfr_set_si(x[1], 1, RND);
        m(f, x, n);
        ad_quotient(f_reciprocal, w1, f, n);
        jet_to_derivs(f_reciprocal, n);
        mpfr_div(delta, f_reciprocal[n - 2], f_reciprocal[n - 1], RND);
        mpfr_mul_si(delta, delta, n - 1, RND);
        mpfr_add(x[0], x[0], delta, RND);
        counter++;
        if (counter > max_it) return 1;
    }
    mpfr_fprintf(stderr, "%3d %20.12RNe %20.12RNe %20.12RNe\n", counter, x[0], f[0], delta);
    return 0;
}

void ad_scale (mpfr_t *S, mpfr_t *U, mpfr_t a, int n) {
    for (int k = 0; k < n; k++) {
        mpfr_mul(S[k], U[k], a, RND);
    }
}

void ad_plus (mpfr_t *P, mpfr_t *V, mpfr_t *U, int n) {
    for (int k = 0; k < n; k++) {
        mpfr_add(P[k], V[k], U[k], RND);
    }
}

void ad_minus (mpfr_t *P, mpfr_t *V, mpfr_t *U, int n) {
    for (int k = 0; k < n; k++) {
        mpfr_sub(P[k], V[k], U[k], RND);
    }
}

void ad_square (mpfr_t *S, mpfr_t *U, int n, mpfr_t *tmp) {
    assert(S != U);
    for (int k = 0; k < n; k++) {
        t_square(tmp, U, k);
        mpfr_set(S[k], *tmp, RND);
    }
}

void ad_product (mpfr_t *P, mpfr_t *V, mpfr_t *U, int n, mpfr_t *tmp) {
    assert(P != U && P != V);
    for (int k = 0; k < n; k++) {
        t_product(tmp, V, U, k);
        mpfr_set(P[k], *tmp, RND);
    }
}

void ad_quotient (mpfr_t *Q, mpfr_t *U, mpfr_t *V, int n) {
    for (int k = 0; k < n; k++) {
        t_quotient(Q, U, V, k);
    }
}

void ad_sqrt (mpfr_t *R, mpfr_t *U, int n) {
    for (int k = 0; k < n; k++) {
        t_sqrt(R, U, k);
    }
}

void ad_power (mpfr_t *P, mpfr_t *U, mpfr_t *a, int n) {
    for (int k = 0; k < n; k++) {
        t_power(P, U, a, k);
    }
}

void ad_exp (mpfr_t *E, mpfr_t *U, int n, mpfr_t *tmp) {
    for (int k = 0; k < n; k++) {
        t_exp(E, U, k, tmp);
    }
}

void ad_sin_cos (mpfr_t *S, mpfr_t *C, mpfr_t *U, int n, mpfr_t *tmp) {
    for (int k = 0; k < n; k++) {
        t_sin_cos(S, C, U, k, tmp);
    }
}

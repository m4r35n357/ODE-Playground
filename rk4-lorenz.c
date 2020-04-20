/*
 * Lorenz System using RK4
 *
 * Example: ./rk4-lorenz-dbg 32 1 .01 10000 -15.8 -17.48 35.64 10 28 8 3
 *
 * NOTE plot interval in place of order!
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

static void dx (mpfr_t *lhs, mpfr_t x, mpfr_t y, mpfr_t s) {
    mpfr_fmms(*lhs, s, y, s, x, RND);
}

static void dy (mpfr_t *lhs, mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t r) {
    mpfr_sub(*lhs, r, z, RND);
    mpfr_fms(*lhs, *lhs, x, y, RND);
}

static void dz (mpfr_t *lhs, mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t b) {
    mpfr_fmms(*lhs, x, y, b, z, RND);
}

static mpfr_t *sum (mpfr_t *b, mpfr_t a1, mpfr_t a2, mpfr_t a3, mpfr_t a4) {
    mpfr_add(*b, a1, a4, RND);
    mpfr_mul_2ui(a2, a2, 1, RND);
    mpfr_add(*b, *b, a2, RND);
    mpfr_mul_2ui(a3, a3, 1, RND);
    mpfr_add(*b, *b, a3, RND);
    mpfr_div_ui(*b, *b, 6, RND);
    return b;
}

int main (int argc, char **argv) {
    long nsteps, interval;
    mpfr_t x, y, z, sigma, rho, beta, h, _, h_2, _x, _y, _z;
    mpfr_t k1, l1, m1, k2, l2, m2, k3, l3, m3, k4, l4, m4;

    assert(argc == 12);
    t_stepper(argv, &interval, &h, &nsteps);
    t_args(argv, argc, &x, &y, &z, &sigma, &rho, &beta, &_);
    mpfr_div(beta, beta, _, RND);
    mpfr_inits(h_2, _x, _y, _z, k1, l1, m1, k2, l2, m2, k3, l3, m3, k4, l4, m4, NULL);

    mpfr_div_2ui(h_2, h, 1, RND);

    t_output(x, y, z, h, 0, _);
    for (long step = 1; step <= nsteps; step++) {
        dx(&k1, x, y, sigma);
        dy(&l1, x, y, z, rho);
        dz(&m1, x, y, z, beta);

        mpfr_fma(_x, k1, h_2, x, RND);
        mpfr_fma(_y, l1, h_2, y, RND);
        mpfr_fma(_z, m1, h_2, z, RND);
        dx(&k2, _x, _y, sigma);
        dy(&l2, _x, _y, _z, rho);
        dz(&m2, _x, _y, _z, beta);

        mpfr_fma(_x, k2, h_2, x, RND);
        mpfr_fma(_y, l2, h_2, y, RND);
        mpfr_fma(_z, m2, h_2, z, RND);
        dx(&k3, _x, _y, sigma);
        dy(&l3, _x, _y, _z, rho);
        dz(&m3, _x, _y, _z, beta);

        mpfr_fma(_x, k3, h, x, RND);
        mpfr_fma(_y, l3, h, y, RND);
        mpfr_fma(_z, m3, h, z, RND);
        dx(&k4, _x, _y, sigma);
        dy(&l4, _x, _y, _z, rho);
        dz(&m4, _x, _y, _z, beta);

        mpfr_fma(x, h, *sum(&_, k1, k2, k3, k4), x, RND);
        mpfr_fma(y, h, *sum(&_, l1, l2, l3, l4), y, RND);
        mpfr_fma(z, h, *sum(&_, m1, m2, m3, m4), z, RND);

        if (step % interval == 0) { t_output(x, y, z, h, step, _); }
    }
    return 0;
}

/*
 * Lorenz System
 *
 * Example: ./rk4-lorenz-dbg 16 16 .01 10000 -15.8 -17.48 35.64 10 28 8 3
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

static void dx (mpfr_t *lhs, mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t s) {
    mpfr_fmms(*lhs, s, y, s, x, RND);
}

static void dy (mpfr_t *lhs, mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t r) {
    mpfr_sub(*lhs, r, z, RND);
    mpfr_fms(*lhs, *lhs, x, y, RND);
}

static void dz (mpfr_t *lhs, mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t b) {
    mpfr_fmms(*lhs, x, y, b, z, RND);
}

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, y0, z0, sigma, rho, beta, h, _, h_2, _x, _y, _z;
    mpfr_t k1, l1, m1, k2, l2, m2, k3, l3, m3, k4, l4, m4;

    assert(argc == 12);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &sigma, &rho, &beta, &_);
    mpfr_div(beta, beta, _, RND);
    mpfr_inits(h_2, _x, _y, _z, k1, l1, m1, k2, l2, m2, k3, l3, m3, k4, l4, m4, NULL);

    mpfr_div_2ui(h_2, h, 1, RND);

    t_xyz_output(x0, y0, z0, t);
    for (long step = 1; step < nsteps + 1; step++) {
        dx(&k1, x0, y0, z0, sigma);
        dy(&l1, x0, y0, z0, rho);
        dz(&m1, x0, y0, z0, beta);

        mpfr_fma(_x, k1, h_2, x0, RND);
        mpfr_fma(_y, l1, h_2, y0, RND);
        mpfr_fma(_z, m1, h_2, z0, RND);
        dx(&k2, _x, _y, _z, sigma);
        dy(&l2, _x, _y, _z, rho);
        dz(&m2, _x, _y, _z, beta);

        mpfr_fma(_x, k2, h_2, x0, RND);
        mpfr_fma(_y, l2, h_2, y0, RND);
        mpfr_fma(_z, m2, h_2, z0, RND);
        dx(&k3, _x, _y, _z, sigma);
        dy(&l3, _x, _y, _z, rho);
        dz(&m3, _x, _y, _z, beta);

        mpfr_fma(_x, k3, h, x0, RND);
        mpfr_fma(_y, l3, h, y0, RND);
        mpfr_fma(_z, m3, h, z0, RND);
        dx(&k4, _x, _y, _z, sigma);
        dy(&l4, _x, _y, _z, rho);
        dz(&m4, _x, _y, _z, beta);

        mpfr_add(_, k1, k4, RND);
        mpfr_mul_2ui(k2, k2, 1, RND);
        mpfr_add(_, _, k2, RND);
        mpfr_mul_2ui(k3, k3, 1, RND);
        mpfr_add(_, _, k3, RND);
        mpfr_div_ui(_, _, 6, RND);
        mpfr_fma(x0, h, _, x0, RND);

        mpfr_add(_, l1, l4, RND);
        mpfr_mul_2ui(l2, l2, 1, RND);
        mpfr_add(_, _, l2, RND);
        mpfr_mul_2ui(l3, l3, 1, RND);
        mpfr_add(_, _, l3, RND);
        mpfr_div_ui(_, _, 6, RND);
        mpfr_fma(y0, h, _, y0, RND);

        mpfr_add(_, m1, m4, RND);
        mpfr_mul_2ui(m2, m2, 1, RND);
        mpfr_add(_, _, m2, RND);
        mpfr_mul_2ui(m3, m3, 1, RND);
        mpfr_add(_, _, m3, RND);
        mpfr_div_ui(_, _, 6, RND);
        mpfr_fma(z0, h, _, z0, RND);

        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x0, y0, z0, t);
    }
    return 0;
}

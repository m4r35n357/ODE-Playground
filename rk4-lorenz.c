/*
 * Lorenz System using RK4
 *
 * Example: ./rk4-lorenz-dbg 32 16 .00001 10000000 -15.8 -17.48 35.64 10 28 8 3 | sed -n '1~1000p' >/tmp/dataB
 *          ./divergence.py d-604-502-1400 /tmp/dataB 3 1.0 >/dev/null
 *          ./compare.py d-604-502-1400 /tmp/dataB 3 >/dev/null &
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

        mpfr_fma(x0, h, *sum(&_, k1, k2, k3, k4), x0, RND);
        mpfr_fma(y0, h, *sum(&_, l1, l2, l3, l4), y0, RND);
        mpfr_fma(z0, h, *sum(&_, m1, m2, m3, m4), z0, RND);

        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x0, y0, z0, t);
    }
    return 0;
}

/*
 * Lorenz System
 *
 * Example: ./rk4-lorenz-dbg NA NA 1 .01 10000 -15.8 -17.48 35.64 10 28 8 3
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

static long double dx (long double x, long double y, long double s) {
    return s * (y - x);
}

static long double dy (long double x, long double y, long double z, long double r) {
    return r * x - y - x * z;
}

static long double dz (long double x, long double y, long double z, long double b) {
    return x * y - b * z;
}

int main (int argc, char **argv) {
    long n, nsteps;
    long double x0, y0, z0, sigma, rho, beta, h, _;

    assert(argc == 13);
    t_stepper(argv, &n, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &sigma, &rho, &beta, &_);
    beta /= _;

    long double k1, l1, m1, k2, l2, m2, k3, l3, m3, k4, l4, m4;

    t_xyz_output(x0, y0, z0, 0.0);
    for (long step = 1; step < nsteps + 1; step++) {
        k1 = dx(x0,                y0,                                   sigma);
        l1 = dy(x0,                y0,                z0,                rho);
        m1 = dz(x0,                y0,                z0,                beta);
        k2 = dx(x0 + 0.5 * k1 * h, y0 + 0.5 * l1 * h,                    sigma);
        l2 = dy(x0 + 0.5 * k1 * h, y0 + 0.5 * l1 * h, z0 + 0.5 * m1 * h, rho);
        m2 = dz(x0 + 0.5 * k1 * h, y0 + 0.5 * l1 * h, z0 + 0.5 * m1 * h, beta);
        k3 = dx(x0 + 0.5 * k2 * h, y0 + 0.5 * l2 * h,                    sigma);
        l3 = dy(x0 + 0.5 * k2 * h, y0 + 0.5 * l2 * h, z0 + 0.5 * m2 * h, rho);
        m3 = dz(x0 + 0.5 * k2 * h, y0 + 0.5 * l2 * h, z0 + 0.5 * m2 * h, beta);
        k4 = dx(x0 +       k3 * h, y0 +       l3 * h,                    sigma);
        l4 = dy(x0 +       k3 * h, y0 +       l3 * h, z0 +       m3 * h, rho);
        m4 = dz(x0 +       k3 * h, y0 +       l3 * h, z0 +       m3 * h, beta);

        x0 += h * (k1 + 2.0 * (k2 + k3) + k4) / 6.0;
        y0 += h * (l1 + 2.0 * (l2 + l3) + l4) / 6.0;
        z0 += h * (m1 + 2.0 * (m2 + m3) + m4) / 6.0;
        t_xyz_output(x0, y0, z0, h * step);
    }
    return 0;
}

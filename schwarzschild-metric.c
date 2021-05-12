/*
 * Schwarzschild metric
 *
 *  Example: ./schwarzschild-metric-dbg 15 10 0.01 10000  0.0 12.0 1.57079632679 0.0  1.0 0.0 0.0 0.0  1.0 0.8
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <math.h>
#include "gr.h"

dual g_t_t (real m, real a, dual r, dual theta) {
    (void)a; (void)theta;
    return d_shift(d_scale(d_inv(r), 2.0L * m), - 1.0L);
}

dual g_t_r (real m, real a, dual r, dual theta) {
    (void)m; (void)a; (void)r; (void)theta;
    return d_dual(0.0L);
}

dual g_t_theta (real m, real a, dual r, dual theta) {
    (void)m; (void)a; (void)r; (void)theta;
    return d_dual(0.0L);
}

dual g_t_phi (real m, real a, dual r, dual theta) {
    (void)m; (void)a; (void)r; (void)theta;
    return d_dual(0.0L);
}

dual g_r_r (real m, real a, dual r, dual theta) {
    (void)a; (void)theta;
    return d_inv(d_shift(d_scale(d_inv(r), - 2.0L * m), 1.0L));
}

dual g_r_theta (real m, real a, dual r, dual theta) {
    (void)m; (void)a; (void)r; (void)theta;
    return d_dual(0.0L);
}

dual g_r_phi (real m, real a, dual r, dual theta) {
    (void)m; (void)a; (void)r; (void)theta;
    return d_dual(0.0L);
}

dual g_theta_theta (real m, real a, dual r, dual theta) {
    (void)m; (void)a; (void)theta;
    return d_sqr(r);
}

dual g_theta_phi (real m, real a, dual r, dual theta) {
    (void)m; (void)a; (void)r; (void)theta;
    return d_dual(0.0L);
}

dual g_phi_phi (real m, real a, dual r, dual theta) {
    (void)m; (void)a;
    return d_mul(d_sqr(r), d_sqr(d_sin(theta)));
}

real i_t_t (real m, real a, real r, real theta) {
    (void)a; (void)theta;
    return - r / (r - 2.0L * m);
}

real i_t_r (real m, real a, real r, real theta) {
    (void)m; (void)a; (void)r; (void)theta;
    return 0.0L;
}

real i_t_theta (real m, real a, real r, real theta) {
    (void)m; (void)a; (void)r; (void)theta;
    return 0.0L;
}

real i_t_phi (real m, real a, real r, real theta) {
    (void)m; (void)a; (void)r; (void)theta;
    return 0.0L;
}

real i_r_r (real m, real a, real r, real theta) {
    (void)a; (void)theta;
    return (r - 2.0L * m) / r;
}

real i_r_theta (real m, real a, real r, real theta) {
    (void)m; (void)a; (void)r; (void)theta;
    return 0.0L;
}

real i_r_phi (real m, real a, real r, real theta) {
    (void)m; (void)a; (void)r; (void)theta;
    return 0.0L;
}

real i_theta_theta (real m, real a, real r, real theta) {
    (void)m; (void)a; (void)theta;
    return 1.0L / (r * r);
}

real i_theta_phi (real m, real a, real r, real theta) {
    (void)m; (void)a; (void)r; (void)theta;
    return 0.0L;
}

real i_phi_phi (real m, real a, real r, real theta) {
    (void)m; (void)a;
    return 1.0L / (r * r *sinl(theta) * sinl(theta));
}

int main (int argc, char **argv) {
    assert(argc == 15);
    tsm4(argc, argv);
    return 0;
}

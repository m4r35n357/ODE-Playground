/*
 * Boyer-Lindquist metric
 *
 *  Compile: c99 -g -O0 -o boyer-lindquist-metric boyer-lindquist-metric.c gr.c dual.c taylor-ode.c -lm
 *
 *  Execute: ./boyer-lindquist-metric 15 10 0.01 10000  0.0 12.0 1.57079632679 0.0  1.0 0.0 0.0 0.0  1.0 0.8
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <math.h>
#include "gr.h"

dual g_t_t (real m, real a, dual r, dual theta) {
    return d_shift(d_div(d_scale(r, 2.0L * m), d_sigma(a, r, theta)), - 1.0L);
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
    return d_scale(d_div(d_mul(d_scale(r, 2.0L * m), d_sqr(d_sin(theta))), d_sigma(a, r, theta)), -a);
}

dual g_r_r (real m, real a, dual r, dual theta) {
    return d_div(d_sigma(a, r, theta), d_delta(m, a, r));
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
    (void)m;
    return d_sigma(a, r, theta);
}

dual g_theta_phi (real m, real a, dual r, dual theta) {
    (void)m; (void)a; (void)r; (void)theta;
    return d_dual(0.0L);
}

dual g_phi_phi (real m, real a, dual r, dual theta) {
    return d_mul(d_sqr(d_sin(theta)),
           d_add(d_div(d_mul(d_scale(r, 2.0L * m * a * a), d_sqr(d_sin(theta))), d_sigma(a, r, theta)), d_ra2(a, r)));
}

real i_t_t (real m, real a, real r, real theta) {
    return (a * a * r_delta(m, a, r) * sinl(theta) * sinl(theta) - r_ra2(a, r) * r_ra2(a, r)) /
           (r_delta(m, a, r) * r_sigma(a, r, theta));
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
    return - 2.0L * m * r * a / (r_delta(m, a, r) * r_sigma(a, r, theta));
}

real i_r_r (real m, real a, real r, real theta) {
    return r_delta(m, a, r) / r_sigma(a, r, theta);
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
    (void)m;
    return 1.0L / r_sigma(a, r, theta);
}

real i_theta_phi (real m, real a, real r, real theta) {
    (void)m; (void)a; (void)r; (void)theta;
    return 0.0L;
}

real i_phi_phi (real m, real a, real r, real theta) {
    return (r_sigma(a, r, theta) - 2.0L * m * r) / (r_delta(m, a, r) * r_sigma(a, r, theta) * sinl(theta) * sinl(theta));
}

int main (int argc, char **argv) {
    assert(argc == 15);
    tsm4(argc, argv);
    return 0;
}

/*
 *  Compile: c99 -g -O0 -o schwarzschild-metric schwarzschild-metric.c gr.c dual.c taylor-ode.c -lm
 *
 *  Execute: ./schwarzschild-metric 15 10 0.01 10000  0.0 12.0 1.57079632679 0.0  1.0 0.0 0.0 0.0  1.0 0.8
 */

#include <stdlib.h>
#include <math.h>
#include "gr.h"

dual g_t_t (real m, real a, dual r, dual theta) {
    return d_shift(d_scale(d_inv(r), 2.0L * m), - 1.0L);
}

dual g_t_r (real m, real a, dual r, dual theta) {
    return d_dual(0.0L);
}

dual g_t_theta (real m, real a, dual r, dual theta) {
    return d_dual(0.0L);
}

dual g_t_phi (real m, real a, dual r, dual theta) {
    return d_dual(0.0L);
}

dual g_r_r (real m, real a, dual r, dual theta) {
    return d_inv(d_shift(d_neg(d_scale(d_inv(r), 2.0L * m)), 1.0L));
}

dual g_r_theta (real m, real a, dual r, dual theta) {
    return d_dual(0.0L);
}

dual g_r_phi (real m, real a, dual r, dual theta) {
    return d_dual(0.0L);
}

dual g_theta_theta (real m, real a, dual r, dual theta) {
    return d_sqr(r);
}

dual g_theta_phi (real m, real a, dual r, dual theta) {
    return d_dual(0.0L);
}

dual g_phi_phi (real m, real a, dual r, dual theta) {
    return d_mul(d_sqr(r), d_sqr(d_sin(theta)));
}

real i_t_t (real m, real a, real r, real theta) {
    return - r / (r - 2.0L * m);
}

real i_t_r (real m, real a, real r, real theta) {
    return 0.0L;
}

real i_t_theta (real m, real a, real r, real theta) {
    return 0.0L;
}

real i_t_phi (real m, real a, real r, real theta) {
    return 0.0L;
}

real i_r_r (real m, real a, real r, real theta) {
    return (r - 2.0L * m) / r;
}

real i_r_theta (real m, real a, real r, real theta) {
    return 0.0L;
}

real i_r_phi (real m, real a, real r, real theta) {
    return 0.0L;
}

real i_theta_theta (real m, real a, real r, real theta) {
    return 1.0L / (r * r);
}

real i_theta_phi (real m, real a, real r, real theta) {
    return 0.0L;
}

real i_phi_phi (real m, real a, real r, real theta) {
    return 1.0L / (r * r *sinl(theta) * sinl(theta));
}

int main (int argc, char **argv) {
    assert(argc == 15);
    tsm4(argc, argv);

    real m = 1.0L;
    real a = 0.8L;
    real r = 12.0L;
    real theta = 0.5L * get_PI();

    matrix4x4 metric, inverse, result;
    real_metric (metric, m, a, r, theta);
    real_inverse (inverse, m, a, r, theta);
    mm_mult(result, metric, inverse);

    return 0;
}

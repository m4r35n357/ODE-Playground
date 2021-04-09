/*
 * Doran metric
 *
 *  Compile: c99 -g -O0 -o doran-metric doran-metric.c gr.c dual.c taylor-ode.c -lm
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
    return d_sqrt(d_div(d_scale(r, 2.0L * m), d_ra2(a, r)));
}

dual g_t_theta (real m, real a, dual r, dual theta) {
    return d_dual(0.0L);
}

dual g_t_phi (real m, real a, dual r, dual theta) {
    return d_div(d_scale(d_mul(d_scale(r, 2.0L * m), d_sqr(d_sin(theta))), - a), d_sigma(a, r, theta));
}

dual g_r_r (real m, real a, dual r, dual theta) {
    return d_div(d_sigma(a, r, theta), d_ra2(a, r));
}

dual g_r_theta (real m, real a, dual r, dual theta) {
    return d_dual(0.0L);
}

dual g_r_phi (real m, real a, dual r, dual theta) {
    return d_scale(d_mul(d_sqrt(d_div(d_scale(r, 2.0L * m), d_ra2(a, r))), d_sqr(d_sin(theta))), - a);
}

dual g_theta_theta (real m, real a, dual r, dual theta) {
    return d_sigma(a, r, theta);
}

dual g_theta_phi (real m, real a, dual r, dual theta) {
    return d_dual(0.0L);
}

dual g_phi_phi (real m, real a, dual r, dual theta) {
    return d_add(d_div(d_mul(d_scale(r, 2.0L * m * a * a), d_sqr(d_sqr(d_sin(theta)))), d_sigma(a, r, theta)),
           d_mul(d_ra2(a, r), d_sqr(d_sin(theta))));
}

real i_t_t (real m, real a, real r, real theta) {
    return - 1.0L;
}

real i_t_r (real m, real a, real r, real theta) {
    return sqrtl(2.0L * m * r * r_ra2(a, r)) / r_sigma(a, r, theta);
}

real i_t_theta (real m, real a, real r, real theta) {
    return 0.0L;
}

real i_t_phi (real m, real a, real r, real theta) {
    return 0.0L;
}

real i_r_r (real m, real a, real r, real theta) {
    return r_delta(m, a, r) / r_sigma(a, r, theta);
}

real i_r_theta (real m, real a, real r, real theta) {
    return 0.0L;
}

real i_r_phi (real m, real a, real r, real theta) {
    return a * sqrtl(2.0L * m * r / r_ra2(a, r)) / r_sigma(a, r, theta);
}

real i_theta_theta (real m, real a, real r, real theta) {
    return 1.0L / r_sigma(a, r, theta);
}

real i_theta_phi (real m, real a, real r, real theta) {
    return 0.0L;
}

real i_phi_phi (real m, real a, real r, real theta) {
    return 1.0L / (r_ra2(a, r) * sinl(theta) * sinl(theta));
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

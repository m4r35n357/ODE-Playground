/*
 *  Compile: c99 -g -O0 -o doran-metric doran-metric.c gr.c dual.c -lm
 */

#include <stdlib.h>
#include <math.h>
#include "gr.h"

dual g_t_t (real m, real a, dual r, dual theta) {
    return d_shift(d_div(d_scale(r, 2.0L * m), sigma(a, r, theta)), - 1.0L);
}

dual g_t_r (real m, real a, dual r, dual theta) {
    return d_sqrt(d_div(d_scale(r, 2.0L * m), ra2(a, r)));
}

dual g_t_theta (real m, real a, dual r, dual theta) {
    return d_dual(0.0L);
}

dual g_t_phi (real m, real a, dual r, dual theta) {
    return d_div(d_scale(d_mul(d_scale(r, 2.0L * m), d_sqr(d_sin(theta))), - a), sigma(a, r, theta));
}

dual g_r_r (real m, real a, dual r, dual theta) {
    return d_div(sigma(a, r, theta), ra2(a, r));
}

dual g_r_theta (real m, real a, dual r, dual theta) {
    return d_dual(0.0L);
}

dual g_r_phi (real m, real a, dual r, dual theta) {
    return d_scale(d_mul(d_sqrt(d_div(d_scale(r, 2.0L * m), ra2(a, r))), d_sqr(d_sin(theta))), - a);
}

dual g_theta_theta (real m, real a, dual r, dual theta) {
    return sigma(a, r, theta);
}

dual g_theta_phi (real m, real a, dual r, dual theta) {
    return d_dual(0.0L);
}

dual g_phi_phi (real m, real a, dual r, dual theta) {
    return d_add(d_div(d_mul(d_scale(r, 2.0L * m * a * a), d_sqr(d_sqr(d_sin(theta)))), sigma(a, r, theta)),
           d_mul(ra2(a, r), d_sqr(d_sin(theta))));
}

real i_t_t (real m, real a, real r, real theta) {
    return - 1.0L;
}

real i_t_r (real m, real a, real r, real theta) {
    return sqrtl(d_scale(d_dual(r), 2.0L * m).val * ra2(a, d_dual(r)).val) / sigma(a, d_dual(r), d_dual(theta)).val;
}

real i_t_theta (real m, real a, real r, real theta) {
    return 0.0L;
}

real i_t_phi (real m, real a, real r, real theta) {
    return 0.0L;
}

real i_r_r (real m, real a, real r, real theta) {
    return delta(m, a, d_dual(r)).val / sigma(a, d_dual(r), d_dual(theta)).val;
}

real i_r_theta (real m, real a, real r, real theta) {
    return 0.0L;
}

real i_r_phi (real m, real a, real r, real theta) {
    return a * sqrtl(d_scale(d_dual(r), 2.0L * m).val / ra2(a, d_dual(r)).val) / sigma(a, d_dual(r), d_dual(theta)).val;
}

real i_theta_theta (real m, real a, real r, real theta) {
    return 1.0L / sigma(a, d_dual(r), d_dual(theta)).val;
}

real i_theta_phi (real m, real a, real r, real theta) {
    return 0.0L;
}

real i_phi_phi (real m, real a, real r, real theta) {
    return 1.0L / (ra2(a, d_dual(r)).val * sinl(theta) * sinl(theta));
}

int main (int argc, char **argv) {
    (void) argc;
    (void) argv;
    real m = 1.0L;
    real a = 0.8L;
    real r = 12.0L;
    real theta = 0.5L;
    matrix4x4 metric, inverse, result;
    real_metric (metric, m, a, r, theta);
    real_inverse (inverse, m, a, r, theta);
    mm_mult(result, metric, inverse);
    return 0;
}

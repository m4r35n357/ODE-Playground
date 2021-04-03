/*
 *  Compile: c99 -g -O0 -o minkowski-metric minkowski-metric.c gr.c -lm
 */

#include <stdlib.h>
#include <math.h>
#include "gr.h"

 real g_t_t (real m, real a, real r, real theta) {
    return 1.0L;
}

 real g_t_r (real m, real a, real r, real theta) {
    return 0.0L;
}

 real g_t_theta (real m, real a, real r, real theta) {
    return 0.0L;
}

 real g_t_phi (real m, real a, real r, real theta) {
    return 0.0L;
}

 real g_r_r (real m, real a, real r, real theta) {
    return 1.0L;
}

 real g_r_theta (real m, real a, real r, real theta) {
    return 0.0L;
}

 real g_r_phi (real m, real a, real r, real theta) {
    return 0.0L;
}

 real g_theta_theta (real m, real a, real r, real theta) {
    return 1.0L;
}

 real g_theta_phi (real m, real a, real r, real theta) {
    return 0.0L;
}

 real g_phi_phi (real m, real a, real r, real theta) {
    return 1.0L;
}

 real i_t_t (real m, real a, real r, real theta) {
    return 1.0L;
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
    return 1.0L;
}

 real i_r_theta (real m, real a, real r, real theta) {
    return 0.0L;
}

 real i_r_phi (real m, real a, real r, real theta) {
    return 0.0L;
}

 real i_theta_theta (real m, real a, real r, real theta) {
    return 1.0L;
}

 real i_theta_phi (real m, real a, real r, real theta) {
    return 0.0L;
}

 real i_phi_phi (real m, real a, real r, real theta) {
    return 1.0L;
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

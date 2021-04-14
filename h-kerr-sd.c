/*
 * Kerr metric geodesics using Wilkins' equations with a "pseudo-Hamiltonian" approach
 *
 * Example:  ./h-kerr-sd-dbg 6 8 .01 10000 0 0.8 1.0 1.0 0.9455050956749083 1.434374509531738 1.0 7.978759958927879 12.0 63.0 >/tmp/$USER/data
 *
 * Example:  gnuplot -p -e "set terminal wxt size 600,450; splot '<cat' with lines" </tmp/$USER/data
 *
 * Example:  gnuplot -p -e "set terminal wxt size 600,450; set yrange [-240:0]; plot '/tmp/$USER/data' using 4:5 with lines, '/tmp/$USER/data' using 4:6 with lines, '/tmp/$USER/data' using 4:7 with lines"
 *
 * Example:  gnuplot -p -e "set terminal wxt size 600,450; plot '/tmp/$USER/data' using 4:8 with lines, '/tmp/$USER/data' using 4:9 with lines"
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "symplectic.h"

typedef struct {
    real M, mu2;  // central mass & particle mass (squared)
    real E, L, Q, K;  // constants of motion
    real a, a2, L2, aL, aE, a2xmu2_E2;  // global constants
    real q_t, q_r, q_theta, q_phi, p_t, p_r, p_theta, p_phi;  // coordinates & velocities
    real r2, ra2, delta, sth, sth2, cth, cth2, sigma;  // global variables
    real P, R, THETA;  // potentials
} parameters;

static void refresh (parameters *p) {
    p->r2 = p->q_r * p->q_r;
    p->ra2 = p->r2 + p->a2;
    p->delta = p->ra2 - 2.0L * p->M * p->q_r;
    p->P = p->E * p->ra2 - p->aL;
    p->R = p->P * p->P - p->delta * (p->r2 * p->mu2 + p->K);
    p->sth = sinl(p->q_theta);
    p->sth2 = p->sth * p->sth;
    p->cth = cosl(p->q_theta);
    p->cth2 = 1.0L - p->sth2;
    p->THETA = p->Q - p->cth2 * (p->a2xmu2_E2 + p->L2 / p->sth2);
    p->p_t = p->a * (p->L - p->aE * p->sth2) + p->ra2 * p->P / p->delta;
    p->p_phi = (p->L / p->sth2 - p->aE) + p->a * p->P / p->delta;
}

static parameters *get_p (int argc, char **argv, int va_begin) {
    parameters *p = malloc(sizeof (parameters));
    real spin, bh_mass, p_mass2, energy, momentum, m_factor, carter, r_0, theta_0;
    t_variables(argv, va_begin, argc, &spin, &bh_mass, &p_mass2, &energy, &momentum, &m_factor, &carter, &r_0, &theta_0);
    p->M = bh_mass;  // constants
    p->mu2 = p_mass2;
    p->E = energy;
    p->L = momentum * m_factor;
    p->Q = carter * m_factor;
    p->a = spin;
    p->a2 = p->a * p->a;
    p->L2 = p->L * p->L;
    p->aL = p->a * p->L;
    p->aE = p->a * p->E;
    p->K = p->Q + (p->L - p->aE) * (p->L - p->aE);
    p->a2xmu2_E2 = p->a2 * (p->mu2 - p->E * p->E);
    p->q_t = 0.0L;  // coordinates
    p->q_r = r_0;
    p->q_theta = elevation_to_colatitude(theta_0);
    p->q_phi = 0.0L;
    refresh(p);  // variables, t & phi velocities
    p->p_r = - sqrtl(p->R >= 0.0L ? p->R : - p->R);
    p->p_theta = - sqrtl(p->THETA >= 0.0L ? p->THETA : - p->THETA);
    return p;
}

void update_q (void *params, real c) {  // dq / dt = dH / dp
    parameters *p = (parameters *)params;
    p->q_t += c * p->p_t;
    p->q_r += c * p->p_r;
    p->q_theta += c * p->p_theta;
    p->q_phi += c * p->p_phi;
    refresh(p);
}

void update_p (void *params, real d) {  // dp / dt = - dH / dq = - (- 0.5 dX / dq) where X is R or THETA
    parameters *p = (parameters *)params;
    p->p_r += d * (2.0L * p->E * p->q_r * p->P - (p->q_r - p->M) * (p->mu2 * p->r2 + p->K) - p->mu2 * p->q_r * p->delta);
    p->p_theta += d * (p->cth * p->sth * (p->a2xmu2_E2 + p->L2 / p->sth2) + p->cth * p->cth2 * p->L2 / (p->sth * p->sth2));
}

static real v_dot_v (real vt, real vr, real vth, real vph, real a, real ra2, real sth2, real sigma, real delta) {  // conserved
    real v1 = a * vt / sigma - ra2 * vph / sigma;
    real v4 = vt / sigma - a * sth2 * vph / sigma;
    return sth2 / sigma * v1 * v1 + vr * vr / delta / sigma + vth * vth / sigma - delta / sigma * v4 * v4;
}

static void plot_path (long dp, void *params, real t) {
    parameters *p = (parameters *)params;
    p->sigma = p->r2 + p->a2 * p->cth2;
    real e4v = error(1.0L + v_dot_v(p->p_t, p->p_r, p->p_theta, p->p_phi, p->a, p->ra2, p->sth2, p->sigma, p->delta));
    real eR = error(0.5L * (p->p_r * p->p_r - p->R));                  // "H" = p_r^2 / 2 + (- R(q_r) / 2) = 0
    real eTHETA = error(0.5L * (p->p_theta * p->p_theta - p->THETA));  // "H" = p_theta^2 / 2 + (- THETA(q_theta) / 2) = 0
    real ra = sqrtl(p->ra2);
    real gamma = p->p_t / p->sigma;
    char fs[128];
    sprintf(fs, "%%+.%ldLe %%+.%ldLe %%+.%ldLe %%.6Le  %%.3Le %%.3Le %%.3Le  %%.3Le %%.3Le\n", dp, dp, dp);
    printf(fs, ra * p->sth * cosl(p->q_phi), ra * p->sth * sinl(p->q_phi), p->q_r * cosl(p->q_theta), t,
           e4v, eR, eTHETA, gamma, sqrtl(1.0L - 1.0L / (gamma * gamma)));
}

static void plot_view (long dp, void *params, real t) {
    parameters *p = (parameters *)params;
    char fs[256];
    sprintf(fs, "%%.6Le 2  %%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.%ldLe  %%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.%ldLe",
            dp, dp, dp, dp, dp, dp, dp, dp);
    printf(fs, t, p->q_r, cosl(p->q_theta), p->q_t, p->q_phi, p->p_r, - sinl(p->q_theta) * p->p_theta, p->p_t, p->p_phi);
    printf("  -1 0 0 0  0 0 0 1  0 1 0 0\n");
}

static void plot_raw (long dp, void *params, real time) {
    (void)dp;
    parameters *p = (parameters *)params;
    printf("%.6Le  %+La %+La %+La %+La  %+La %+La %+La %+La\n",
           time, p->q_t, p->q_r, p->q_theta, p->q_phi, p->p_t, p->p_r, p->p_theta, p->p_phi);
}

int main (int argc, char **argv) {
    assert(argc == 15);
    int plot_type_position = 5;
    long plot_type = strtol(argv[plot_type_position], NULL, 10);
    plotter plot;
    switch (plot_type) {
        case 0: plot = plot_path; break;
        case 1: plot = plot_view; break;
        case 2: plot = plot_raw; break;
        default:
            printf("Plot type is {%ld} but should be 0 (x,y,z,error,speed), 1 (view), or 2 (raw)\n", plot_type);
            exit(1);
    }
    solve(argv, get_p(argc, argv, plot_type_position + 1), update_q, update_p, plot);
    return 0;
}

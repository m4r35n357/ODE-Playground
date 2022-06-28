/*
 * Kerr metric geodesics using Wilkins' equations with a "pseudo-Hamiltonian" approach with automatic differentiation
 * Separation of R and THETA equations is enabled by using non-affine (Mino) time
 *
 * Example:  ./h-kerr-dbg 6 8 .01 10000 0 0.8 1.0 0.9455050956749083 1.434374509531738 1.0 7.978759958927879 12.0 63.0 >/tmp/$USER/data
 *
 gnuplot -p << EOF
set terminal wxt size 600,450
splot '/tmp/$USER/data' with lines
EOF
 *
 gnuplot -p << EOF
set key left
set terminal wxt size 600,450
set yrange [-240:0]
set xlabel 'time'
set ylabel 'errors'
plot '/tmp/$USER/data' using 4:5 title 'v4' with lines, '' u 4:6 t 'r' w l, '' u 4:7 t 'theta' w l
EOF
 *
 gnuplot -p << EOF
set key left
set terminal wxt size 600,450
set xlabel 'time'
set ylabel 'gamma & speed'
plot '/tmp/$USER/data' using 4:8 title 'gamma' with lines, '' u 4:9 t 'speed' w l
EOF
 *
 ./h-kerr-static $(yad --columns=2 --title="Kerr Orbit" --form --separator=" " --align=right \
    --field="Display Places":NUM --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="Plot type":CB --field="BH spin" --field="particle mass" \
    --field="particle energy" --field="particle momentum" --field="momentum factor" --field="Carter's constant" \
    --field="r0" --field="theta0" \
    -- '6!0..36!3' '4!2..10!2' '0.01!0.001..1!0.001!3' '10000!1..100000!1' '0!1!2' "0.8" "1.0" \
       "0.9455050956749083" "1.434374509531738" "1.0" "7.978759958927879" "12.0" "63.0")
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "symplectic.h"
#include "dual.h"

typedef struct {
    real mu2;  // central mass & particle mass (squared)
    real E, L, Q, K;  // constants of motion
    real a, a2, L2, aL, aE, a2xmu2_E2;  // global constants
    real q_t, q_r, q_theta, q_phi, p_t, p_r, p_theta, p_phi;  // coordinates & velocities
    dual ra2, delta, sth2, R, THETA;  // global variables & potentials
} parameters;

static void refresh (parameters *p) {
    dual r = d_var(p->q_r);
    dual r2 = d_sqr(r);
    p->ra2 = d_shift(r2, p->a2);
    p->delta = d_sub(p->ra2, d_scale(r, 2.0L));
    dual P = d_shift(d_scale(p->ra2, p->E), - p->aL);
    p->R = d_sub(d_sqr(P), d_mul(p->delta, d_shift(d_scale(r2, p->mu2), p->K)));
    p->sth2 = d_sqr(d_sin(d_var(p->q_theta)));
    p->THETA = d_shift(d_mul(d_shift(d_scale(d_inv(p->sth2), p->L2), p->a2xmu2_E2), d_shift(p->sth2, - 1.0L)), p->Q);
    p->p_t = p->a * (p->L - p->aE * p->sth2.val) + p->ra2.val * P.val / p->delta.val;
    p->p_phi = (p->L / p->sth2.val - p->aE) + p->a * P.val / p->delta.val;
}

void *get_p (int argc, char **argv, int va_begin) {
    parameters *p = malloc(sizeof (parameters));
    real spin, p_mass, energy, momentum, m_factor, carter, r_0, theta_0;
    t_variables(argv, va_begin, argc, &spin, &p_mass, &energy, &momentum, &m_factor, &carter, &r_0, &theta_0);
    p->mu2 = p_mass * p_mass;  // constants
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
    p->p_r = - sqrtl(p->R.val >= 0.0L ? p->R.val : - p->R.val);
    p->p_theta = - sqrtl(p->THETA.val >= 0.0L ? p->THETA.val : - p->THETA.val);
    return p;
}

void update_q (void *params, real c) {  // dq / dt = d"H" / dp
    parameters *p = (parameters *)params;
    p->q_t += c * p->p_t;
    p->q_r += c * p->p_r;
    p->q_theta += c * p->p_theta;
    p->q_phi += c * p->p_phi;
    refresh(p);
}

void update_p (void *params, real d) {  // dp / dt = - d"H" / dq = - (- 0.5 dX / dq) where X is R or THETA
    parameters *p = (parameters *)params;
    p->p_r += d * 0.5L * p->R.dot;
    p->p_theta += d * 0.5L * p->THETA.dot;
}

static real v_dot_v (real vt, real vr, real vth, real vph, real a, real ra2, real sth2, real sigma, real delta) {  // conserved
    real v1 = a * vt / sigma - ra2 * vph / sigma;
    real v4 = vt / sigma - a * sth2 * vph / sigma;
    return sth2 / sigma * v1 * v1 + vr * vr / delta / sigma + vth * vth / sigma - delta / sigma * v4 * v4;
}

static void plot_path (int dp, void *params, real t) {
    parameters *p = (parameters *)params;
    real sigma = p->q_r * p->q_r + p->a * p->a * (1.0L - p->sth2.val);
    real ra_sth = sqrtl(p->ra2.val) * sinl(p->q_theta);
    real gamma = p->p_t / sigma;
    printf("%+.*Le %+.*Le %+.*Le  %.6Le %+.*Le%+.*Le %+.*Le  %+.*Le %+.*Le\n",
           dp, ra_sth * cosl(p->q_phi), dp, ra_sth * sinl(p->q_phi), dp, p->q_r * cosl(p->q_theta), t,
           dp, error(1.0L + v_dot_v(p->p_t, p->p_r, p->p_theta, p->p_phi, p->a, p->ra2.val, p->sth2.val, sigma, p->delta.val)),
           dp, error(0.5L * (p->p_r * p->p_r - p->R.val)),              // "H" = p_r^2 / 2 + (- R(q_r) / 2) = 0
           dp, error(0.5L * (p->p_theta * p->p_theta - p->THETA.val)),  // "H" = p_theta^2 / 2 + (- THETA(q_theta) / 2) = 0
           dp, gamma, dp, sqrtl(1.0L - 1.0L / (gamma * gamma)));
}

static void plot_view (int dp, void *params, real t) {
    parameters *p = (parameters *)params;
    printf("%.6Le 2  %+.*Le %+.*Le %+.*Le %+.*Le  %+.*Le %+.*Le %+.*Le %+.*Le  -1 0 0 0  0 0 0 1  0 1 0 0\n",
           t, dp, p->q_r, dp, cosl(p->q_theta), dp, p->q_t, dp, p->q_phi,
           dp, p->p_r, dp, - sinl(p->q_theta) * p->p_theta, dp, p->p_t, dp, p->p_phi);
}

static void plot_raw (int dp, void *params, real time) {
    (void)dp;
    parameters *p = (parameters *)params;
    printf("%.6Le  %+La %+La %+La %+La  %+La %+La %+La %+La\n",
           time, p->q_t, p->q_r, p->q_theta, p->q_phi, p->p_t, p->p_r, p->p_theta, p->p_phi);
}

int main (int argc, char **argv) {
    assert(argc == 14);
    int plot_type_position = 5;
    long plot_type = strtol(argv[plot_type_position], NULL, 10);
    plotter plot;
    switch (plot_type) {
        case 0: plot = plot_path; break;  // for plot3d.py
        case 1: plot = plot_view; break;  // for kerr-image
        case 2: plot = plot_raw; break;   // for debugging
        default:
            printf("Plot type is {%ld} but should be 0 (x,y,z,error,speed), 1 (view), or 2 (raw)\n", plot_type);
            exit(2);
    }
    solve(argv, get_p(argc, argv, plot_type_position + 1), plot);
    return 0;
}

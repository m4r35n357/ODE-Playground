/*
 * Kerr metric geodesics using Wilkins' equations with a "pseudo-Hamiltonian" approach with automatic differentiation
 *
 * Example:  ./h-kerr-dbg  6 4 .01 10000  .8 1 1 0.94550509567490792 1.4343745095317371 7.9787599589278697 12 .001 >/tmp/$USER/data
 * Example:  gnuplot -p -e "set terminal wxt size 600,450; set yrange [-240:0]; plot '/tmp/$USER/data' using 4:5 with lines, '/tmp/$USER/data' using 4:6 with lines, '/tmp/$USER/data' using 4:7 with lines"
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "dual.h"

static real four_V (real Ut, real Ur, real Uth, real Uph, real a, real ra2, real sth2, real sigma, real delta) {  // conserved
    real U1 = a * Ut / sigma - ra2 * Uph / sigma;
    real U4 = Ut / sigma - a * sth2 * Uph / sigma;
    return sth2 / sigma * U1 * U1 + Ur * Ur / delta / sigma + Uth * Uth / sigma - delta / sigma * U4 * U4;
}

typedef struct {
    real m;  // central mass
    real mu2;  // particle mass (squared)
    real E, L, Q, K;  // constants of motion
    real a, a2, L2, aL, aE, a2xmu2_E2;  // global constants
    real q_t, q_r, q_theta, q_phi;  // coordinates
    dual ra2, delta, sth, sth2;  // global variables
    real sigma;  // global variables
    real p_t, p_r, p_theta, p_phi;  // velocities
    dual R, THETA;  // potentials
} parameters;

static void refresh (parameters *p) {
    dual r = d_var(p->q_r);
    dual r2 = d_sqr(r);
    p->ra2 = d_shift(r2, p->a2);
    p->delta = d_sub(p->ra2, d_scale(r, 2.0L * p->m));
    dual P = d_shift(d_scale(p->ra2, p->E), - p->aL);
    p->R = d_sub(d_sqr(P), d_mul(p->delta, d_shift(d_scale(r2, p->mu2), p->K)));
    dual theta = d_var(p->q_theta);
    p->sth = d_sin(theta);
    p->sth2 = d_sqr(p->sth);
    dual cth2 = d_shift(d_neg(p->sth2), 1.0L);
    p->THETA = d_shift(d_neg(d_mul(d_shift(d_scale(d_inv(p->sth2), p->L2), p->a2xmu2_E2), cth2)), p->Q);
    p->p_t = p->a * (p->L - p->aE * p->sth2.val) + p->ra2.val * P.val / p->delta.val;
    p->p_phi = (p->L / p->sth2.val - p->aE) + p->a * P.val / p->delta.val;
    p->sigma = r2.val + p->a2 * cth2.val;
}

static parameters *get_p (int argc, char **argv) {
    parameters *p = malloc(sizeof (parameters));
    real spin, bh_mass, p_mass2, energy, momentum, carter, r_0, theta_0;
    t_variables(argv, argc, &spin, &bh_mass, &p_mass2, &energy, &momentum, &carter, &r_0, &theta_0);
    p->m = bh_mass;  // constants
    p->mu2 = p_mass2;
    p->E = energy;
    p->L = momentum;
    p->Q = carter;
    p->a = spin;
    p->a2 = p->a * p->a;
    p->L2 = p->L * p->L;
    p->aL = p->a * p->L;
    p->aE = p->a * p->E;
    p->K = p->Q + (p->L - p->aE) * (p->L - p->aE);
    p->a2xmu2_E2 = p->a2 * (p->mu2 - p->E * p->E);
    p->q_t = 0.0L;  // coordinates
    p->q_r = r_0;
    p->q_theta = (90.0L - theta_0) * M_PI  / 180.0L;
    p->q_phi = 0.0L;
    refresh(p);  // variables, t & phi velocities
    p->p_r = sqrtl(p->R.val >= 0.0L ? p->R.val : - p->R.val);
    p->p_theta = sqrtl(p->THETA.val >= 0.0L ? p->THETA.val : - p->THETA.val);
    return p;
}

static void update_q (void *params, real c) {  // dq / dt = dH / dp
    parameters *p = (parameters *)params;
    p->q_t += c * p->p_t;
    p->q_r += c * p->p_r;
    p->q_theta += c * p->p_theta;
    p->q_phi += c * p->p_phi;
    refresh(p);
}

static void update_p (void *params, real d) {  // dp / dt = - dH / dq = - (- 0.5 dX / dq) where X is R or THETA
    parameters *p = (parameters *)params;
    p->p_r += 0.5L * d * p->R.dot;
    p->p_theta += 0.5L * d * p->THETA.dot;
}

static void plot (long dp, void *params, real t) {
    parameters *p = (parameters *)params;
    real e4v = error(p->mu2 + four_V(p->p_t, p->p_r, p->p_theta, p->p_phi, p->a, p->ra2.val, p->sth2.val, p->sigma, p->delta.val));
    real eR = error(0.5L * (p->p_r * p->p_r - p->R.val) / (p->sigma * p->sigma));  // "H" = p_r^2 / 2 + (- R(q_r) / 2) = 0
    real eTHETA = error(0.5L * (p->p_theta * p->p_theta - p->THETA.val) / (p->sigma * p->sigma)); // similar to above
    real ra = sqrtl(p->ra2.val);
    char fs[128];
    sprintf(fs, "%%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.6Le %%+.3Le %%+.3Le %%+.3Le\n", dp, dp, dp);
    printf(fs, ra * p->sth.val * cosl(p->q_phi), ra * p->sth.val * sinl(p->q_phi), p->q_r * cosl(p->q_theta), t, e4v, eR, eTHETA);
}

int main (int argc, char **argv) {
    assert(argc == 13);
    solve(argv, get_p(argc, argv), update_q, update_p, plot);
    return 0;
}

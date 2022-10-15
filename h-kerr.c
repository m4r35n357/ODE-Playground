/*
 * Kerr metric geodesics using Wilkins' equations with a "pseudo-Hamiltonian" approach with automatic differentiation
 * Separation of R and THETA equations is enabled by using non-affine (Mino) time
 *
 * Example:  ./h-kerr-dbg 6 8 .01 10000 0 0.8 1.0 0.9455050956749083 1.434374509531738 1.0 7.978759958927879 12.0 63.0 >/tmp/$USER/data
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "symplectic.h"
#include "h-kerr.h"

real elevation_to_colatitude (real elevation) {
    return (90.0L - elevation) * MY_PI  / 180.0L;
}

real sigma (parameters *p) {
    return p->q_r * p->q_r + p->a * p->a * (1.0L - p->sth2.val);
}

pair gamma_v (parameters *p, real sigma) {
    real g = p->p_t / sigma;
    return (pair){g, sqrtl(1.0L - 1.0L / (g * g))};
}

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

parameters *get_p_kerr (int argc, char **argv) {
    fprintf(stderr, "[ "); for (int i = 0; i < argc; i++) fprintf(stderr, "%s ", argv[i]); fprintf(stderr, "]\n");
    assert(argc == 14);
    parameters *p = malloc(sizeof (parameters));
    p->step = strtold(argv[3], NULL);  // constants
    p->a = strtold(argv[6], NULL);
    real p_mass = strtold(argv[7], NULL);
    p->mu2 = p_mass * p_mass;
    p->E = strtold(argv[8], NULL);
    real m_factor = strtold(argv[10], NULL);
    p->L = strtold(argv[9], NULL) * m_factor;
    p->Q = strtold(argv[11], NULL) * m_factor;
    p->horizon = 1.0F + (float)sqrtl(1.0L - p->a * p->a);
    p->a2 = p->a * p->a;
    p->L2 = p->L * p->L;
    p->aL = p->a * p->L;
    p->aE = p->a * p->E;
    p->K = p->Q + (p->L - p->aE) * (p->L - p->aE);
    p->a2xmu2_E2 = p->a2 * (p->mu2 - p->E * p->E);
    p->q_t = p->tau = 0.0L;  // coordinates & proper time
    p->q_r = strtold(argv[12], NULL);;
    p->q_theta = elevation_to_colatitude(strtold(argv[13], NULL));
    p->q_phi = 0.0L;
    refresh(p);  // update variables, t & phi velocities
    p->p_r = - sqrtl(p->R.val >= 0.0L ? p->R.val : - p->R.val);  // potentials
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

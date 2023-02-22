/*
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "symplectic.h"
#include "h-kerr.h"

real elevation_to_colatitude (real elevation) {
    return (90.0L - elevation) * acosl(-1.0L)  / 180.0L;
}

real sigma (kerr *k) {
    return k->q_r * k->q_r + k->a * k->a * (1.0L - k->sth2.val);
}

pair gamma_v (kerr *k, real sigma) {
    real g = k->v_t / sigma;
    return (pair){g, sqrtl(1.0L - 1.0L / (g * g))};
}

static void refresh (kerr *k) {
    dual r = d_var(k->q_r);
    dual r2 = d_sqr(r);
    k->ra2 = d_shift(r2, k->a2);
    k->D = d_sub(k->ra2, d_scale(r, 2.0L));
    dual P = d_shift(d_scale(k->ra2, k->E), - k->aL);
    k->R = d_sub(d_sqr(P), d_mul(k->D, d_shift(d_scale(r2, k->mu2), k->K)));
    k->sth2 = d_sqr(d_sin(d_var(k->q_theta)));
    k->TH = d_shift(d_mul(d_shift(d_scale(d_inv(k->sth2), k->L2), k->a2xmu2_E2), d_shift(k->sth2, - 1.0L)), k->Q);
    k->v_t = k->a * (k->L - k->aE * k->sth2.val) + k->ra2.val * P.val / k->D.val;
    k->v_phi = (k->L / k->sth2.val - k->aE) + k->a * P.val / k->D.val;
}

kerr *get_p_kerr (int argc, char **argv, real step_size) {
    assert(argc == 13);
    kerr *k = malloc(sizeof (kerr));
    k->step_size = step_size;
    k->a = strtold(argv[5], NULL);  // constants
    real p_mass = strtold(argv[6], NULL);
    k->mu2 = p_mass * p_mass;
    k->E = strtold(argv[7], NULL);
    real m_factor = strtold(argv[9], NULL);
    k->L = strtold(argv[8], NULL) * m_factor;
    k->Q = strtold(argv[10], NULL) * m_factor;
    k->horizon = 1.0F + (float)sqrtl(1.0L - k->a * k->a);
    k->a2 = k->a * k->a;
    k->L2 = k->L * k->L;
    k->aL = k->a * k->L;
    k->aE = k->a * k->E;
    k->K = k->Q + (k->L - k->aE) * (k->L - k->aE);
    k->a2xmu2_E2 = k->a2 * (k->mu2 - k->E * k->E);
    k->q_t = k->tau = 0.0L;  // coordinates & proper time
    k->q_r = strtold(argv[11], NULL);;
    k->q_theta = elevation_to_colatitude(strtold(argv[12], NULL));
    k->q_phi = 0.0L;
    refresh(k);  // update variables, t & phi velocities
    k->v_r = - sqrtl(k->R.val >= 0.0L ? k->R.val : - k->R.val);  // potentials
    k->v_theta = - sqrtl(k->TH.val >= 0.0L ? k->TH.val : - k->TH.val);
    return k;
}

void update_q (void *params, real c) {  // dq / dt = d"H" / dp
    kerr *k = (kerr *)params;
    k->q_t += c * k->v_t;
    k->q_r += c * k->v_r;
    k->q_theta += c * k->v_theta;
    k->q_phi += c * k->v_phi;
    refresh(k);
}

void update_p (void *params, real d) {  // dp / dt = - d"H" / dq = - (- 0.5 dX / dq) where X is R or THETA
    kerr *k = (kerr *)params;
    k->v_r += d * 0.5L * k->R.dot;
    k->v_theta += d * 0.5L * k->TH.dot;
}

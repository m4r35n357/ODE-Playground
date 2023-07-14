/*
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "symplectic.h"
#include "h-kerr.h"

real elevation_to_colatitude (real elevation) {
    return (90.0L - elevation) * acosl(-1.0L)  / 180.0L;
}

real sigma (parameters *k) {
    return k->q_r * k->q_r + k->a * k->a * (1.0L - k->sth2.val);
}

pair gamma_v (parameters *k, real sigma) {
    real g = k->v_t / sigma;
    return (pair){g, sqrtl(1.0L - 1.0L / (g * g))};
}

static void refresh (parameters *k) {
    dual r = d_var(k->q_r);
    dual r2 = d_sqr(r);
    k->ra2 = d_shift(r2, k->a2);
    k->D = d_sub(k->ra2, d_scale(r, 2.0L));
    dual P = d_shift(d_scale(k->ra2, k->E), - k->aL);
    k->R = d_sub(d_sqr(P), d_mul(k->D, d_shift(d_scale(r2, k->mu2), k->K)));
    k->sth2 = d_sqr(d_sin(d_var(k->q_th)));
    k->TH = d_shift(d_mul(d_shift(d_scale(d_inv(k->sth2), k->L2), k->a2xmu2_E2), d_shift(k->sth2, - 1.0L)), k->Q);
    k->v_t = k->a * (k->L - k->aE * k->sth2.val) + k->ra2.val * P.val / k->D.val;
    k->v_ph = (k->L / k->sth2.val - k->aE) + k->a * P.val / k->D.val;
}

parameters *kerr_get_p (int argc, char **argv, real step_size) {
    CHECK(argc == 13);
    parameters *k = malloc(sizeof (parameters)); CHECK(k);
    k->step_size = step_size;
    k->a = strtold(argv[5], NULL);          CHECK(k->a >= -1.0L && k->a <= 1.0L);  // constants
    k->mu2 = strtold(argv[6], NULL) == 0.0L ? 0.0L : 1.0L;
    k->E = strtold(argv[7], NULL);          CHECK(k->E >= 0.0L);
    real m_factor = strtold(argv[9], NULL); CHECK(m_factor >= 0.0L && m_factor <= 1.0L);
    k->L = strtold(argv[8], NULL) * m_factor;
    k->Q = strtold(argv[10], NULL) * m_factor;
    k->a2 = k->a * k->a;
    k->horizon = 1.0F + (float)sqrtl(1.0L - k->a2);
    k->L2 = k->L * k->L;
    k->aL = k->a * k->L;
    k->aE = k->a * k->E;
    k->K = k->Q + (k->L - k->aE) * (k->L - k->aE);
    k->a2xmu2_E2 = k->a2 * (k->mu2 - k->E * k->E);
    k->q_t = k->tau = 0.0L;  // coordinates & proper time
    k->q_r = strtold(argv[11], NULL);;
    k->q_th = elevation_to_colatitude(strtold(argv[12], NULL));
    k->q_ph = 0.0L;
    refresh(k);  // update variables, t & phi velocities
    k->v_r = - sqrtl(k->R.val >= 0.0L ? k->R.val : - k->R.val);  // potentials
    k->v_th = - sqrtl(k->TH.val >= 0.0L ? k->TH.val : - k->TH.val);
    return k;
}

void update_q (parameters *p, real c) {  // dq/dt = d"H"/dp
    p->q_t  += c * p->v_t;
    p->q_r  += c * p->v_r;
    p->q_th += c * p->v_th;
    p->q_ph += c * p->v_ph;
    refresh(p);
}

void update_p (parameters *p, real d) {  // dp/dt = - d"H"/dq = - (- 0.5 dX/dq) where X is R or THETA
    p->v_r  += 0.5L * d * p->R.dot;
    p->v_th += 0.5L * d * p->TH.dot;
}

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

real sigma (model *_) {
    return _->q_r * _->q_r + _->a * _->a * (1.0L - _->sth2.val);
}

pair gamma_v (model *_, real sigma) {
    real g = _->v_t / sigma;
    return (pair){g, sqrtl(1.0L - 1.0L / (g * g))};
}

static void refresh (model *_) {
    dual r = d_var(_->q_r);
    dual r2 = d_sqr(r);
    _->ra2 = d_shift(r2, _->a2);
    _->D = d_sub(_->ra2, d_scale(r, 2.0L));
    dual P = d_shift(d_scale(_->ra2, _->E), - _->aL);
    _->R = d_sub(d_sqr(P), d_mul(_->D, d_shift(d_scale(r2, _->mu2), _->K)));
    _->sth2 = d_sqr(d_sin(d_var(_->q_th)));
    _->TH = d_shift(d_mul(d_shift(d_scale(d_inv(_->sth2), _->L2), _->a2xmu2_E2), d_shift(_->sth2, - 1.0L)), _->Q);
    _->v_t = _->a * (_->L - _->aE * _->sth2.val) + _->ra2.val * P.val / _->D.val;
    _->v_ph = (_->L / _->sth2.val - _->aE) + _->a * P.val / _->D.val;
}

model *kerr_get_p (int argc, char **argv, real step_size) {
    CHECK(argc == 13);
    model *_ = malloc(sizeof (model)); CHECK(_);
    _->step_size = step_size;
    _->a = strtold(argv[5], NULL);          CHECK(_->a >= -1.0L && _->a <= 1.0L);  // constants
    _->mu2 = strtold(argv[6], NULL) == 0.0L ? 0.0L : 1.0L;
    _->E = strtold(argv[7], NULL);          CHECK(_->E >= 0.0L);
    real m_factor = strtold(argv[9], NULL); CHECK(m_factor >= 0.0L && m_factor <= 1.0L);
    _->L = strtold(argv[8], NULL) * m_factor;
    _->Q = strtold(argv[10], NULL) * m_factor;
    _->a2 = _->a * _->a;
    _->horizon = 1.0F + (float)sqrtl(1.0L - _->a2);
    _->L2 = _->L * _->L;
    _->aL = _->a * _->L;
    _->aE = _->a * _->E;
    _->K = _->Q + (_->L - _->aE) * (_->L - _->aE);
    _->a2xmu2_E2 = _->a2 * (_->mu2 - _->E * _->E);
    _->q_t = _->tau = 0.0L;  // coordinates & proper time
    _->q_r = strtold(argv[11], NULL);;
    _->q_th = elevation_to_colatitude(strtold(argv[12], NULL));
    _->q_ph = 0.0L;
    refresh(_);  // update variables, t & phi velocities
    _->v_r = - sqrtl(_->R.val >= 0.0L ? _->R.val : - _->R.val);  // potentials
    _->v_th = - sqrtl(_->TH.val >= 0.0L ? _->TH.val : - _->TH.val);
    return _;
}

void update_q (model *_, real c) {  // dq/dt = d"H"/dp
    _->q_t  += c * _->v_t;
    _->q_r  += c * _->v_r;
    _->q_th += c * _->v_th;
    _->q_ph += c * _->v_ph;
    refresh(_);
}

void update_p (model *_, real d) {  // dp/dt = - d"H"/dq = - (- 0.5 dX/dq) where X is R or THETA
    _->v_r  += 0.5L * d * _->R.dot;
    _->v_th += 0.5L * d * _->TH.dot;
}

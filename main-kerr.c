/*
 * Kerr metric geodesics using Wilkins' equations with a "pseudo-Hamiltonian" approach with automatic differentiation
 * Separation of R and THETA equations is enabled by using non-affine (Mino) time
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "symplectic.h"
#include "h-kerr.h"

static real v_dot_v (real vt, real vr, real vth, real vph, real a, real ra2, real sth2, real sigma, real delta) {  // conserved
    real v1 = a * vt / sigma - ra2 * vph / sigma;
    real v4 = vt / sigma - a * sth2 * vph / sigma;
    return sth2 / sigma * v1 * v1 + vr * vr / delta / sigma + vth * vth / sigma - delta / sigma * v4 * v4;
}

static void plot (int dp, void *params, real mino) {
    kerr *p = (kerr *)params;
    real S = sigma(p);
    pair Y = gamma_v(p, S);
    p->tau += p->step_size * S;
    real ra_sth = sqrtl(p->ra2.val) * sinl(p->q_theta);
    printf("%+.*Le %+.*Le %+.*Le  %.6Le %+.*Le %+.*Le %+.*Le  %+.*Le %+.*Le  %.6Le %.6Le\n",
           dp, ra_sth * cosl(p->q_phi), dp, ra_sth * sinl(p->q_phi), dp, p->q_r * cosl(p->q_theta), mino,
           dp, error(1.0L + v_dot_v(p->p_t, p->p_r, p->p_theta, p->p_phi, p->a, p->ra2.val, p->sth2.val, S, p->delta.val)),
           dp, error(0.5L * (p->p_r * p->p_r - p->R.val)),              // "H" = p_r^2 / 2 + (- R(r) / 2) = 0
           dp, error(0.5L * (p->p_theta * p->p_theta - p->THETA.val)),  // "H" = p_theta^2 / 2 + (- THETA(theta) / 2) = 0
           dp, Y.a, dp, Y.b, p->tau, p->q_t);
}

int main (int argc, char **argv) {
    controls *c = get_c_symp(argc, argv);
    solve(argv, c, get_p_kerr(argc, argv, c->step_size), plot);
    return 0;
}

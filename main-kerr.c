/*
 * Kerr metric geodesics using Wilkins' equations with a "pseudo-Hamiltonian" approach with automatic differentiation
 * Separation of R and THETA equations is enabled by using non-affine (Mino) time
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "symplectic.h"
#include "h-kerr.h"

static real v2 (real vt, real vr, real vth, real vph, real a, real ra2, real sth2, real S, real D) {
    real va = (a * vt - ra2 * vph) / S;
    real vb = (vt - a * sth2 * vph) / S;
    return SQR(va) * sth2 / S + SQR(vr) / D / S + SQR(vth) / S -  SQR(vb) * D / S;
}

static void plot (int dp, model *p, real mino) {
    real S = sigma(p);
    pair Y = gamma_v(p, S);
    p->tau += p->step_size * S;
    real ra_sth = sqrtl(p->ra2.val) * sinl(p->q_th);
    printf("% .*Le % .*Le % .*Le  %.6Le % .*Le % .*Le % .*Le  % .*Le % .*Le  %.6Le %.6Le\n",
           dp, ra_sth * cosl(p->q_ph), dp, ra_sth * sinl(p->q_ph), dp, p->q_r * cosl(p->q_th), mino,
           dp, error(1.0L + v2(p->v_t, p->v_r, p->v_th, p->v_ph, p->a, p->ra2.val, p->sth2.val, S, p->D.val)),
           dp, error(0.5L * (SQR(p->v_r) - p->R.val)),      // "H" = p_r^2 / 2 + (- R(r) / 2) = 0
           dp, error(0.5L * (SQR(p->v_th) - p->TH.val)),   // "H" = p_th^2 / 2 + (- TH(th) / 2) = 0
           dp, Y.a, dp, Y.b, p->tau, p->q_t);
}

int main (int argc, char **argv) {
    controls *c = symp_get_c(argc, argv);
    solve(c, kerr_get_p(argc, argv, c->h), plot);
    return 0;
}

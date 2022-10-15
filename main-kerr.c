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
#include <math.h>
#include "symplectic.h"
#include "h-kerr.h"

static real v_dot_v (real vt, real vr, real vth, real vph, real a, real ra2, real sth2, real sigma, real delta) {  // conserved
    real v1 = a * vt / sigma - ra2 * vph / sigma;
    real v4 = vt / sigma - a * sth2 * vph / sigma;
    return sth2 / sigma * v1 * v1 + vr * vr / delta / sigma + vth * vth / sigma - delta / sigma * v4 * v4;
}

static void plot_path (int dp, void *params, real mino) {
    parameters *p = (parameters *)params;
    real S = sigma(p);
    real ra_sth = sqrtl(p->ra2.val) * sinl(p->q_theta);
    pair Y = gamma_v(p, S);
    p->tau += p->step * S;
    printf("%+.*Le %+.*Le %+.*Le  %.6Le %+.*Le %+.*Le %+.*Le  %+.*Le %+.*Le  %.6Le %.6Le\n",
           dp, ra_sth * cosl(p->q_phi), dp, ra_sth * sinl(p->q_phi), dp, p->q_r * cosl(p->q_theta), mino,
           dp, error(1.0L + v_dot_v(p->p_t, p->p_r, p->p_theta, p->p_phi, p->a, p->ra2.val, p->sth2.val, S, p->delta.val)),
           dp, error(0.5L * (p->p_r * p->p_r - p->R.val)),              // "H" = p_r^2 / 2 + (- R(r) / 2) = 0
           dp, error(0.5L * (p->p_theta * p->p_theta - p->THETA.val)),  // "H" = p_theta^2 / 2 + (- THETA(theta) / 2) = 0
           dp, Y.a, dp, Y.b, p->tau, p->q_t);
}

static void plot_view (int dp, void *params, real mino) {
    parameters *p = (parameters *)params;
    printf("%.6Le 2  %+.*Le %+.*Le %+.*Le %+.*Le  %+.*Le %+.*Le %+.*Le %+.*Le  -1 0 0 0  0 0 0 1  0 1 0 0\n",
           mino, dp, p->q_r, dp, cosl(p->q_theta), dp, p->q_t, dp, p->q_phi,
           dp, p->p_r, dp, - sinl(p->q_theta) * p->p_theta, dp, p->p_t, dp, p->p_phi);
}

static void plot_raw (int dp, void *params, real mino) { (void)dp;
    parameters *p = (parameters *)params;
    printf("%.6Le  %+La %+La %+La %+La  %+La %+La %+La %+La\n",
           mino, p->q_t, p->q_r, p->q_theta, p->q_phi, p->p_t, p->p_r, p->p_theta, p->p_phi);
}

int main (int argc, char **argv) {
    long plot_type = strtol(argv[5], NULL, BASE);
    plotter plot;
    switch (plot_type) {
        case 0: plot = plot_path; break;  // for plot3d.py
        case 1: plot = plot_view; break;  // for kerr-image
        case 2: plot = plot_raw; break;   // for debugging
        default:
            printf("Plot type is {%ld} but should be 0 (x,y,z,error,speed), 1 (view), or 2 (raw)\n", plot_type);
            exit(2);
    }
    solve(argv, get_c_symp(argv), get_p_kerr(argc, argv), plot);
    return 0;
}

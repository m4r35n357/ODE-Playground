/*
 * Newtonian central value problem using Hamilton's equations with automatic differentiation
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "symplectic.h"
#include "dual.h"

static dual hamiltonian (real gm, dual q_r, dual p_r, dual p_phi) {
    return d_sub(d_scale(d_add(d_sqr(p_r), d_div(d_sqr(p_phi), d_sqr(q_r))), 0.5L), d_scale(d_inv(q_r), gm));
}

struct Parameters {
    real m, q_r, p_r, q_phi, p_phi, h0;
};

parameters *symp_init_p (int argc, char **argv) { (void)argc;
    CHECK(argc == 8);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    _->m = strtold(argv[5], NULL);
    _->q_r = strtold(argv[6], NULL);
    _->q_phi = _->p_r = 0.0L;
    _->p_phi = strtold(argv[7], NULL) * _->m * sqrtl(_->q_r); // arg 7 = 1.0 gives a circular orbit
    _->h0 = hamiltonian(_->m, d_dual(_->q_r), d_dual(_->p_r), d_dual(_->p_phi)).val;
    return _;
}

void update_q (parameters *_, real c) {
    _->q_r   += c * hamiltonian(_->m, d_dual(_->q_r),  d_var(_->p_r), d_dual(_->p_phi)).dot;
    _->q_phi += c * hamiltonian(_->m, d_dual(_->q_r), d_dual(_->p_r),  d_var(_->p_phi)).dot;
}

void update_p (parameters *_, real d) {
    _->p_r -= d * hamiltonian(_->m, d_var(_->q_r), d_dual(_->p_r), d_dual(_->p_phi)).dot;
}

static void plot (int dp, parameters *_, real t) {
    real h_now = hamiltonian(_->m, d_dual(_->q_r), d_dual(_->p_r), d_dual(_->p_phi)).val;
    printf("%+.*Le %+.*Le %+.3Lf %.6Le %+.*Le %+.*Le\n",
           dp, _->q_r * sinl(_->q_phi), dp, _->q_r * cosl(_->q_phi), 0.0L, t, dp, error(h_now - _->h0), dp, h_now);
}

int main (int argc, char **argv) {
    solve(argv, symp_get_c(argc, argv), symp_init_p(argc, argv), plot);
    return 0;
}

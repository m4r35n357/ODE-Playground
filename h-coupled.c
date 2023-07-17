/*
 * Coupled oscillator system using Hamilton's equations with automatic differentiation
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "symplectic.h"
#include "dual.h"

static dual hamiltonian (real m, real k, dual q1, dual p1, dual q2, dual p2) {
    dual H = d_add(d_scale(d_sqr(p1), 0.5L / m), d_scale(d_sqr(p2), 0.5L / m));
    H = d_add(H, d_scale(d_sqr(q1), 0.5L * k));
    H = d_add(H, d_scale(d_sqr(q2), 0.5L * k));
    return d_add(H, d_scale(d_sqr(d_sub(q2, q1)), 0.5L * k));
}

struct Parameters {
    real m, k, q1, p1, q2, p2, h0;
};

parameters *symp_init_p (int argc, char **argv) { (void)argc;
    CHECK(argc == 9);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    _->m = strtold(argv[5], NULL);
    _->k = strtold(argv[6], NULL);
    _->q1 = strtold(argv[7], NULL) - 1.0L;
    _->q2 = strtold(argv[8], NULL) - 2.0L;
    _->p1 = _->p2 = 0.0L;
    _->h0 = hamiltonian(_->m, _->k, d_dual(_->q1), d_dual(_->p1), d_dual(_->q2), d_dual(_->p2)).val;
    return _;
}

void update_q (parameters *_, real c) {
    _->q1 += c * hamiltonian(_->m, _->k, d_dual(_->q1),  d_var(_->p1), d_dual(_->q2), d_dual(_->p2)).dot;
    _->q2 += c * hamiltonian(_->m, _->k, d_dual(_->q1), d_dual(_->p1), d_dual(_->q2),  d_var(_->p2)).dot;
}

void update_p (parameters *_, real d) {
    _->p1 -= d * hamiltonian(_->m, _->k,  d_var(_->q1), d_dual(_->p1), d_dual(_->q2), d_dual(_->p2)).dot;
    _->p2 -= d * hamiltonian(_->m, _->k, d_dual(_->q1), d_dual(_->p1),  d_var(_->q2), d_dual(_->p2)).dot;
}

static void plot (int dp, parameters *_, real t) {
    real h_now = hamiltonian(_->m, _->k, d_dual(_->q1), d_dual(_->p1), d_dual(_->q2), d_dual(_->p2)).val;
    printf("%+.*Le %+.*Le %+.*Le %+.*Le %.6Le %+.*Le %+.*Le\n",
            dp, _->q1 + 1.0L, dp, _->p1, dp, _->q2 + 2.0L, dp, _->p2, t, dp, error(h_now - _->h0), dp, h_now);
}

int main (int argc, char **argv) {
    solve(argv, symp_get_c(argc, argv), symp_init_p(argc, argv), plot);
    return 0;
}

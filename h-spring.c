/*
 * Mass-spring system using Hamilton's equations with automatic differentiation
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "symplectic.h"
#include "dual.h"

static dual hamiltonian (real m, real k, dual q, dual p) {
    return d_add(d_scale(d_sqr(p), 0.5L / m), d_scale(d_sqr(q), 0.5L * k));
}

struct Parameters {
    real m, k, q, p, h0;
};

parameters *symp_init_p (int argc, char **argv) { (void)argc;
    CHECK(argc == 8);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    _->m = strtold(argv[5], NULL);
    _->k = strtold(argv[6], NULL);
    _->q = strtold(argv[7], NULL) - 1.0L;
    _->p = 0.0L;
    _->h0 = hamiltonian(_->m, _->k, d_dual(_->q), d_dual(_->p)).val;
    return _;
}

void update_q (parameters *_, real c) {
    _->q += c * hamiltonian(_->m, _->k, d_dual(_->q), d_var(_->p)).dot;
}

void update_p (parameters *_, real d) {
    _->p -= d * hamiltonian(_->m, _->k, d_var(_->q), d_dual(_->p)).dot;
}

static void plot (int dp, parameters *_, real t) {
    real h_now = hamiltonian(_->m, _->k, d_dual(_->q), d_dual(_->p)).val;
    printf("%+.*Le %+.*Le %.6Le %+.*Le %+.*Le\n", dp, _->q + 1.0L, dp, _->p, t, dp, error(h_now - _->h0), dp, h_now);
}

int main (int argc, char **argv) {
    solve(argv, symp_get_c(argc, argv), symp_init_p(argc, argv), plot);
    return 0;
}

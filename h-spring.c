/*
 * Mass-spring system using Hamilton's equations with automatic differentiation
 *
 * Example:  ./h-spring-std  6 8 1 10000  1 1 1 >/tmp/$USER/data
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "symplectic.h"
#include "dual.h"

static dual hamiltonian (real M, real k, real l, dual q, dual p) {
    return d_add(d_scale(d_sqr(p), 0.5L / M), d_scale(d_sqr(d_shift(q, -l)), 0.5L * k));
}

struct Parameters {
    real m, k, l, q, p, h0;  // mass, spring constant, length, coordinate, momentum, initial hamiltonian
};

parameters *symp_init_p (int argc, char **argv) { (void)argc;
    CHECK(argc == 8);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    _->m = strtold(argv[5], NULL);
    _->k = strtold(argv[6], NULL);
    _->l = strtold(argv[7], NULL);
    _->q = _->l + 1.0L;
    _->p = 0.0L;
    _->h0 = hamiltonian(_->m, _->k, _->l, d_dual(_->q), d_dual(_->p)).val;
    return _;
}

void update_q (parameters *_, real c) {
    _->q += c * hamiltonian(_->m, _->k, _->l, d_dual(_->q), d_var(_->p)).dot;
}

void update_p (parameters *_, real d) {
    _->p -= d * hamiltonian(_->m, _->k, _->l, d_var(_->q), d_dual(_->p)).dot;
}

static void plot (int dp, parameters *_, real t) {
    real h_now = hamiltonian(_->m, _->k, _->l, d_dual(_->q), d_dual(_->p)).val;
    printf("%+.*Le %+.*Le %+.3Lf %.6Le %+.*Le %+.*Le\n", dp, _->q, dp, _->p, 0.0L, t, dp, error(h_now - _->h0), dp, h_now);
}

int main (int argc, char **argv) {
    solve(argv, symp_get_c(argc, argv), symp_init_p(argc, argv), plot);
    return 0;
}

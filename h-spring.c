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
    parameters *p = malloc(sizeof (parameters)); CHECK(p);
    p->m = strtold(argv[5], NULL);
    p->k = strtold(argv[6], NULL);
    p->l = strtold(argv[7], NULL);
    p->q = p->l + 1.0L;
    p->p = 0.0L;
    p->h0 = hamiltonian(p->m, p->k, p->l, d_dual(p->q), d_dual(p->p)).val;
    return p;
}

void update_q (parameters *p, real c) {
    p->q += c * hamiltonian(p->m, p->k, p->l, d_dual(p->q), d_var(p->p)).dot;
}

void update_p (parameters *p, real d) {
    p->p -= d * hamiltonian(p->m, p->k, p->l, d_var(p->q), d_dual(p->p)).dot;
}

static void plot (int dp, parameters *p, real t) {
    real h_now = hamiltonian(p->m, p->k, p->l, d_dual(p->q), d_dual(p->p)).val;
    printf("%+.*Le %+.*Le %+.3Lf %.6Le %+.*Le %+.*Le\n", dp, p->q, dp, p->p, 0.0L, t, dp, error(h_now - p->h0), dp, h_now);
}

int main (int argc, char **argv) {
    solve(argv, symp_get_c(argc, argv), symp_init_p(argc, argv), plot);
    return 0;
}

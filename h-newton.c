/*
 * Newtonian central value problem using Hamilton's equations with automatic differentiation
 *
 * Example:  ./h-newton-std  6 8 1 10000  1 12 .6 >/tmp/$USER/data
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "symplectic.h"
#include "dual.h"

static dual hamiltonian (real gm, dual q_r, dual p_r, dual p_phi) {
    return d_sub(d_scale(d_add(d_sqr(p_r), d_div(d_sqr(p_phi), d_sqr(q_r))), 0.5L), d_scale(d_inv(q_r), gm));
}

typedef struct Parameters {
    real m, q_r, p_r, q_phi, p_phi, h0;  // mass, coordinates, momenta, initial hamiltonian
} parameters;

void *get_p (int argc, char **argv) { (void)argc;
    CHECK(argc == 8);
    parameters *p = malloc(sizeof (parameters));
    p->m = strtold(argv[5], NULL);
    p->q_r = strtold(argv[6], NULL);
    p->q_phi = p->p_r = 0.0L;
    p->p_phi = strtold(argv[7], NULL) * p->m * sqrtl(p->q_r); // arg 7 = 1.0 gives a circular orbit
    p->h0 = hamiltonian(p->m, d_dual(p->q_r), d_dual(p->p_r), d_dual(p->p_phi)).val;
    return p;
}

void update_q (void *params, real c) {
    parameters *p = (parameters *)params;
    p->q_r += c * hamiltonian(p->m, d_dual(p->q_r), d_var(p->p_r), d_dual(p->p_phi)).dot;
    p->q_phi += c * hamiltonian(p->m, d_dual(p->q_r), d_dual(p->p_r), d_var(p->p_phi)).dot;
}

void update_p (void *params, real d) {
    parameters *p = (parameters *)params;
    p->p_r -= d * hamiltonian(p->m, d_var(p->q_r), d_dual(p->p_r), d_dual(p->p_phi)).dot;
}

static void plot (int dp, void *params, real t) {
    parameters *p = (parameters *)params;
    real h_now = hamiltonian(p->m, d_dual(p->q_r), d_dual(p->p_r), d_dual(p->p_phi)).val;
    printf("%+.*Le %+.*Le %+.3Lf %.6Le %+.*Le %+.*Le\n",
           dp, p->q_r * sinl(p->q_phi), dp, p->q_r * cosl(p->q_phi), 0.0L, t, dp, error(h_now - p->h0), dp, h_now);
}

int main (int argc, char **argv) {
    solve(argv, get_c_symp(argc, argv), get_p(argc, argv), plot);
    return 0;
}

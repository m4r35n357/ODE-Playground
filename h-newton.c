/*
 * Newtonian central value problem using Hamilton's equations with automatic differentiation
 *
 * Example:  ./h-newton-dbg  6 8 1 10000  1 12 .6 >/tmp/$USER/data
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "symplectic.h"
#include "dual.h"

static dual h (real gm, dual q_r, dual p_r, dual p_phi) {
    return d_sub(d_scale(d_add(d_sqr(p_r), d_div(d_sqr(p_phi), d_sqr(q_r))), 0.5L), d_scale(d_inv(q_r), gm));
}

typedef struct Parameters {
    real m;  // central mass
    real q_r, p_r, q_phi, p_phi;  // coordinates & momenta
    real h0;  // stored initial value of Hamiltonian
} parameters;

void *get_p (int argc, char **argv, int va_begin) { (void)argc; (void)va_begin;
    fprintf(stderr, "[ "); for (int i = 0; i < argc; i++) fprintf(stderr, "%s ", argv[i]); fprintf(stderr, "]\n");
    assert(argc == 8);
    parameters *p = malloc(sizeof (parameters));
    p->m = strtold(argv[5], NULL);
    p->q_r = strtold(argv[6], NULL);
    p->p_r = 0.0L;
    p->q_phi = 0.0L;
    p->p_phi = strtold(argv[7], NULL) * p->m * sqrtl(p->q_r);
    p->h0 = h(p->m, d_dual(p->q_r), d_dual(p->p_r), d_dual(p->p_phi)).val;
    return p;
}

void update_q (void *params, real c) {
    parameters *p = (parameters *)params;
    real q_r = c * h(p->m, d_dual(p->q_r), d_var(p->p_r), d_dual(p->p_phi)).dot;
    real q_phi = c * h(p->m, d_dual(p->q_r), d_dual(p->p_r), d_var(p->p_phi)).dot;
    p->q_r += q_r;
    p->q_phi += q_phi;
}

void update_p (void *params, real d) {
    parameters *p = (parameters *)params;
    p->p_r -= d * h(p->m, d_var(p->q_r), d_dual(p->p_r), d_dual(p->p_phi)).dot;
}

static void plot (int dp, void *params, real t) {
    parameters *p = (parameters *)params;
    real h_now = h(p->m, d_dual(p->q_r), d_dual(p->p_r), d_dual(p->p_phi)).val;
    printf("%+.*Le %+.*Le %+.3Lf %.6Le %+.*Le %+.*Le\n",
           dp, p->q_r * sinl(p->q_phi), dp, p->q_r * cosl(p->q_phi), 0.0L, t, dp, error(h_now - p->h0), dp, h_now);
}

int main (int argc, char **argv) {
    solve(argv, get_c(argv), get_p(argc, argv, 0), plot);
    return 0;
}

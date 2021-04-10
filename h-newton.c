/*
 * Newtonian central value problem using Hamilton's equations with automatic differentiation
 *
 * Example:  ./h-newton-dbg  6 8 1 10000  1 12 .6 >/tmp/$USER/data
 *
 * Example:  gnuplot -p -e "set terminal wxt size 600,450; splot '<cat' with lines" </tmp/$USER/data
 *
 * Example:  gnuplot -p -e "set terminal wxt size 600,450; set yrange [-240:0]; plot '/tmp/$USER/data' using 4:5 with lines"
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "symplectic.h"

static dual h (real m, real gm, dual q_r, dual p_r, dual p_phi) {
    return d_sub(d_scale(d_add(d_sqr(p_r), d_div(d_sqr(p_phi), d_sqr(q_r))), 0.5L / m), d_scale(d_inv(q_r), gm));
}

typedef struct {
    real m, gm;  // central mass, m X g
    real q_r, p_r, q_phi, p_phi;  // coordinates & momenta
    real h0;  // stored initial value of Hamiltonian
} parameters;

static parameters *get_p (int argc, char **argv, int va_begin) {
    parameters *p = malloc(sizeof (parameters));
    real r0, l_fac;
    t_variables(argv, va_begin, argc, &p->m, &r0, &l_fac);
    p->gm = p->m;  //p->gm = 6.6743e-11L * p->m;
    p->q_r = r0;
    p->p_r = 0.0L;
    p->q_phi = 0.0L;
    p->p_phi = l_fac * p->m * sqrtl(r0);
    p->h0 = h(p->m, p->gm, d_dual(p->q_r), d_dual(p->p_r), d_dual(p->p_phi)).val;
    return p;
}

void update_q (void *params, real c) {
    parameters *p = (parameters *)params;
    real q_r = c * h(p->m, p->gm, d_dual(p->q_r), d_var(p->p_r), d_dual(p->p_phi)).dot;
    real q_phi = c * h(p->m, p->gm, d_dual(p->q_r), d_dual(p->p_r), d_var(p->p_phi)).dot;
    p->q_r += q_r;
    p->q_phi += q_phi;
}

void update_p (void *params, real d) {
    parameters *p = (parameters *)params;
    p->p_r -= d * h(p->m, p->gm, d_var(p->q_r), d_dual(p->p_r), d_dual(p->p_phi)).dot;
}

static void plot (long dp, void *params, real t) {
    parameters *p = (parameters *)params;
    char fs[128];
    sprintf(fs, "%%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.6Le %%+.3Le %%+.3Le\n", dp, dp, dp);
    real h_now = h(p->m, p->gm, d_dual(p->q_r), d_dual(p->p_r), d_dual(p->p_phi)).val;
    printf(fs, p->q_r * sinl(p->q_phi), p->q_r * cosl(p->q_phi), 0.0L, t, error(h_now - p->h0), h_now);
}

int main (int argc, char **argv) {
    assert(argc == 8);
    solve(argv, get_p(argc, argv, 5), update_q, update_p, plot);
    return 0;
}
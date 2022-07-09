/*
 * Newtonian central value problem using Hamilton's equations with automatic differentiation
 *
 * Example:  ./h-newton-dbg  6 8 1 10000  1 12 .6 >/tmp/$USER/data
 *
 * Example:  gnuplot -p -e "set terminal wxt size 600,450; splot '/tmp/$USER/data' with lines"
 *
 * Example:  gnuplot -p -e "set terminal wxt size 600,450; plot '/tmp/$USER/data' using 1:2 with lines"
 *
 * Example:  gnuplot -p -e "set terminal wxt size 600,450; set yrange [-240:0]; plot '/tmp/$USER/data' using 4:5 with lines"
 *
 ./h-newton-static $(yad --title="Newtonian Orbit" --form --separator=" " --align=right \
    --field="Display Places":NUM --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="mass" --field="r0" --field="L factor":NUM \
    -- '6!3..64!3' '4!2..10!2' '0.2!0.1..1!0.1!1' '10000!1..100000!1' "1" "12" '0.6!0.0..1.0!0.1!1') >/tmp/$USER/data
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

typedef struct {
    real m;  // central mass
    real q_r, p_r, q_phi, p_phi;  // coordinates & momenta
    real h0;  // stored initial value of Hamiltonian
} parameters;

void *get_p (int argc, char **argv, int va_begin) {
    parameters *p = malloc(sizeof (parameters));
    real r0, l_fac;
    t_variables(argv, va_begin, argc, &p->m, &r0, &l_fac);
    p->q_r = r0;
    p->p_r = 0.0L;
    p->q_phi = 0.0L;
    p->p_phi = l_fac * p->m * sqrtl(r0);
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
    parameters *p;
    
    assert(argc == 8);
    p = get_p(argc, argv, 5);
    solve(argv, p, plot);
    while ((p = generate(argv, p))) {
        fprintf(stderr, "");
    }
    return 0;
}

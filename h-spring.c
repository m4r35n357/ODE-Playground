/*
 * Mass-spring system using Hamilton's equations with automatic differentiation
 *
 * Example:  ./h-spring-dbg  6 8 1 10000  1 1 1 >/tmp/$USER/data
 *
 * Example:  gnuplot -p -e "set terminal wxt size 600,450; splot '/tmp/$USER/data' with lines"
 *
 * Example:  gnuplot -p -e "set terminal wxt size 600,450; plot '/tmp/$USER/data' using 1:2 with lines"
 *
 * Example:  gnuplot -p -e "set terminal wxt size 600,450; set yrange [-240:0]; plot '/tmp/$USER/data' using 4:5 with lines"
 *
 ./h-spring-static $(yad --title="Mass-Spring System" --form --separator=" " --align=right \
    --field="Display Places":NUM \
    --field="Order":NUM \
    --field="Step Size":NUM \
    --field="Steps":NUM \
    --field="mass" \
    --field="spring constant" \
    --field="length" \
    -- '6!3..64!3' '4!2..10!2' '0.1!0.1..1!0.1!1' '10000!1..100000!1' "1" "1" "1") >/tmp/$USER/data
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "symplectic.h"
#include "dual.h"

static dual hamiltonian (real M, real k, real l, dual q, dual p) {
    return d_add(d_scale(d_sqr(p), 0.5L / M), d_scale(d_sqr(d_shift(q, -l)), 0.5L * k));
}

typedef struct {
    real m, k, l;  // mass, spring constant & length
    real q, p;  // coordinate & momentum
    real h0;  // stored initial value of Hamiltonian
} parameters;

void *get_p (int argc, char **argv, int va_begin) {
    parameters *p = malloc(sizeof (parameters));
    t_variables(argv, va_begin, argc, &p->m, &p->k, &p->l);
    p->q = p->l + 1.0L;
    p->p = 0.0L;
    p->h0 = hamiltonian(p->m, p->k, p->l, d_dual(p->q), d_dual(p->p)).val;
    return p;
}

void update_q (void *params, real c) {
    parameters *p = (parameters *)params;
    p->q += c * hamiltonian(p->m, p->k, p->l, d_dual(p->q), d_var(p->p)).dot;
}

void update_p (void *params, real d) {
    parameters *p = (parameters *)params;
    p->p -= d * hamiltonian(p->m, p->k, p->l, d_var(p->q), d_dual(p->p)).dot;
}

static void plot (long dp, void *params, real t) {
    parameters *p = (parameters *)params;
    char fs[128];
    sprintf(fs, "%%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.6Le %%+.3Le %%+.3Le\n", dp, dp, dp);
    real h_now = hamiltonian(p->m, p->k, p->l, d_dual(p->q), d_dual(p->p)).val;
    printf(fs, p->q, p->p, 0.0L, t, error(h_now - p->h0), h_now);
}

int main (int argc, char **argv) {
    assert(argc == 8);
    solve(argv, get_p(argc, argv, 5), plot);
    return 0;
}

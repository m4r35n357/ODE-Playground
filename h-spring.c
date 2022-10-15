/*
 * Mass-spring system using Hamilton's equations with automatic differentiation
 *
 * Example:  ./h-spring-dbg  6 8 1 10000  1 1 1 >/tmp/$USER/data
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "symplectic.h"
#include "dual.h"

static dual h (real M, real k, real l, dual q, dual p) {
    return d_add(d_scale(d_sqr(p), 0.5L / M), d_scale(d_sqr(d_shift(q, -l)), 0.5L * k));
}

typedef struct Parameters {
    real m, k, l;  // mass, spring constant & length
    real q, p;  // coordinate & momentum
    real h0;  // stored initial value of Hamiltonian
} parameters;

void *get_p (int argc, char **argv) { (void)argc;
    fprintf(stderr, "[ "); for (int i = 0; i < argc; i++) fprintf(stderr, "%s ", argv[i]); fprintf(stderr, "]\n");
    assert(argc == 8);
    parameters *p = malloc(sizeof (parameters));
    p->m = strtold(argv[5], NULL);
    p->k = strtold(argv[6], NULL);
    p->l = strtold(argv[7], NULL);
    p->q = p->l + 1.0L;
    p->p = 0.0L;
    p->h0 = h(p->m, p->k, p->l, d_dual(p->q), d_dual(p->p)).val;
    return p;
}

void update_q (void *params, real c) {
    parameters *p = (parameters *)params;
    p->q += c * h(p->m, p->k, p->l, d_dual(p->q), d_var(p->p)).dot;
}

void update_p (void *params, real d) {
    parameters *p = (parameters *)params;
    p->p -= d * h(p->m, p->k, p->l, d_var(p->q), d_dual(p->p)).dot;
}

static void plot (int dp, void *params, real t) {
    parameters *p = (parameters *)params;
    real h_now = h(p->m, p->k, p->l, d_dual(p->q), d_dual(p->p)).val;
    printf("%+.*Le %+.*Le %+.3Lf %.6Le %+.*Le %+.*Le\n",
           dp, p->q, dp, p->p, 0.0L, t, dp, error(h_now - p->h0), dp, h_now);
}

int main (int argc, char **argv) {
    solve(argv, get_c_symp(argv), get_p(argc, argv), plot);
    return 0;
}

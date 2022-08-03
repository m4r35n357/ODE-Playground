/*
 * Symplectic integration analyzer
 *
 * Example:  ./h-analysis-dbg  6 8 1 1 >/tmp/$USER/data
  *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "symplectic.h"

typedef struct Parameters { real c; real d; } parameters;

void *get_p (int argc, char **argv, int va_begin) { (void)argc; (void)argv; (void)va_begin;
    parameters *p = malloc(sizeof (parameters));
    p->c = 0.0L;
    p->d = 0.0L;
    return p;
}

static void plot (int dp, void *params, real t) { (void)dp; (void)t;
    parameters *p = (parameters *)params;
    printf("%+.3Le %+.3Le\n", p->c, p->d);
}

void update_q (void *params, real c) {
    parameters *p = (parameters *)params;
    p->c += c;
    plot(0L, p, 0.0L);
}

void update_p (void *params, real d) {
    parameters *p = (parameters *)params;
    p->d += d;
    plot(0L, p, 0.0L);
}

int main (int argc, char **argv) {
    assert(argc == 5);
    solve(argv, get_c(argv), get_p(argc, argv, 5), plot);
    return 0;
}

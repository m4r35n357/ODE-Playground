/*
 * Symplectic integration analyzer
 *
 * Example:  ./h-analysis-std  6 8 1 1 >/tmp/$USER/data
  *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "symplectic.h"

typedef struct Parameters { real c, d; } parameters;

void *get_p (int argc, char **argv) { (void)argc; (void)argv;
    CHECK(argc == 5);
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
    solve(argv, get_c_symp(argc, argv), get_p(argc, argv), plot);
    return 0;
}

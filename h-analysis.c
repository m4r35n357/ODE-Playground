/*
 * Symplectic integration analyzer
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "symplectic.h"

struct Parameters { real c, d; };

parameters *symp_init_p (int argc, char **argv) { (void)argc; (void)argv;
    CHECK(argc == 5);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    _->c = 0.0L;
    _->d = 0.0L;
    return _;
}

static void plot (int dp, parameters *p, real t) { (void)dp; (void)t;
    printf("%+.3Le %+.3Le\n", p->c, p->d);
}

void update_q (parameters *p, real c) {
    p->c += c;
    plot(0L, p, 0.0L);
}

void update_p (parameters *p, real d) {
    p->d += d;
    plot(0L, p, 0.0L);
}

int main (int argc, char **argv) {
    solve(argv, symp_get_c(argc, argv), symp_init_p(argc, argv), plot);
    return 0;
}

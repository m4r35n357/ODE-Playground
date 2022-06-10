/*
 * Symplectic integration analyzer
 *
 * Example:  ./h-analysis-dbg  6 8 1 1 >/tmp/$USER/data
 *
 * Example:  gnuplot -p -e "set terminal wxt size 600,450; plot '/tmp/$USER/data' using 2:3 with linespoints pt 7 ps 0"
  *
 ./h-analysis-static $(yad --title="Step Analysis" --form --separator=" " --align=right \
    --field="Display Places":NUM \
    --field="Order":NUM \
    --field="Step Size":RO \
    --field="Steps":RO \
    -- '6!3..64!3' '4!2..10!2' "1.0" "1") >/tmp/$USER/data
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "symplectic.h"

typedef struct {
    real c;
    real d;
    long count;
} parameters;

void *get_p (int argc, char **argv, int va_begin) {
    (void)argc; (void)argv; (void)va_begin;
    parameters *p = malloc(sizeof (parameters));
    p->c = 0.0L;
    p->d = 0.0L;
    p->count = 0L;
    return p;
}

static void plot (long dp, void *params, real t) {
    (void)dp; (void)t;
    parameters *p = (parameters *)params;
    printf("%+.6ld %+.3Le %+.3Le\n", p->count, p->c, p->d);
}

void update_q (void *params, real c) {
    parameters *p = (parameters *)params;
    p->c += c;
    p->count += 1;
    plot(0L, p, 0.0L);
}

void update_p (void *params, real d) {
    parameters *p = (parameters *)params;
    p->d += d;
    p->count += 1;
    plot(0L, p, 0.0L);
}

int main (int argc, char **argv) {
    assert(argc == 5);
    solve(argv, get_p(argc, argv, 5), plot);
    return 0;
}

/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include "ode-common.h"

controls *get_c (char **argv) {
    controls *c = malloc(sizeof (controls));
    c->order = (int)strtol(argv[2], NULL, BASE); assert(c->order >= 2 && c->order <= 64);
    c->step_size = strtold(argv[3], NULL); assert(c->step_size > 0.0L);
    c->steps = (int)strtol(argv[4], NULL, BASE); assert(c->steps >= 0 && c->steps <= 1000000);
    return c;
}

void t_params (char **argv, int argc, ...) {
    fprintf(stderr, "[ ");
    for (int i = 0; i < argc; i++) {
        fprintf(stderr, "%s ", argv[i]);
    }
    fprintf(stderr, "]\n");
    va_list model;
    va_start(model, argc);
    for (int i = 8; i < argc; i++) {
        *va_arg(model, real *) = strtold(argv[i], NULL);
    }
    va_end(model);
}

void t_out (int dp, real x, real y, real z, real t, char *x_tag, char *y_tag, char *z_tag, clock_t since) {
    real cpu = (real)(clock() - since) / CLOCKS_PER_SEC;
    if (dp == 0) {
        printf("%+La %+La %+La %.6Le %s %s %s %.3Lf\n", x, y, z, t, x_tag, y_tag, z_tag, cpu);
    } else {
        printf("%+.*Le %+.*Le %+.*Le %.6Le %s %s %s %.3Lf\n", dp, x, dp, y, dp, z, t, x_tag, y_tag, z_tag, cpu);
    }
}

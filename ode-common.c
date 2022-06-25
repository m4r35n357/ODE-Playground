/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include "taylor-ode.h"

const int BASE = 10;

void t_params (char **argv, int argc, ...) {
    va_list model;
    va_start(model, argc);
    for (int i = 8; i < argc; i++) {
        *va_arg(model, real *) = strtold(argv[i], NULL);
    }
    va_end(model);
}

void t_out (int dp, real x, real y, real z, real t, char *x_tag, char *y_tag, char *z_tag, real cpu) {
    printf("%+.*Le %+.*Le %+.*Le %.6Le %s %s %s %.3Lf\n", dp, x, dp, y, dp, z, t, x_tag, y_tag, z_tag, cpu);
}

real cpu (clock_t since) {
    return (real)(clock() - since) / CLOCKS_PER_SEC;
}

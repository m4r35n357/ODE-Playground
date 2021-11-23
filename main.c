/*
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    int display_precision = (int)strtol(argv[1], NULL, BASE);
    long double internal_precision = strtod(argv[2], NULL) * 3.322L;
    mpfr_set_default_prec((int)internal_precision);
    fprintf(stderr, " MPFR default precision: %lu bits\n", mpfr_get_default_prec());
    t_init(display_precision);

    mpfr_t step_size, x0, y0, z0;

    int order = (int)strtol(argv[3], NULL, BASE);
    mpfr_init_set_str(step_size, argv[4], BASE, RND);
    int steps = (int)strtol(argv[5], NULL, BASE);

    mpfr_init_set_str(x0, argv[6], BASE, RND);
    mpfr_init_set_str(y0, argv[7], BASE, RND);
    mpfr_init_set_str(z0, argv[8], BASE, RND);

    tsm(argc, argv, order, step_size, steps, x0, y0, z0);

    return 0;
}

/*
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long dp = strtol(argv[1], NULL, BASE);
    double precision = strtod(argv[2], NULL) * 3.322;
    mpfr_set_default_prec((int)precision);
    fprintf(stderr, " MPFR default precision: %lu bits\n", mpfr_get_default_prec());
    long n = strtol(argv[3], NULL, BASE);
    mpfr_t h, x0, y0, z0;
    mpfr_init_set_str(h, argv[4], BASE, RND);
    long steps = strtol(argv[5], NULL, BASE);
    mpfr_init_set_str(x0, argv[6], BASE, RND);
    mpfr_init_set_str(y0, argv[7], BASE, RND);
    mpfr_init_set_str(z0, argv[8], BASE, RND);

    t_tempvars(dp);
    tsm(argc, argv, n, h, steps, x0, y0, z0);

    return 0;
}

/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    PRINT_ARGS(argc, argv);

    int display_precision = (int)strtol(argv[1], NULL, BASE); CHECK(display_precision >= 2);
    int mpfr_precision_bits = (int)strtol(argv[2], NULL, BASE); CHECK(mpfr_precision_bits >= 53);
    mpfr_set_default_prec((mpfr_prec_t)mpfr_precision_bits);
    fprintf(stderr, " MPFR default precision: %lu bits\n", mpfr_get_default_prec());
    int order = (int)strtol(argv[3], NULL, BASE); CHECK(order >= 2);

    mpfr_t step_size, x0, y0, z0;
    mpfr_init_set_str(step_size, argv[4], BASE, RND); CHECK(mpfr_sgn(step_size) > 0);
    int steps = (int)strtol(argv[5], NULL, BASE); CHECK(steps >= 0 && steps <= 1000000);

    mpfr_init_set_str(x0, argv[6], BASE, RND);
    mpfr_init_set_str(y0, argv[7], BASE, RND);
    mpfr_init_set_str(z0, argv[8], BASE, RND);

    series3 *jets = malloc(sizeof (series3)); CHECK(jets);
    jets->x = t_const(order + 1, x0);
    jets->y = t_const(order + 1, y0);
    jets->z = t_const(order + 1, z0);

    tsm_init(display_precision);
    tsm(order, step_size, steps, jets, get_p(argc, argv, order), clock());

    return 0;
}

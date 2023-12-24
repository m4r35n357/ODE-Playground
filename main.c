/*
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    PRINT_ARGS(argc, argv); CHECK(argc > 9);

    int display_precision = (int)strtol(argv[1], NULL, BASE); CHECK(display_precision >= 2);
    int precision_in_bits = (int)strtol(argv[2], NULL, BASE); CHECK(precision_in_bits >= 53);
    mpfr_set_default_prec((mpfr_prec_t)precision_in_bits);
    fprintf(stderr, " %sMPFR default precision:%s %lu bits\n", GRY, NRM, mpfr_get_default_prec());
    int order = (int)strtol(argv[3], NULL, BASE); CHECK(order >= 2);

    real step_size, x0, y0, z0;
    mpfr_init_set_str(step_size, argv[4], BASE, RND); CHECK(mpfr_sgn(step_size) > 0);
    int steps = (int)strtol(argv[5], NULL, BASE); CHECK(steps >= 0 && steps <= 1000000);

    mpfr_init_set_str(x0, argv[6], BASE, RND);
    mpfr_init_set_str(y0, argv[7], BASE, RND);
    mpfr_init_set_str(z0, argv[8], BASE, RND);

    triplet *lhs = malloc(sizeof (triplet)); CHECK(lhs);
    mpfr_inits(lhs->x, lhs->y, lhs->z, NULL);

    xyz *jets = malloc(sizeof (xyz)); CHECK(jets);
    jets->x = tsm_jet(order + 1); mpfr_set(jets->x[0], x0, RND);
    jets->y = tsm_jet(order + 1); mpfr_set(jets->y[0], y0, RND);
    jets->z = tsm_jet(order + 1); mpfr_set(jets->z[0], z0, RND);

    tsm_init(display_precision);
    tsm(order, step_size, steps, lhs, jets, tsm_init_p(argc, argv, order), clock());

    return 0;
}

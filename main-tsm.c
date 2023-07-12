/*
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    int display_precision = (int)strtol(argv[1], NULL, BASE); CHECK(display_precision >= 0 && display_precision <= 32);
    controls *c = tsm_get_c(argc, argv);

    tsm_stdout(display_precision, c, tsm_init_xyz(argv, c->order), tsm_init_p(argc, argv, c->order), clock());

    return 0;
}

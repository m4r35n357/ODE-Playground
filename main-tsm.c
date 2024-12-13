/*
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    CHECK(argc > 8);

    controls *c = tsm_get_c(argc, argv);
    tsm_out(c, tsm_init(argv, c->order), tsm_init_p(argc, argv, c->order), clock());

    return 0;
}

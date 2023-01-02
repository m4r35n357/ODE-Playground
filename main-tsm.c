/*
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    int display_precision = (int)strtol(argv[1], NULL, BASE); assert(display_precision >= 0 && display_precision <= 32);
    controls *c = get_c_tsm(argc, argv);

    tsm_stdout(display_precision, c, initial_values(argv, c->order), get_p(argc, argv, c->order), clock());

    return 0;
}

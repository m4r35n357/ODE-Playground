/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    int display_precision = (int)strtol(argv[1], NULL, BASE); assert(display_precision >= 0 && display_precision <= 32);
    controls *c = get_c(argv);

    tsm(display_precision, c,
        strtold(argv[5], NULL), strtold(argv[6], NULL), strtold(argv[7], NULL), get_p(argc, argv, c->order), clock());

    return 0;
}

/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    int display_precision = (int)strtol(argv[1], NULL, BASE); assert(display_precision >= 0 && display_precision <= 32);
    int order = (int)strtol(argv[2], NULL, BASE); assert(order >= 2 && order <= 64);

    real step_size = strtold(argv[3], NULL); assert(step_size > 0.0L);
    int steps = (int)strtol(argv[4], NULL, BASE); assert(steps >= 1 && steps <= 1000000);

    tsm(display_precision, order, step_size, steps,
        strtold(argv[5], NULL), strtold(argv[6], NULL), strtold(argv[7], NULL), get_p(argc, argv, order), clock());

    return 0;
}

/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "runge-kutta-ode.h"

int main (int argc, char **argv) {
    int display_precision = (int)strtol(argv[1], NULL, BASE); assert(display_precision >= 0 && display_precision <= 32);
    controls *c = get_c(argv);
    components xyz = (components) { strtold(argv[5], NULL), strtold(argv[6], NULL), strtold(argv[7], NULL) };

    rk4(display_precision, c, &xyz, get_p(argc, argv), clock());

    return 0;
}

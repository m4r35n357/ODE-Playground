/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "runge-kutta-ode.h"

int main (int argc, char **argv) {
    int display_precision = (int)strtol(argv[1], NULL, BASE); assert(display_precision >= 1 && display_precision <= 99);
    int plot_interval = (int)strtol(argv[2], NULL, BASE); assert(plot_interval >= 1 && plot_interval <= 1000);

    real step_size = strtold(argv[3], NULL); assert(step_size > 0.0L);
    int steps = (int)strtol(argv[4], NULL, BASE); assert(steps >= 1 && steps <= 1000000);

    rk4(display_precision, plot_interval, step_size, steps,
        strtold(argv[5], NULL), strtold(argv[6], NULL), strtold(argv[7], NULL), get_p(argc, argv));

    return 0;
}

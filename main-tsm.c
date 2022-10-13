/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    int display_precision = (int)strtol(argv[1], NULL, BASE); assert(display_precision >= 0 && display_precision <= 32);
    controls *c = get_c_tsm(argv);

    series3 *jets = malloc(sizeof (series3));
    jets->x = t_jet(c->order + 1); jets->x[0] = strtold(argv[5], NULL);
    jets->y = t_jet(c->order + 1); jets->y[0] = strtold(argv[6], NULL);
    jets->z = t_jet(c->order + 1); jets->z[0] = strtold(argv[7], NULL);

    tsm_stdout(display_precision, c, jets, get_p(argc, argv, c->order), clock());

    return 0;
}

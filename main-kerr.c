/*
 * Kerr metric geodesics using Wilkins' equations with a "pseudo-Hamiltonian" approach with automatic differentiation
 * Separation of R and THETA equations is enabled by using non-affine (Mino) time
 *
 * Example:  ./h-kerr-dbg 6 8 .01 10000 0 0.8 1.0 0.9455050956749083 1.434374509531738 1.0 7.978759958927879 12.0 63.0 >/tmp/$USER/data
 *
 gnuplot -p << EOF
set terminal wxt size 600,450
splot '/tmp/$USER/data' with lines
EOF
 *
 gnuplot -p << EOF
set key left
set terminal wxt size 600,450
set yrange [-240:0]
set xlabel 'time'
set ylabel 'errors'
plot '/tmp/$USER/data' using 4:5 title 'v4' with lines, '' u 4:6 t 'r' w l, '' u 4:7 t 'theta' w l
EOF
 *
 gnuplot -p << EOF
set key left
set terminal wxt size 600,450
set xlabel 'time'
set ylabel 'gamma & speed'
plot '/tmp/$USER/data' using 4:8 title 'gamma' with lines, '' u 4:9 t 'speed' w l
EOF
 *
 ./h-kerr-static $(yad --columns=2 --title="Kerr Orbit" --form --separator=" " --align=right \
    --field="Display Places":NUM --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="Plot type":CB --field="BH spin" --field="particle mass" \
    --field="particle energy" --field="particle momentum" --field="momentum factor" --field="Carter's constant" \
    --field="r0" --field="theta0" \
    -- '6!0..36!3' '4!2..10!2' '0.01!0.001..1!0.001!3' '10000!1..100000!1' '0!1!2' "0.8" "1.0" \
       "0.9455050956749083" "1.434374509531738" "1.0" "7.978759958927879" "12.0" "63.0") >/tmp/$USER/data
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "symplectic.h"
#include "h-kerr.h"


int main (int argc, char **argv) {
    int plot_type_position = 5;
    long plot_type = strtol(argv[plot_type_position], NULL, BASE);
    plotter plot;
    switch (plot_type) {
        case 0: plot = plot_path; break;  // for plot3d.py
        case 1: plot = plot_view; break;  // for kerr-image
        case 2: plot = plot_raw; break;   // for debugging
        default:
            printf("Plot type is {%ld} but should be 0 (x,y,z,error,speed), 1 (view), or 2 (raw)\n", plot_type);
            exit(2);
    }
    solve(argv, get_c(argv), get_p(argc, argv, plot_type_position + 1), plot);
    return 0;
}

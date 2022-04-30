/*
 * Interface for ODE types and I/O functions
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include "real.h"

/*
 * The numerical base for string IO conversions
 */
extern const int BASE;

/*
 * Prints a line of data to stdout, with turning point markers
 */
void t_output (int dp, real x, real y, real z, real t, char *x_label, char *y_label, char *z_label);

/*
 * Retrieves ODE parameters from the tail of the command (arguments 8 onwards)
 */
void t_params (char **argv, int count, ...);

/*
 * For returning x, y, z velocities from the model
 */
typedef struct {
    real x;
    real y;
    real z;
} components;

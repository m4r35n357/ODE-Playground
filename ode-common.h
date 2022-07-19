/*
 * Interface for ODE types and I/O functions
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <time.h>
#include "real.h"

/*
 * The numerical base for string IO conversions
 */
extern const int BASE;

/*
 * For returning x, y, z velocities from the model
 */
typedef struct Components {
    real x, y, z;
} components;

/*
 * For returning integrator control parameters
 */
typedef struct Controls {
    long order, step, steps;
    real step_size;
} controls;

typedef struct Line {
    int newest;
    components *buffer;
} line;

typedef struct Particle {
    components *coordinates;
    components colour;
    line *track;
    float size, view_radius, view_latitude, view_longitude;
} particle;

/*
 * Retrieves integrator control parameters
 */
controls *get_c (char **argv);

/*
 * Retrieves ODE parameters from the tail of the command (arguments 8 onwards)
 */
void t_params (char **argv, int count, ...);

/*
 * Prints a line of data to stdout, with turning point markers
 */
void t_out (int dp, real x, real y, real z, real t, char *x_label, char *y_label, char *z_label, clock_t since);

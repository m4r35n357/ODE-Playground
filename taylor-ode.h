/*
 * Interface for solving systems of Ordinary Differential Equations
 * using the Taylor Series Method (TSM)
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include "real.h"

/*
 * Diffentiate real numbers from jets
 */
typedef long double *series;

/*
 * Prints an index column, and x, y, z columns, into a single line
 */
void t_output (long dp, real x, real y, real z, real t, char *x_label, char *y_label, char *z_label) ;

/*
 * Bulk set ODE parameters from the command line arguments (5 onwards)
 */
void t_params (char **argv, int count, ...);

/*
 * Creates an zeroed jet of the specified size
 */
series t_jet (long size);

/*
 * Returns a jet of the specified size, with element zero set to value and the rest zeroed (represents a constant in an ODE)
 */
series t_jet_c (long size, real value);

/*
 * Returns kth element of the absolute value of U, no user-supplied jet storage needed
 */
real t_abs (series u, int k);

/*
 * Returns kth element of the product of U and V, no user-supplied jet storage needed
 */
real t_prod (series U, series V, int k);

/*
 * Returns kth element of the exponential of U, results stored in user-supplied jet E
 */
real t_exp (series E, series U, int k);

/*
 * Selects either a trigonometric or hyperbolic version of the function
 */
typedef enum {TRIG, HYP} geometry;

/*
 * For returning combined recurrence values
 */
typedef struct {
    real a;
    real b;
} pair;

/*
 * Returns kth elements of both sine and cosine of U, results stored in user-supplied jets S and C
 */
pair t_sin_cos (series S, series C, series U, int k, geometry g);

/*
 * Returns kth elements of both tangent and squared secant of U, results stored in user-supplied jets T and S2
 */
pair t_tan_sec2 (series T, series S2, series U, int k, geometry g);

/*
 * For returning x, y, z values
 */
typedef struct {
    real x;
    real y;
    real z;
} components;

/*
 * Obligatory client method signatures
 */
void *get_p (int argc, char **argv, long order);

components ode (series x, series y, series z, void *params, int k);

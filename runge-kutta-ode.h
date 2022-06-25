/*
 * Interface for solving systems of Ordinary Differential Equations
 * using the Runge-Kutta-4 Method (RK4)
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include "ode-common.h"

/*
 * RK4 algorithm
 */
void rk4 (int dp, int n, real h, int steps, real x0, real y0, real z0, void *P, clock_t since);

/*
 * Obligatory client method signatures
 */

/*
 * Get a blob of parameter data from the model to be passed back in from ode()
 */
void *get_p (int argc, char **argv);

/*
 * Calculate velocity field
 */
components ode (real X, real Y, real Z, void *P);

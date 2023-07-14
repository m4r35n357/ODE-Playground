/*
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#pragma once

#include <stdbool.h>
#include "real.h"

/*
 * Retrieves integrator control parameters
 */
controls *symp_get_c (int argc, char **argv);

/*
 * Get model data from the command
 */
parameters *symp_init_p (int argc, char **argv);

/*
 * "Logarithmic" error function
 */
real error (real e);

/*
 * To pass a plotter as parameter
 */
typedef void (*plotter)(int dp, parameters *params, real t);

/*
 * To pass an integrator as parameter
 */
typedef void (*integrator)(controls *cont, parameters *params, real h);

/*
 * Coordinate updater dq = (dH/dp).dt
 */
void update_q (parameters *p, real c);

/*
 * Momentum updater dp = -(dH/dq).dt
 */
void update_p (parameters *p, real d);

/*
 * Call the symplectic solver
 */
void solve (char **argv, controls *cont, parameters *p, plotter output);

/*
 * Call the symplectic generator
 */
bool generate (controls *cont, parameters *p);

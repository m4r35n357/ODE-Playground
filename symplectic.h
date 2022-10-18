/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#pragma once

#include "real.h"

/*
 * Retrieves integrator control parameters
 */
controls *get_c_symp (char **argv);

/*
 * Get model data from the command to be passed into solve()
 */
void *get_p (int argc, char **argv);

/*
 * "Logarithmic" error function
 */
real error (real e);

/*
 * To pass a plotter as parameter
 */
typedef void (*plotter)(int dp, void *params, real t);

/*
 * To pass an integrator as parameter
 */
typedef void (*integrator)(controls *cont, void *params, real h);

/*
 * Coordinate updater dq = (dH/dp).dt
 */
void update_q (void *params, real c);

/*
 * Momentum updater dp = -(dH/dq).dt
 */
void update_p (void *params, real d);

/*
 * Call the symplectic solver
 */
void solve (char **argv, controls *cont, void *p, plotter output);

/*
 * Call the symplectic generator
 */
int generate (controls *cont, void *p);

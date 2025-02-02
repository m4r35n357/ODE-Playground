/*
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#pragma once
#include "real.h"

/*
 * Retrieves integrator control parameters
 */
controls *symp_get_c (int argc, char **argv);

/*
 * Get model data from the command
 */
model *symp_init_p (int argc, char **argv);

/*
 * "Logarithmic" error function
 */
real error (real e);

/*
 * To pass a plotter as parameter
 */
typedef void (*plotter)(int dp, model *p, real t);

/*
 * Coordinate updater dq = (dH/dp).dt
 */
void update_q (model *p, real c);

/*
 * Momentum updater dp = -(dH/dq).dt
 */
void update_p (model *p, real d);

/*
 * Call the symplectic solver
 */
void solve (controls *c, model *p, plotter output);

/*
 * Call the symplectic generator
 */
bool generate (controls *c, model *p);

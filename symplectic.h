
#include "real.h"

/*
 * Get a blob of model data from the command to be passed into solve()
 */
void *get_p (int argc, char **argv, int va_begin);

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
 * Call the integrator
 */
void solve (char **argv, controls *cont, void *p, plotter output);

int generate (controls *cont, void *p);

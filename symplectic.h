
#include "real.h"

/*
 * Bulk set model variables from the command line arguments (start onwards)
 */
void t_variables (char **argv, int start, int count, ...);

/*
 * Get a blob of parameter data from the model to be passed into solve()
 */
void *get_p (int argc, char **argv, int va_begin);

/*
 * "Logarithmic" error function
 */
real error (real e);

/*
 * To pass an updater as parameter
 */
typedef void (*updater)(void *params, real c_d);

/*
 * To pass a plotter as parameter
 */
typedef void (*plotter)(long dp, void *params, real t);

/*
 * To pass an integrator as parameter
 */
typedef void (*integrator)(void *params, updater uq, updater up, real h);

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
void solve (char **argv, void *p, updater uq, updater up, plotter output);

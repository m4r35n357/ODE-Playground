
#include "real.h"

typedef struct {
    long order;
    real step_size;
    long steps;
} controls;

/*
 * Get integrator control data from the command to be passed into solve()
 */
controls *get_c (char **argv);

/*
 * Bulk set model variables from the command line arguments (start onwards)
 */
void t_variables (char **argv, int start, int count, ...);

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
typedef void (*integrator)(void *params, real h);

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
void solve (char **argv, void *p, plotter output);

void *generate (controls *cont, void *p);

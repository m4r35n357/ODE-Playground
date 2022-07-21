
#include "real.h"

typedef struct Weights {
    real fwd, rev;
} weights;

typedef struct Controls {
    long order, step, steps;
    real step_size;
    weights r1, r2, r3, r4;
} controls;

/*
 * Get integrator control data from the command to be passed into solve()
 */
controls *get_c (char **argv);

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
void solve (char **argv, void *p, plotter output);

void *generate (controls *cont, void *p);

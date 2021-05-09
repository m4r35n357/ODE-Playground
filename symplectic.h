
#include "dual.h"

/*
 * Sets control parameters from the command line arguments (1 to 4)
 */
void t_stepper (char **argv, long *dp, long *method, real *h, long *nsteps);

/*
 * Bulk set model variables from the command line arguments (start onwards)
 */
void t_variables (char **argv, int start, int count, ...);

real error (real e);

typedef void (*updater)(void *, real);

typedef void (*base)(void *, updater, updater, real);

typedef void (*plotter)(long, void *, real);

typedef void (*integrator)(void *, updater, updater, real);

void update_q (void *params, real c);

void update_p (void *params, real d);

void solve (char **argv, void *p, updater uq, updater up, plotter output);


#include "real.h"

/*
 * Bulk set model variables from the command line arguments (start onwards)
 */
void t_variables (char **argv, int start, int count, ...);

real error (real e);

typedef void (*updater)(void *, real);

typedef void (*plotter)(long, void *, real);

typedef void (*integrator)(void *, updater, updater, real);

void update_q (void *params, real c);

void update_p (void *params, real d);

void solve (char **argv, void *p, updater uq, updater up, plotter output);

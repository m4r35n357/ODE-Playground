
#include <assert.h>
#include <math.h>

#define M_PI 3.14159265358979323846

typedef long double real;

/*
 * Sets control parameters from the command line arguments (1 to 4)
 */
void t_stepper (char **argv, long *dp, long *method, real *h, long *nsteps);

/*
 * Bulk set model variables from the command line arguments (5 onwards)
 */
void t_variables (char **argv, int count, ...);

real error (real e);

typedef void (*updater)(void *, real cd);

typedef void (*plotter)(long dp, void *, real h);

void solve (char **argv, void *p, updater uq, updater up, plotter output);

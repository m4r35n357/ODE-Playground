
#include <assert.h>
#include <math.h>

#define M_PI 3.14159265358979323846

typedef long double real;

/*
 * For returning dual numbers
 */
typedef struct {
    real val;
    real dot;
} dual;

dual d_dual (real a);

dual d_var (real a);

dual d_abs (dual a);

dual d_neg (dual b);

dual d_inv (dual a);

dual d_sqr (dual a);

dual d_shift (dual a, real b);

dual d_scale (dual a, real b);

dual d_add (dual a, dual b);

dual d_sub (dual a, dual b);

dual d_mul (dual a, dual b);

dual d_div (dual a, dual b);

dual d_exp (dual a);

dual d_log (dual a);

dual d_pow (dual a, real b);

dual d_sin (dual a);

dual d_cos (dual a);

dual d_tan (dual a);

dual d_sinh (dual a);

dual d_cosh (dual a);

dual d_tanh (dual a);

/*
 * Sets control parameters from the command line arguments (1 to 4)
 */
void t_stepper (char **argv, long *dp, long *method, real *h, long *nsteps);

/*
 * Bulk set model variables from the command line arguments (5 onwards)
 */
void t_variables (char **argv, int count, ...);

typedef void (*updater)(void *, real cd);

typedef void (*plotter)(long dp, void *, real h);

typedef void (*integrator)(void *, long size, real cd[], updater q, updater p);

void solve (char **argv, void *p, updater uq, updater up, plotter output);

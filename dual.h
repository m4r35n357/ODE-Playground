
#include <assert.h>
#include <math.h>
#include "real.h"

real get_PI (void);

real elevation_to_colatitude (real elevation);

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

dual d_sqrt (dual a);

dual d_pow (dual a, real b);

dual d_sin (dual a);

dual d_cos (dual a);

dual d_tan (dual a);

dual d_sinh (dual a);

dual d_cosh (dual a);

dual d_tanh (dual a);
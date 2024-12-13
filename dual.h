/*
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#pragma once
#include "real.h"

/*
 * For returning dual numbers
 */
typedef struct Dual { real val, dot; } dual;

dual d_dual (real a);

dual d_var (real a);

dual d_abs (const dual a);

dual d_inv (const dual a);

dual d_sqr (const dual a);

dual d_shift (const dual a, real b);

dual d_scale (const dual a, real b);

dual d_add (const dual a, dual b);

dual d_sub (const dual a, dual b);

dual d_mul (const dual a, dual b);

dual d_div (const dual a, dual b);

dual d_exp (const dual a);

dual d_ln (const dual a);

dual d_sqrt (const dual a);

dual d_pow (const dual a, real b);

dual d_sin (const dual a);

dual d_cos (const dual a);

dual d_tan (const dual a);

dual d_sinh (const dual a);

dual d_cosh (const dual a);

dual d_tanh (const dual a);

dual d_asin (const dual a);

dual d_acos (const dual a);

dual d_atan (const dual a);

dual d_asinh (const dual a);

dual d_acosh (const dual a);

dual d_atanh (const dual a);

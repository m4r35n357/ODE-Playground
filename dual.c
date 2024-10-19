/*
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dual.h"

dual d_dual (real a) {
    return (dual){.val = a, .dot = 0.0L};
}

dual d_var (real a) {
    return (dual){.val = a, .dot = 1.0L};
}

dual d_abs (const dual a) {
    CHECK(a.val != 0.0L);
    return (dual){.val = a.val < 0.0L ? - a.val : a.val, .dot = a.val < 0.0L ? - a.dot : a.dot};
}

dual d_inv (const dual b) {
    CHECK(b.val != 0.0L);
    return (dual){.val = 1.0L / b.val, .dot = - b.dot / SQR(b.val)};
}

dual d_sqr (const dual a) {
    return (dual){.val = SQR(a.val), .dot = 2.0L * a.val * a.dot};
}

dual d_shift (const dual a, real b) {
    return (dual){.val = a.val + b, .dot = a.dot};
}

dual d_scale (const dual a, real b) {
    return (dual){.val = a.val * b, .dot = a.dot * b};
}

dual d_add (const dual a, dual b) {
    return (dual){.val = a.val + b.val, .dot = a.dot + b.dot};
}

dual d_sub (const dual a, dual b) {
    return (dual){.val = a.val - b.val, .dot = a.dot - b.dot};
}

dual d_mul (const dual a, dual b) {
    return (dual){.val = a.val * b.val, .dot = a.dot * b.val + a.val * b.dot};
}

dual d_div (const dual a, dual b) {
    CHECK(b.val != 0.0L);
    return (dual){.val = a.val / b.val, .dot = (a.dot * b.val - a.val * b.dot) / SQR(b.val)};
}

dual d_exp (const dual a) {
    real exp = expl(a.val);
    return (dual){.val = exp, .dot = a.dot * exp};
}

dual d_ln (const dual a) {
    CHECK(a.val > 0.0L);
    return (dual){.val = logl(a.val), .dot = a.dot / a.val};
}

dual d_sqrt (const dual a) {
    CHECK(a.val > 0.0L);
    real sqrt = sqrtl(a.val);
    return (dual){.val = sqrt, .dot = a.dot * 0.5L / sqrt};
}

dual d_pow (const dual a, real b) {
    CHECK(a.val > 0.0L);
    real pow = powl(a.val, b);
    return (dual){.val = pow, .dot = a.dot * b * pow / a.val};
}

dual d_sin (const dual a) {
    return (dual){.val = sinl(a.val), .dot = a.dot * cosl(a.val)};
}

dual d_cos (const dual a) {
    return (dual){.val = cosl(a.val), .dot = - a.dot * sinl(a.val)};
}

dual d_tan (const dual a) {
    real tan = tanl(a.val);
    return (dual){.val = tan, .dot = a.dot * (1.0L + SQR(tan))};
}

dual d_sinh (const dual a) {
    return (dual){.val = sinhl(a.val), .dot = a.dot * coshl(a.val)};
}

dual d_cosh (const dual a) {
    return (dual){.val = coshl(a.val), .dot = a.dot * sinhl(a.val)};
}

dual d_tanh (const dual a) {
    real tanh = tanhl(a.val);
    return (dual){.val = tanh, .dot = a.dot * (1.0L - SQR(tanh))};
}

dual d_asin (const dual a) {
    CHECK(a.val > -1.0L && a.val < 1.0L);
    return (dual){.val = asinl(a.val), .dot = a.dot / sqrtl(1.0L - SQR(a.val))};
}

dual d_acos (const dual a) {
    CHECK(a.val > -1.0L && a.val < 1.0L);
    return (dual){.val = acosl(a.val), .dot = - a.dot / sqrtl(1.0L - SQR(a.val))};
}

dual d_atan (const dual a) {
    return (dual){.val = atanl(a.val), .dot = a.dot / (1.0L + SQR(a.val))};
}

dual d_asinh (const dual a) {
    return (dual){.val = asinhl(a.val), .dot = a.dot / sqrtl(SQR(a.val) + 1.0L)};
}

dual d_acosh (const dual a) {
    CHECK(a.val > 1.0L);
    return (dual){.val = acoshl(a.val), .dot = a.dot / sqrtl(SQR(a.val) - 1.0L)};
}

dual d_atanh (const dual a) {
    CHECK(a.val > -1.0L && a.val < 1.0L);
    return (dual){.val = atanhl(a.val), .dot = a.dot / (1.0L - SQR(a.val))};
}

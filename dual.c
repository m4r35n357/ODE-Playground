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

dual d_abs (dual a) {
    return (dual){.val = a.val < 0.0L ? - a.val : a.val, .dot = a.val < 0.0L ? - a.dot : a.dot};
}

dual d_inv (dual b) {
    CHECK(b.val != 0.0L);
    return (dual){.val = 1.0L / b.val, .dot = - b.dot / (b.val * b.val)};
}

dual d_sqr (dual a) {
    return (dual){.val = a.val * a.val, .dot = 2.0L * a.val * a.dot};
}

dual d_shift (dual a, real b) {
    return (dual){.val = a.val + b, .dot = a.dot};
}

dual d_scale (dual a, real b) {
    return (dual){.val = a.val * b, .dot = a.dot * b};
}

dual d_add (dual a, dual b) {
    return (dual){.val = a.val + b.val, .dot = a.dot + b.dot};
}

dual d_sub (dual a, dual b) {
    return (dual){.val = a.val - b.val, .dot = a.dot - b.dot};
}

dual d_mul (dual a, dual b) {
    return (dual){.val = a.val * b.val, .dot = a.val * b.dot + a.dot * b.val};
}

dual d_div (dual a, dual b) {
    CHECK(b.val != 0.0L);
    return (dual){.val = a.val / b.val, .dot = (a.dot * b.val - a.val * b.dot) / (b.val * b.val)};
}

dual d_exp (dual a) {
    real exp_val = expl(a.val);
    return (dual){.val = exp_val, .dot = a.dot * exp_val};
}

dual d_ln (dual a) {
    CHECK(a.val > 0.0L);
    return (dual){.val = logl(a.val), .dot = a.dot / a.val};
}

dual d_sqrt (dual a) {
    CHECK(a.val > 0.0L);
    real root_val = sqrtl(a.val);
    return (dual){.val = root_val, .dot = a.dot * 0.5L / root_val};
}

dual d_pow (dual a, real b) {
    CHECK(a.val > 0.0L);
    real pow_val = powl(a.val, b);
    return (dual){.val = pow_val, .dot = a.dot * b * pow_val / a.val};
}

dual d_sin (dual a) {
    return (dual){.val = sinl(a.val), .dot = a.dot * cosl(a.val)};
}

dual d_cos (dual a) {
    return (dual){.val = cosl(a.val), .dot = - a.dot * sinl(a.val)};
}

dual d_tan (dual a) {
    real tan_val = tanl(a.val);
    return (dual){.val = tan_val, .dot = a.dot * (1.0L + tan_val * tan_val)};
}

dual d_sinh (dual a) {
    return (dual){.val = sinhl(a.val), .dot = a.dot * coshl(a.val)};
}

dual d_cosh (dual a) {
    return (dual){.val = coshl(a.val), .dot = a.dot * sinhl(a.val)};
}

dual d_tanh (dual a) {
    real tanh_val = tanhl(a.val);
    return (dual){.val = tanh_val, .dot = a.dot * (1.0L - tanh_val * tanh_val)};
}

dual d_asin (dual a) {
    CHECK(a.val >= -1.0L && a.val <= 1.0L);
    return (dual){.val = asinl(a.val), .dot = a.dot / sqrtl(1.0L - a.val * a.val)};
}

dual d_acos (dual a) {
    CHECK(a.val >= -1.0L && a.val <= 1.0L);
    return (dual){.val = acosl(a.val), .dot = - a.dot / sqrtl(1.0L - a.val * a.val)};
}

dual d_atan (dual a) {
    return (dual){.val = atanl(a.val), .dot = a.dot / (1.0L + a.val * a.val)};
}

dual d_asinh (dual a) {
    return (dual){.val = asinhl(a.val), .dot = a.dot / sqrtl(1.0L + a.val * a.val)};
}

dual d_acosh (dual a) {
    CHECK(a.val >= 1.0L);
    return (dual){.val = acoshl(a.val), .dot = a.dot / sqrtl(1.0L - a.val * a.val)};
}

dual d_atanh (dual a) {
    CHECK(a.val >= -1.0L && a.val <= 1.0L);
    return (dual){.val = atanhl(a.val), .dot = a.dot / (1.0L - a.val * a.val)};
}

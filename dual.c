/*
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <math.h>
#include "dual.h"

dual d_dual (real a) {
    return (dual){a, 0.0L};
}

dual d_var (real a) {
    return (dual){a, 1.0L};
}

dual d_abs (dual a) {
    return (dual){a.val < 0.0L ? - a.val : a.val, a.val < 0.0L ? - a.dot : a.dot};
}

dual d_inv (dual b) {
    assert(b.val != 0.0L);
    return (dual){1.0L / b.val, - b.dot / (b.val * b.val)};
}

dual d_sqr (dual a) {
    return (dual){a.val * a.val, 2.0L * a.val * a.dot};
}

dual d_shift (dual a, real b) {
    return (dual){a.val + b, a.dot};
}

dual d_scale (dual a, real b) {
    return (dual){a.val * b, a.dot * b};
}

dual d_add (dual a, dual b) {
    return (dual){a.val + b.val, a.dot + b.dot};
}

dual d_sub (dual a, dual b) {
    return (dual){a.val - b.val, a.dot - b.dot};
}

dual d_mul (dual a, dual b) {
    return (dual){a.val * b.val, a.val * b.dot + a.dot * b.val};
}

dual d_div (dual a, dual b) {
    assert(b.val != 0.0L);
    return (dual){a.val / b.val, (a.dot * b.val - a.val * b.dot) / (b.val * b.val)};
}

dual d_exp (dual a) {
    real exp_val = expl(a.val);
    return (dual){exp_val, a.dot * exp_val};
}

dual d_ln (dual a) {
    assert(a.val > 0.0L);
    return (dual){logl(a.val), a.dot / a.val};
}

dual d_sqrt (dual a) {
    assert(a.val > 0.0L);
    real root_val = sqrtl(a.val);
    return (dual){root_val, a.dot * 0.5L / root_val};
}

dual d_pow (dual a, real b) {
    assert(a.val > 0.0L);
    real pow_val = powl(a.val, b);
    return (dual){pow_val, a.dot * b * pow_val / a.val};
}

dual d_sin (dual a) {
    return (dual){sinl(a.val), a.dot * cosl(a.val)};
}

dual d_cos (dual a) {
    return (dual){cosl(a.val), - a.dot * sinl(a.val)};
}

dual d_tan (dual a) {
    real tan_val = tanl(a.val);
    return (dual){tan_val, a.dot * (1.0L + tan_val * tan_val)};
}

dual d_sinh (dual a) {
    return (dual){sinhl(a.val), a.dot * coshl(a.val)};
}

dual d_cosh (dual a) {
    return (dual){coshl(a.val), a.dot * sinhl(a.val)};
}

dual d_tanh (dual a) {
    real tanh_val = tanhl(a.val);
    return (dual){tanh_val, a.dot * (1.0L - tanh_val * tanh_val)};
}

dual d_asin (dual a) {
    assert(a.val >= -1.0L && a.val <= 1.0L);
    return (dual){asinl(a.val), a.dot / sqrtl(1.0L - a.val * a.val)};
}

dual d_acos (dual a) {
    assert(a.val >= -1.0L && a.val <= 1.0L);
    return (dual){acosl(a.val), - a.dot / sqrtl(1.0L - a.val * a.val)};
}

dual d_atan (dual a) {
    return (dual){atanl(a.val), a.dot / (1.0L + a.val * a.val)};
}

dual d_asinh (dual a) {
    return (dual){asinhl(a.val), a.dot / sqrtl(1.0L + a.val * a.val)};
}

dual d_acosh (dual a) {
    assert(a.val >= 1.0L);
    return (dual){acoshl(a.val), a.dot / sqrtl(1.0L - a.val * a.val)};
}

dual d_atanh (dual a) {
    assert(a.val >= -1.0L && a.val <= 1.0L);
    return (dual){atanhl(a.val), a.dot / (1.0L - a.val * a.val)};
}

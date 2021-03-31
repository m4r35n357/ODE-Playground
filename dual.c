
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include "dual.h"

static const real MY_PI = 3.1415926535897932384626433832795029L;

real get_PI (void) {
    return MY_PI;
}

real elevation_to_colatitude (real elevation) {
    return (90.0L - elevation) * MY_PI  / 180.0L;
}

dual d_dual (real a) {
    return (dual) { .val = a, .dot = 0.0L };
}

dual d_var (real a) {
    return (dual) { .val = a, .dot = 1.0L };
}

dual d_abs (dual a) {
    return (dual) { .val = a.val < 0.0L ? - a.val : a.val, .dot = a.val < 0.0L ? - a.dot : a.dot };
}

dual d_neg (dual b) {
    return (dual) { .val = - b.val, .dot = - b.dot };
}

dual d_inv (dual b) {
    assert(b.val != 0.0L);
    return (dual) { .val = 1.0L / b.val, .dot = - b.dot / (b.val * b.val) };
}

dual d_sqr (dual a) {
    return (dual) { .val = a.val * a.val, .dot = 2.0L * a.val * a.dot };
}

dual d_shift (dual a, real b) {
    return (dual) { .val = a.val + b, .dot = a.dot };
}

dual d_scale (dual a, real b) {
    return (dual) { .val = a.val * b, .dot = a.dot * b };
}

dual d_add (dual a, dual b) {
    return (dual) { .val = a.val + b.val, .dot = a.dot + b.dot };
}

dual d_sub (dual a, dual b) {
    return (dual) { .val = a.val - b.val, .dot = a.dot - b.dot };
}

dual d_mul (dual a, dual b) {
    return (dual) { .val = a.val * b.val, .dot = a.val * b.dot + a.dot * b.val };
}

dual d_div (dual a, dual b) {
    assert(b.val != 0.0L);
    return (dual) { .val = a.val / b.val, .dot = (a.dot * b.val - a.val * b.dot) / (b.val * b.val) };
}

dual d_exp (dual a) {
    real exp_val = expl(a.val);
    return (dual) { .val = exp_val, .dot = a.dot * exp_val };
}

dual d_log (dual a) {
    assert(a.val > 0.0L);
    return (dual) { .val = logl(a.val), .dot = a.dot / a.val };
}

dual d_pow (dual a, real b) {
    assert(a.val > 0.0L);
    return (dual) { .val = powl(a.val, b), .dot = a.dot * b * powl(a.val, (b - 1.0L)) };
}

dual d_sin (dual a) {
    return (dual) { .val = sinl(a.val), .dot = a.dot * cosl(a.val) };
}

dual d_cos (dual a) {
    return (dual) { .val = cosl(a.val), .dot = - a.dot * sinl(a.val) };
}

dual d_tan (dual a) {
    real tan_val = tanl(a.val);
    return (dual) { .val = tan_val, .dot = a.dot * (1.0L + tan_val * tan_val) };
}

dual d_sinh (dual a) {
    return (dual) { .val = sinhl(a.val), .dot = a.dot * coshl(a.val) };
}

dual d_cosh (dual a) {
    return (dual) { .val = coshl(a.val), .dot = a.dot * sinhl(a.val) };
}

dual d_tanh (dual a) {
    real tanh_val = tanhl(a.val);
    return (dual) { .val = tanh_val, .dot = a.dot * (1.0L - tanh_val * tanh_val) };
}

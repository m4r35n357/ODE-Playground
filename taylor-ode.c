/*
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include "taylor-ode.h"

void t_xyz_output (real x, real y, real z, real t) {
    printf("%+.12Le %+.12Le %+.12Le %+.6Le\n", x, y, z, t);
}

void t_stepper (char **argv, long *n, real *h, long *nsteps) {
    *n = strtol(argv[3], NULL, 10);
    *h = strtold(argv[4], NULL);
    *nsteps = strtol(argv[5], NULL, 10);
}

void t_args (char **argv, int count, ...) {
    va_list vars;
    va_start(vars, count);
    for (int i = 6; i < count; i++) {
        *va_arg(vars, real *) = strtold(argv[i], NULL);
    }
    va_end(vars);
}

series t_jet (int n) {
    assert(n > 0);
    return malloc(sizeof (real) * n);
}

series t_jet_c (int n, real value) {
    series jet = t_jet(n);
    jet[0] = value;
    for (int i = 1; i < n; i++) {
        jet[i] = 0.0;
    }
    return jet;
}

real t_horner (series jet, int n, real h) {
    assert(n > 0);
    real sum = 0.0;
    for (int i = n; i >= 0; i--) {
        sum = sum * h + jet[i];
    }
    if (isnan(sum) || isinf(sum)) { fprintf(stderr, "OVERFLOW !\n"); exit(1); }
    return jet[0] = sum;
}

real t_abs (series u, int k) {
    assert(k >= 0);
    return u[0] < 0.0 ? - u[k] : u[k];
}

static real cauchy (series a, series b, int k, int lower, int upper) {
    real c = 0.0;
    for (int j = lower; j <= upper; j++) {
        c += a[j] * b[k - j];
    }
    return c;
}

real t_prod (series u, series v, int k) {
    assert(k >= 0);
    return cauchy(u, v, k, 0, k);
}

real t_quot (series q, series u, series v, int k) {
    assert(v[0] != 0.0);
    assert(q != u && q != v && u != v);
    assert(k >= 0);
    return q[k] = (k == 0 ? u[0] : u[k] - cauchy(q, v, k, 0, k - 1)) / v[0];
}

real t_inv (series i, series v, int k) {
    assert(v[0] != 0.0);
    assert(i != v);
    assert(k >= 0);
    return i[k] = (k == 0 ? 1.0 : - cauchy(i, v, k, 0, k - 1)) / v[0];
}

real t_sqr (series u, int k) {
    assert(k >= 0);
    return cauchy(u, u, k, 0, k);
}

real t_sqrt (series r, series u, int k) {
    assert(u[0] > 0.0);
    assert(r != u);
    assert(k >= 0);
    return r[k] = k == 0 ? sqrt(u[0]) : 0.5 * (u[k] - cauchy(r, r, k, 1, k - 1)) / r[0];
}

static real d_cauchy (series h, series u, int k, int lower, int upper, real factor) {
    real f = 0.0;
    for (int j = lower; j <= upper; j++) {
        f += h[j] * (k - j) * u[k - j];
    }
    return factor * f / k;
}

real t_exp (series e, series u, int k) {
    assert(e != u);
    assert(k >= 0);
    return e[k] = k == 0 ? exp(u[0]) : d_cauchy(e, u, k, 0, k - 1, 1.0);
}

pair t_sin_cos (series s, series c, series u, int k, geometry g) {
    assert(s != c && s != u && c != u);
    assert(k >= 0);
    if (k == 0) {
        return (pair){s[k] = g == TRIG ? sin(u[0]) : sinh(u[0]), c[k] = g == TRIG ? cos(u[0]) : cosh(u[0])};
    } else {
        return (pair){s[k] = d_cauchy(c, u, k, 0, k - 1, 1.0), c[k] = d_cauchy(s, u, k, 0, k - 1, g == TRIG ? - 1.0 : 1.0)};
    }
}

pair t_tan_sec2 (series t, series s, series u, int k, geometry g) {
    assert(t != s && t != u && s != u);
    assert(k >= 0);
    if (k == 0) {
        return (pair){t[k] = g == TRIG ? tan(u[0]) : tanh(u[0]), s[k] = g == TRIG ? 1.0 + t[0] * t[0] : 1.0 - t[0] * t[0]};
    } else {
        return (pair){t[k] = d_cauchy(s, u, k, 0, k - 1, 1.0), s[k] = d_cauchy(t, t, k, 0, k - 1, g == TRIG ? 2.0 : - 2.0)};
    }
}

real t_pwr (series p, series u, real a, int k) {
    assert(u[0] > 0.0);
    assert(p != u);
    assert(k >= 0);
    return p[k] = k == 0 ? pow(u[0], a) : (d_cauchy(p, u, k, 0, k - 1, a) - d_cauchy(u, p, k, 1, k - 1, 1.0)) / u[0];
}

real t_ln (series l, series u, int k) {
    assert(u[0] > 0.0);
    assert(l != u);
    assert(k >= 0);
    return l[k] = k == 0 ? log(u[0]) : (u[k] - d_cauchy(u, l, k, 1, k - 1, 1.0)) / u[0];
}

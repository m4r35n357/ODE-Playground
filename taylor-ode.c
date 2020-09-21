/*
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include "taylor-ode.h"

void t_output (real x, real y, real z, real t) {
    printf("%+.12Le %+.12Le %+.12Le %+.6Le\n", x, y, z, t);
}

void t_stepper (char **argv, long *n, real *h, long *nsteps, real *x, real *y, real *z) {
    *n = strtol(argv[3], NULL, BASE);
    assert(*n > 0);
    *h = strtold(argv[4], NULL);
    assert(*h > 0.0);
    *nsteps = strtol(argv[5], NULL, BASE);
    assert(*nsteps > 0);
    *x = strtold(argv[6], NULL);
    *y = strtold(argv[7], NULL);
    *z = strtold(argv[8], NULL);
}

void t_args (char **argv, int count, ...) {
    va_list vars;
    va_start(vars, count);
    for (int i = 9; i < count; i++) {
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
    if (isnan(sum) || isinf(sum)) {
        fprintf(stderr, "OVERFLOW !\n");
        exit(1);
    }
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
    if (k == 0) return (pair) {
        .a = s[k] = g == TRIG ? sin(u[0]) : sinh(u[0]),
        .b = c[k] = g == TRIG ? cos(u[0]) : cosh(u[0])
    };
    return (pair) {
        .a = s[k] = d_cauchy(c, u, k, 0, k - 1, 1.0),
        .b = c[k] = d_cauchy(s, u, k, 0, k - 1, g == TRIG ? - 1.0 : 1.0)
    };
}

pair t_tan_sec2 (series t, series s, series u, int k, geometry g) {
    assert(t != s && t != u && s != u);
    assert(k >= 0);
    if (k == 0) return (pair) {
        .a = t[k] = g == TRIG ? tan(u[0]) : tanh(u[0]),
        .b = s[k] = g == TRIG ? 1.0 + t[0] * t[0] : 1.0 - t[0] * t[0]
    };
    return (pair) {
        .a = t[k] = d_cauchy(s, u, k, 0, k - 1, 1.0),
        .b = s[k] = d_cauchy(t, t, k, 0, k - 1, g == TRIG ? 2.0 : - 2.0)
    };
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

void tsm (int argc, char **argv, tsm_model ode, ode_params get_p, ode_inters get_i) {
    long n, steps;
    real x0, y0, z0, h;
    t_stepper(argv, &n, &h, &steps, &x0, &y0, &z0);
    series x = t_jet_c(n + 1, x0), y = t_jet_c(n + 1, y0), z = t_jet_c(n + 1, z0);
    void *p = get_p(argc, argv, n);
    void *i = get_i == NULL ? NULL : get_i(n);
    t_output(x[0], y[0], z[0], 0.0);
    for (long step = 1; step < steps + 1; step++) {
        for (int k = 0; k < n; k++) {
            components c = ode(x, y, z, p, i, k);
            x[k + 1] = c.x / (k + 1);
            y[k + 1] = c.y / (k + 1);
            z[k + 1] = c.z / (k + 1);
        }
        t_output(t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h), h * step);
    }
}

void rk4 (int argc, char **argv, rk4_model ode, ode_params get_p) {
    long interval, steps;
    real x, y, z, h;
    t_stepper(argv, &interval, &h, &steps, &x, &y, &z);
    void *p = get_p(argc, argv, 0);
    t_output(x, y, z, 0.0);
    for (long step = 1; step < steps + 1; step++) {
        components k1 = ode(x, y, z, p);
        components k2 = ode(x + 0.5 * k1.x * h, y + 0.5 * k1.y * h, z + 0.5 * k1.z * h, p);
        components k3 = ode(x + 0.5 * k2.x * h, y + 0.5 * k2.y * h, z + 0.5 * k2.z * h, p);
        components k4 = ode(x + k3.x * h, y + k3.y * h, z + k3.z * h, p);
        x += h * (k1.x + 2.0 * (k2.x + k3.x) + k4.x) / 6.0;
        y += h * (k1.y + 2.0 * (k2.y + k3.y) + k4.y) / 6.0;
        z += h * (k1.z + 2.0 * (k2.z + k3.z) + k4.z) / 6.0;
        if (step % interval == 0) {
            t_output(x, y, z, h * step);
        }
    }
}

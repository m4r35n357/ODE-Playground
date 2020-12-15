/*
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include "taylor-ode.h"

void t_output (long dp, real x, real y, real z, real t, char *x_label, char *y_label, char *z_label) {
    char fs[128];
    sprintf(fs, "%%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.6Le %s %s %s\n", dp, dp, dp, x_label, y_label, z_label);
    printf(fs, x, y, z, t);
}

void t_control (char **argv, long *dp, long *n, real *h, long *nsteps, real *x, real *y, real *z) {
    *dp = strtol(argv[1], NULL, 10);
    *n = strtol(argv[2], NULL, 10);
    assert(*n >= 2 && *n <= 34);
    *h = strtold(argv[3], NULL);
    assert(*h > 0.0L && *h <= 1.0L);
    *nsteps = strtol(argv[4], NULL, 10);
    assert(*nsteps >= 1 && *nsteps <= 100000);
    *x = strtold(argv[5], NULL);
    *y = strtold(argv[6], NULL);
    *z = strtold(argv[7], NULL);
}

void t_params (char **argv, int argc, ...) {
    va_list vars;
    va_start(vars, argc);
    for (int i = 8; i < argc; i++) {
        *va_arg(vars, real *) = strtold(argv[i], NULL);
    }
    va_end(vars);
}

series t_jet (long n) {
    return calloc((size_t)n, sizeof (real));
}

series t_jet_c (long n, real value) {
    series jet = t_jet(n);
    jet[0] = value;
    return jet;
}

real t_horner (series jet, long n, real h) {
    real sum = 0.0L;
    for (long i = n; i >= 0; i--) {
        sum = sum * h + jet[i];
    }
    if (isnan(sum) || isinf(sum)) exit(1);
    return jet[0] = sum;
}

real t_abs (series u, int k) {
    return u[0] < 0.0L ? - u[k] : u[k];
}

real t_prod (series a, series b, int k) {
    real p = 0.0L;
    for (int j = 0; j <= k; j++) {
        p += a[j] * b[k - j];
    }
    return p;
}

static real cauchy (series h, series u, int k, real factor) {
    real f = 0.0L;
    for (int j = 0; j < k; j++) {
        f += h[j] * (k - j) * u[k - j];
    }
    return factor * f / k;
}

real t_exp (series e, series u, int k) {
    assert(e != u);
    return e[k] = k == 0 ? expl(u[0]) : cauchy(e, u, k, 1.0L);
}

pair t_sin_cos (series s, series c, series u, int k, geometry g) {
    assert(s != c && s != u && c != u);
    if (k == 0) return (pair) {
        .a = s[k] = g == TRIG ? sinl(u[0]) : sinhl(u[0]),
        .b = c[k] = g == TRIG ? cosl(u[0]) : coshl(u[0])
    };
    return (pair) {
        .a = s[k] = cauchy(c, u, k, 1.0L),
        .b = c[k] = cauchy(s, u, k, g == TRIG ? - 1.0L : 1.0L)
    };
}

pair t_tan_sec2 (series t, series s, series u, int k, geometry g) {
    assert(t != s && t != u && s != u);
    if (k == 0) return (pair) {
        .a = t[k] = g == TRIG ? tanl(u[0]) : tanhl(u[0]),
        .b = s[k] = g == TRIG ? 1.0L + t[0] * t[0] : 1.0L - t[0] * t[0]
    };
    return (pair) {
        .a = t[k] = cauchy(s, u, k, 1.0L),
        .b = s[k] = cauchy(t, t, k, g == TRIG ? 2.0L : - 2.0L)
    };
}

void tsm (int argc, char **argv, tsm_model ode, tsm_params get_p, tsm_inters get_i) {
    long n, steps, dp;
    real x0, y0, z0, h;
    t_control(argv, &dp, &n, &h, &steps, &x0, &y0, &z0);
    series x = t_jet_c(n + 1, x0), y = t_jet_c(n + 1, y0), z = t_jet_c(n + 1, z0);
    void *p = get_p(argc, argv, n);
    void *i = get_i == NULL ? NULL : get_i(n);
    t_output(dp, x[0], y[0], z[0], 0.0L, "_", "_", "_");
    for (long step = 1; step <= steps; step++) {
        real xdot = x[1], ydot = y[1], zdot = z[1];
        for (int k = 0; k < n; k++) {
            components c = ode(x, y, z, p, i, k);
            x[k + 1] = c.x / (k + 1);
            y[k + 1] = c.y / (k + 1);
            z[k + 1] = c.z / (k + 1);
        }
        t_output(dp, t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h), h * step,
                 x[1] * xdot < 0.0L ? (x[2] > 0.0L ? "x" : (x[2] < 0.0L ? "X" : "_")) : "_",
                 y[1] * ydot < 0.0L ? (y[2] > 0.0L ? "y" : (y[2] < 0.0L ? "Y" : "_")) : "_",
                 z[1] * zdot < 0.0L ? (z[2] > 0.0L ? "z" : (z[2] < 0.0L ? "Z" : "_")) : "_");
    }
}

/*
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include "taylor-ode.h"

const int BASE = 10;

void t_params (char **argv, int argc, ...) {
    va_list model_params;
    va_start(model_params, argc);
    for (int i = 8; i < argc; i++) {
        *va_arg(model_params, real *) = strtold(argv[i], NULL);
    }
    va_end(model_params);
}

series t_jet (long n) {
    series s = calloc((size_t)n, sizeof (real));
    if (s == NULL) {
        fprintf(stderr, "Allocation failure!\n");
        exit(1);
    }
    for (long i = 0; i < n; i++) {
        s[i] = 0.0L;
    }
    return s;
}

void t_output (long dp, real x, real y, real z, real t, char *x_label, char *y_label, char *z_label) {
    char fs[128];
    sprintf(fs, "%%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.6Le %s %s %s\n", dp, dp, dp, x_label, y_label, z_label);
    printf(fs, x, y, z, t);
}

real t_horner (series jet, long n, real h) {
    real sum = 0.0L;
    for (long i = n; i >= 0; i--) {
        sum = sum * h + jet[i];
    }
    if (isnan(sum) || isinf(sum)) {
        fprintf(stderr, "Value error!\n");
        exit(2);
    }
    return sum;
}

void tsm (int argc, char **argv, long dp, long n, real h, long steps, real x0, real y0, real z0) {
    series x = t_jet(n + 1); x[0] = x0;
    series y = t_jet(n + 1); y[0] = y0;
    series z = t_jet(n + 1); z[0] = z0;
    void *p = get_p(argc, argv, n);
    components slope = (components) {0.0L, 0.0L, 0.0L};
    for (long step = 0; step < steps; step++) {
        for (int k = 0; k < n; k++) {
            components v = ode(x, y, z, p, k);
            x[k + 1] = v.x / (k + 1);
            y[k + 1] = v.y / (k + 1);
            z[k + 1] = v.z / (k + 1);
        }
        t_output(dp, x[0], y[0], z[0], h * step, x[1] * slope.x < 0.0L ? (x[2] > 0.0L ? "x" : "X") : "_",
                                                 y[1] * slope.y < 0.0L ? (y[2] > 0.0L ? "y" : "Y") : "_",
                                                 z[1] * slope.z < 0.0L ? (z[2] > 0.0L ? "z" : "Z") : "_");
        slope = (components) {x[1], y[1], z[1]};
        x[0] = t_horner(x, n, h);
        y[0] = t_horner(y, n, h);
        z[0] = t_horner(z, n, h);
    }
    t_output(dp, x[0], y[0], z[0], h * steps, "_", "_", "_");
}

real t_const (real a, int k) {
    return k == 0 ? a : 0.0L;
}

real t_abs (series u, int k) {
    return u[0] < 0.0L ? - u[k] : u[k];
}

static real cauchy (series a, series b, int k, int j_lower, int j_upper) {
    real sum = 0.0L;
    for (int j = j_lower; j <= j_upper; j++) {
        sum += a[j] * b[k - j];
    }
    return sum;
}

real t_prod (series u, series v, int k) {
    return cauchy(u, v, k, 0, k);
}

real t_sqr (series u, int k) {
    return cauchy(u, u, k, 0, k);
}

real t_quot (series q, series u, series v, int k) {
    assert(v[0] != 0.0L);
    assert(q != u && q != v);
    return q[k] = (k == 0 ? u[0] : u[k] - cauchy(q, v, k, 0, k - 1)) / v[0];
}

real t_inv (series i, series v, int k) {
    assert(v[0] != 0.0L);
    assert(i != v);
    return i[k] = (k == 0 ? 1.0L : - cauchy(i, v, k, 0, k - 1)) / v[0];
}

real t_sqrt (series r, series u, int k) {
    assert(u[0] > 0.0L);
    assert(r != u);
    return r[k] = k == 0 ? sqrtl(u[0]) : 0.5L * (u[k] - cauchy(r, r, k, 1, k - 1)) / r[0];
}

static real f_k (series df_du, series u, int k, int j_lower, int j_upper) {
    real sum = 0.0L;
    for (int j = j_lower; j <= j_upper; j++) {
        sum += df_du[j] * (k - j) * u[k - j];
    }
    return sum / k;
}

real t_exp (series e, series u, int k) {
    assert(e != u);
    return e[k] = k == 0 ? expl(u[0]) : f_k(e, u, k, 0, k - 1);
}

pair t_sin_cos (series s, series c, series u, int k, geometry g) {
    assert(s != c && s != u && c != u);
    return k == 0 ? (pair) {
        .a = s[k] = g == TRIG ? sinl(u[0]) : sinhl(u[0]),
        .b = c[k] = g == TRIG ? cosl(u[0]) : coshl(u[0])
    } : (pair) {
        .a = s[k] = f_k(c, u, k, 0, k - 1),
        .b = c[k] = f_k(s, u, k, 0, k - 1) * (g == TRIG ? - 1.0L : 1.0L)
    };
}

pair t_tan_sec2 (series t, series s, series u, int k, geometry g) {
    assert(t != s && t != u && s != u);
    return k == 0 ? (pair) {
        .a = t[k] = g == TRIG ? tanl(u[0]) : tanhl(u[0]),
        .b = s[k] = g == TRIG ? 1.0L + t[0] * t[0] : 1.0L - t[0] * t[0]
    } : (pair) {
        .a = t[k] = f_k(s, u, k, 0, k - 1),
        .b = s[k] = f_k(t, t, k, 0, k - 1) * (g == TRIG ? 2.0L : - 2.0L)
    };
}

real t_pwr (series p, series u, real a, int k) {
    assert(u[0] > 0.0L);
    assert(p != u);
    return p[k] = k == 0 ? powl(u[0], a) : (f_k(p, u, k, 0, k - 1) * a - f_k(u, p, k, 1, k - 1)) / u[0];
}

real t_ln (series l, series u, int k) {
    assert(u[0] > 0.0L);
    assert(l != u);
    return l[k] = k == 0 ? logl(u[0]) : (u[k] - f_k(u, l, k, 1, k - 1)) / u[0];
}

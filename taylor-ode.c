/*
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
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

void t_params (char **argv, int argc, ...) {
    va_list vars;
    va_start(vars, argc);
    for (int i = 8; i < argc; i++) {
        *va_arg(vars, real *) = strtold(argv[i], NULL);
    }
    va_end(vars);
}

series t_jet (long n) {
    series s = calloc((size_t)n, sizeof (real));
    if (s == NULL) exit(2);
    for (long i = 0; i < n; i++) {
        s[i] = 0.0L;
    }
    return s;
}

real t_horner (series jet, long n, real h) {
    real sum = 0.0L;
    for (long i = n; i >= 0; i--) {
        sum = sum * h + jet[i];
    }
    if (isnan(sum) || isinf(sum)) exit(1);
    return jet[0] = sum;
}

real t_const (real a, int k) {
	return k == 0 ? a : 0.0L;
}

real t_abs (series u, int k) {
    return u[0] < 0.0L ? - u[k] : u[k];
}

static real cauchy (series a, series b, int k, int lower, int upper) {
    real sum = 0.0L;
    for (int j = lower; j <= upper; j++) {
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

static real f_k (series df_du, series u, int k, int lower, int upper, real factor) {
    real sum = 0.0L;
    for (int j = lower; j <= upper; j++) {
        sum += df_du[j] * (k - j) * u[k - j];
    }
    return factor * sum / k;
}

real t_exp (series e, series u, int k) {
    assert(e != u);
    return e[k] = k == 0 ? expl(u[0]) : f_k(e, u, k, 0, k - 1, 1.0L);
}

pair t_sin_cos (series s, series c, series u, int k, geometry g) {
    assert(s != c && s != u && c != u);
    if (k == 0) return (pair) {
        .a = s[k] = g == TRIG ? sinl(u[0]) : sinhl(u[0]),
        .b = c[k] = g == TRIG ? cosl(u[0]) : coshl(u[0])
    };
    return (pair) {
        .a = s[k] = f_k(c, u, k, 0, k - 1, 1.0L),
        .b = c[k] = f_k(s, u, k, 0, k - 1, g == TRIG ? - 1.0L : 1.0L)
    };
}

pair t_tan_sec2 (series t, series s, series u, int k, geometry g) {
    assert(t != s && t != u && s != u);
    if (k == 0) return (pair) {
        .a = t[k] = g == TRIG ? tanl(u[0]) : tanhl(u[0]),
        .b = s[k] = g == TRIG ? 1.0L + t[0] * t[0] : 1.0L - t[0] * t[0]
    };
    return (pair) {
        .a = t[k] = f_k(s, u, k, 0, k - 1, 1.0L),
        .b = s[k] = f_k(t, t, k, 0, k - 1, g == TRIG ? 2.0L : - 2.0L)
    };
}
real t_pwr (series p, series u, real a, int k) {
    assert(u[0] > 0.0L);
    assert(p != u);
    return p[k] = k == 0 ? powl(u[0], a) : (f_k(p, u, k, 0, k - 1, a) - f_k(u, p, k, 1, k - 1, 1.0L)) / u[0];
}

real t_ln (series l, series u, int k) {
    assert(u[0] > 0.0L);
    assert(l != u);
    return l[k] = k == 0 ? logl(u[0]) : (u[k] - f_k(u, l, k, 1, k - 1, 1.0L)) / u[0];
}

/*
 * Arbitrary-Order Taylor Series Integrator, together with a well-tested set of recurrence relations
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "taylor-ode.h"

controls *get_c_tsm (int argc, char **argv) {
    PRINT_ARGS(argc, argv);
    controls *c = malloc(sizeof (controls)); CHECK(c);
    c->order = (int)strtol(argv[2], NULL, BASE); CHECK(c->order >= 2 && c->order <= 64);
    c->step_size = strtold(argv[3], NULL);       CHECK(c->step_size > 0.0L);
    c->steps = (int)strtol(argv[4], NULL, BASE); CHECK(c->steps >= 0 && c->steps <= 1000000);
    return c;
}

series3 *initial_values (char **argv, int order) {
    series3 *jets = malloc(sizeof (series3)); CHECK(jets);
    jets->x = t_const(order + 1, strtold(argv[5], NULL));
    jets->y = t_const(order + 1, strtold(argv[6], NULL));
    jets->z = t_const(order + 1, strtold(argv[7], NULL));
    return jets;
}

void t_params (char **argv, int argc, ...) {
    va_list model;
    va_start(model, argc);
    for (int i = 8; i < argc; i++) {
        *va_arg(model, real *) = strtold(argv[i], NULL);
    }
    va_end(model);
}

void t_out (int dp, real x, real y, real z, real t, char *x_tag, char *y_tag, char *z_tag, clock_t since) {
    real cpu = (real)(clock() - since) / CLOCKS_PER_SEC;
    if (!dp) {
        printf("%+La %+La %+La %.6Le %s %s %s %.3Lf\n", x, y, z, t, x_tag, y_tag, z_tag, cpu);
    } else {
        printf("%+.*Le %+.*Le %+.*Le %.6Le %s %s %s %.3Lf\n", dp, x, dp, y, dp, z, t, x_tag, y_tag, z_tag, cpu);
    }
}

series t_jet (int n) {
    CHECK(n > 0);
    series s = malloc((size_t)n * sizeof (real)); CHECK(s);
    return s;
}

series t_const (int n, real a) {
    series c = t_jet(n);
    for (int k = 0; k < n; k++) {
        c[k] = !k ? a : 0.0L;
    }
    return c;
}

real t_horner (series s, int n, real h) {
    real sum = 0.0L;
    for (int i = n; i >= 0; i--) {
        sum = sum * h + s[i];
    }
    CHECK(isfinite(sum));
    return sum;
}

static char *tag (series jet, real slope, char *min, char *max) {
    return jet[1] * slope < 0.0L ? (jet[2] > 0.0L ? min : max) : "_";
}

static void derivatives (series3 *j, void *p, int n) {
    for (int k = 0; k < n; k++) {
        triplet v = ode(j->x, j->y, j->z, p, k);
        j->x[k + 1] = v.x / (k + 1);
        j->y[k + 1] = v.y / (k + 1);
        j->z[k + 1] = v.z / (k + 1);
    }
}

static void update (series3 *j, int n, real h) {
    j->x[0] = t_horner(j->x, n, h);
    j->y[0] = t_horner(j->y, n, h);
    j->z[0] = t_horner(j->z, n, h);
}

void tsm_stdout (int dp, controls *c, series3 *jets, void *p, clock_t t0) {
    triplet slope = {0.0L, 0.0L, 0.0L};
    for (int step = 0; step < c->steps; step++) {
        derivatives(jets, p, c->order);
        t_out(dp, jets->x[0], jets->y[0], jets->z[0], c->step_size * step,
              tag(jets->x, slope.x, "x", "X"), tag(jets->y, slope.y, "y", "Y"), tag(jets->z, slope.z, "z", "Z"), t0);
        slope.x = jets->x[1]; slope.y = jets->y[1]; slope.z = jets->z[1];
        update(jets, c->order, c->step_size);
    }
    t_out(dp, jets->x[0], jets->y[0], jets->z[0], c->step_size * c->steps, "_", "_", "_", t0);
}

bool tsm_gen (controls *c, series3 *jets, void *p) {
    static bool looping = false;
    if (looping) goto resume; else looping = true;
    for (c->step = 0; c->step < c->steps; c->step++) {
        derivatives(jets, p, c->order);
        update(jets, c->order, c->step_size);
        return true;
        resume: ;
    }
    return looping = false;
}

real t_abs (series u, int k) {
    return u[0] < 0.0L ? -u[k] : u[k];
}

static real fa (series a, series b, int k0, int k1, int k) {
    real _ = 0.0L;
    for (int j = k0; j < k1; j++) {
        _ += a[j] * b[k - j];
    }
    return _;
}

real t_mul (series u, series v, int k) {
    return fa(u, v, 0, k + 1, k);
}

real t_div (series q, series u, series v, int k) {
    CHECK(v[0] != 0.0L); CHECK(q != u && q != v);
    return q[k] = (!k ? (u ? u[0] : 1.0L) : (u ? u[k] : 0.0L) - fa(q, v, 0, k, k)) / v[0];
}

static int half (int k) {
    return 1 + (k - (k % 2 ? 1 : 2)) / 2;
}

static real rem (series a, int k) {
    return k % 2 ? 0.0L : a[k / 2] * a[k / 2];
}

real t_sqr (series u, int k) {
    return 2.0L * fa(u, u, 0, half(k), k) + rem(u, k);
}

real t_sqrt (series r, series u, int k) {
    CHECK(u[0] > 0.0L); CHECK(r != u);
    return r[k] = !k ? sqrtl(u[0]) : 0.5L * (u[k] - 2.0L * fa(r, r, 1, half(k), k) - rem(r, k)) / r[0];
}

static real fb (series df_du, series u, int k) {
    real _ = 0.0L;
    for (int j = 0; j < k; j++) {
        _ += df_du[j] * (k - j) * u[k - j];
    }
    return _ / k;
}

real t_exp (series e, series u, int k) {
    CHECK(e != u);
    return e[k] = !k ? expl(u[0]) : fb(e, u, k);
}

pair t_sin_cos (series s, series c, series u, int k, bool trig) {
    CHECK(s != c && s != u && c != u);
    return !k ? (pair){
        .a = s[0] = trig ? sinl(u[0]) : sinhl(u[0]),
        .b = c[0] = trig ? cosl(u[0]) : coshl(u[0])
    } : (pair){
        .a = s[k] = fb(c, u, k),
        .b = c[k] = fb(s, u, k) * (trig ? -1.0L : 1.0L)
    };
}

pair t_tan_sec2 (series t, series s2, series u, int k, bool trig) {
    CHECK(t != s2 && t != u && s2 != u);
    return !k ? (pair){
        .a =  t[0] = trig ? tanl(u[0]) : tanhl(u[0]),
        .b = s2[0] = trig ? 1.0L + t[0] * t[0] : 1.0L - t[0] * t[0]
    } : (pair){
        .a =  t[k] = fb(s2, u, k),
        .b = s2[k] = fb(t, t, k) * (trig ? 2.0L : -2.0L)
    };
}

real t_pwr (series p, series u, real a, int k) {
    CHECK(u[0] > 0.0L); CHECK(p != u);
    if (!k) return p[0] = powl(u[0], a);
    real _ = 0.0L;
    for (int j = 0; j < k; j++) {
        _ += (a * (k - j) - j) * p[j] * u[k - j];
    }
    return p[k] = _ / (k * u[0]);
}

static real fc (series f, series du_df, series u, int k, bool neg) {
    real _ = 0.0L;
    for (int j = 1; j < k; j++) {
        _ += j * f[j] * du_df[k - j];
    }
    return (u[k] + (neg ? _ : -_) / k) / du_df[0];
}

real t_ln (series l, series u, int k) {
    CHECK(u[0] > 0.0L); CHECK(l != u);
    return l[k] = !k ? logl(u[0]) : fc(l, u, u, k, false);
}

pair t_asin (series as, series g, series u, int k, bool trig) {
    CHECK(trig ? u[0] >= -1.0L && u[0] <= 1.0L : 1); CHECK(as != g && as != u && g != u);
    return !k ? (pair){
        .a = as[0] = trig ? asinl(u[0]) : asinhl(u[0]),
        .b =  g[0] = trig ? sqrtl(1.0L - u[0] * u[0]) : sqrtl(1.0L + u[0] * u[0])
    } : (pair){
        .a = as[k] = fc(as, g, u, k, false),
        .b =  g[k] = fb(u, as, k) * (trig ? -1.0L : 1.0L)
    };
}

pair t_acos (series ac, series g, series u, int k, bool trig) {
    CHECK(trig ? u[0] >= -1.0L && u[0] <= 1.0L : u[0] >= 1.0L); CHECK(ac != g && ac != u && g != u);
    return !k ? (pair){
        .a = ac[0] = trig ? acosl(u[0]) : acoshl(u[0]),
        .b =  g[0] = trig ? - sqrtl(1.0L - u[0] * u[0]) : sqrtl(u[0] * u[0] - 1.0L)
    } : (pair){
        .a = ac[k] = fc(ac, g, u, k, trig),
        .b =  g[k] = fb(u, ac, k)
    };
}

pair t_atan (series at, series g, series u, int k, bool trig) {
    CHECK(trig ? 1 : u[0] >= -1.0L && u[0] <= 1.0L); CHECK(at != g && at != u && g != u);
    return !k ? (pair){
        .a = at[0] = trig ? atanl(u[0]) : atanhl(u[0]),
        .b =  g[0] = trig ? 1.0L + u[0] * u[0] : 1.0L - u[0] * u[0]
    } : (pair){
        .a = at[k] = fc(at, g, u, k, false),
        .b =  g[k] = fb(u, u, k) * (trig ? 2.0L : -2.0L)
    };
}

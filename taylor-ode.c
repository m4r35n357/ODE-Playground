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

static void _out (int dp, real x, real y, real z, real t, char *x_tag, char *y_tag, char *z_tag, clock_t since) {
    real _ = (real)(clock() - since) / CLOCKS_PER_SEC;
    if (dp) printf("%+.*Le %+.*Le %+.*Le %.6Le %s %s %s %.3Lf\n", dp, x, dp, y, dp, z, t, x_tag, y_tag, z_tag, _);
    else printf("%+La %+La %+La %.6Le %s %s %s %.3Lf\n", x, y, z, t, x_tag, y_tag, z_tag, _);
}

controls *tsm_get_c (int argc, char **argv) {
    PRINT_ARGS(argc, argv);
    controls *_ = malloc(sizeof (controls)); CHECK(_);
    _->looping = false;
    _->order = (int)strtol(argv[2], NULL, BASE); CHECK(_->order >= 2 && _->order <= 64);
    _->step_size = strtold(argv[3], NULL);       CHECK(_->step_size > 0.0L);
    _->steps = (int)strtol(argv[4], NULL, BASE); CHECK(_->steps >= 0 && _->steps <= 1000000);
    return _;
}

series3 *tsm_init_xyz (char **argv, int order) {
    series3 *_ = malloc(sizeof (series3)); CHECK(_);
    _->x = t_const(order + 1, strtold(argv[5], NULL));
    _->y = t_const(order + 1, strtold(argv[6], NULL));
    _->z = t_const(order + 1, strtold(argv[7], NULL));
    return _;
}

void tsm_get_p (char **argv, int argc, ...) {
    va_list _;
    va_start(_, argc);
    for (int i = 8; i < argc; i++) *va_arg(_, real *) = strtold(argv[i], NULL);
    va_end(_);
}

series t_jet (int n) {
    CHECK(n > 0);
    series _ = malloc((size_t)n * sizeof (real)); CHECK(_);
    return _;
}

series t_const (int n, real a) {
    series _ = t_jet(n);
    for (int k = 0; k < n; k++) _[k] = !k ? a : 0.0L;
    return _;
}

real t_horner (series s, int n, real h) {
    real _ = 0.0L;
    for (int i = n; i >= 0; i--) _ = _ * h + s[i];
    CHECK(isfinite(_));
    return _;
}

static void _diff (series3 *j, parameters *p, int n) {
    for (int k = 0; k < n; k++) {
        triplet _ = ode(j->x, j->y, j->z, p, k);
        j->x[k + 1] = _.x / (k + 1);
        j->y[k + 1] = _.y / (k + 1);
        j->z[k + 1] = _.z / (k + 1);
    }
}

static void _next (series3 *j, int n, real h) {
    j->x[0] = t_horner(j->x, n, h);
    j->y[0] = t_horner(j->y, n, h);
    j->z[0] = t_horner(j->z, n, h);
}

static char *_tag (series jet, real slope, char *min, char *max) {
    return jet[1] * slope < 0.0L ? (jet[2] > 0.0L ? min : max) : "_";
}

void tsm_stdout (int dp, controls *c, series3 *jets, parameters *p, clock_t t0) {
    real slope_x = 0.0L, slope_y = 0.0L, slope_z = 0.0L;
    for (int step = 0; step < c->steps; step++) {
        _diff(jets, p, c->order);
        _out(dp, jets->x[0], jets->y[0], jets->z[0], c->step_size * step,
             _tag(jets->x, slope_x, "x", "X"), _tag(jets->y, slope_y, "y", "Y"), _tag(jets->z, slope_z, "z", "Z"), t0);
        slope_x = jets->x[1]; slope_y = jets->y[1]; slope_z = jets->z[1];
        _next(jets, c->order, c->step_size);
    }
    _out(dp, jets->x[0], jets->y[0], jets->z[0], c->step_size * c->steps, "_", "_", "_", t0);
}

bool tsm_gen (controls *c, series3 *jets, parameters *p) {
    if (c->looping) goto resume; else c->looping = true;
    for (c->step = 0; c->step < c->steps; c->step++) {
        _diff(jets, p, c->order);
        _next(jets, c->order, c->step_size);
        return true;
        resume: ;
    }
    return c->looping = false;
}

static int _half_ (int k) {
    return 1 + (k - (k % 2 ? 1 : 2)) / 2;
}

static real _rem_ (series a, int k) {
    return k % 2 ? 0.0L : a[k / 2] * a[k / 2];
}

static real _prod_ (series a, series b, int k, int k0, int k1) {
    real _ = 0.0L;
    for (int j = k0; j < k1; j++) _ += a[j] * b[k - j];
    return _;
}

static real _exp_ (series df_du, series u, int k) {
    real _ = 0.0L;
    for (int j = 0; j < k; j++) _ += df_du[j] * (k - j) * u[k - j];
    return _ / k;
}

static real _log_ (series f, series du_df, series u, int k, bool neg) {
    real _ = 0.0L;
    for (int j = 1; j < k; j++) _ += du_df[j] * (k - j) * f[k - j];
    return (u[k] + (neg ? _ : -_) / k) / du_df[0];
}

real t_abs (series u, int k) {
    CHECK(u[0] != 0.0L);
    return u[0] < 0.0L ? -u[k] : u[k];
}

real t_mul (series u, series v, int k) {
    return _prod_(u, v, k, 0, k + 1);
}

real t_div (series q, series u, series v, int k) {
    CHECK(v[0] != 0.0L); CHECK(q != u && q != v);
    return q[k] = (!k ? (u ? u[k] : 1.0L) : (u ? u[k] : 0.0L) - _prod_(q, v, k, 0, k)) / v[0];
}

real t_sqr (series u, int k) {
    return 2.0L * _prod_(u, u, k, 0, _half_(k)) + _rem_(u, k);
}

real t_sqrt (series r, series u, int k) {
    CHECK(u[0] > 0.0L); CHECK(r != u);
    return r[k] = !k ? sqrtl(u[k]) : 0.5L * (u[k] - 2.0L * _prod_(r, r, k, 1, _half_(k)) - _rem_(r, k)) / r[0];
}

real t_exp (series e, series u, int k) {
    CHECK(e != u);
    return e[k] = !k ? expl(u[k]) : _exp_(e, u, k);
}

pair t_sin_cos (series s, series c, series u, int k, bool trig) {
    CHECK(s != c && s != u && c != u);
    return !k ? (pair){
        .a = s[k] = trig ? sinl(u[k]) : sinhl(u[k]),
        .b = c[k] = trig ? cosl(u[k]) : coshl(u[k])
    } : (pair){
        .a = s[k] = _exp_(c, u, k),
        .b = c[k] = _exp_(s, u, k) * (trig ? -1.0L : 1.0L)
    };
}

pair t_tan_sec2 (series t, series s2, series u, int k, bool trig) {
    CHECK(t != s2 && t != u && s2 != u);
    return !k ? (pair){
        .a =  t[k] = trig ? tanl(u[k]) : tanhl(u[k]),
        .b = s2[k] = trig ? 1.0L + t[k] * t[k] : 1.0L - t[k] * t[k]
    } : (pair){
        .a =  t[k] = _exp_(s2, u, k),
        .b = s2[k] = _exp_(t, t, k) * (trig ? 2.0L : -2.0L)
    };
}

real t_pwr (series p, series u, real a, int k) {
    CHECK(u[0] > 0.0L); CHECK(p != u);
    if (!k) return p[k] = powl(u[k], a);
    real _ = 0.0L;
    for (int j = 0; j < k; j++) _ += (a * (k - j) - j) * p[j] * u[k - j];
    return p[k] = _ / (k * u[0]);
}

real t_ln (series l, series u, int k) {
    CHECK(u[0] > 0.0L); CHECK(l != u);
    return l[k] = !k ? logl(u[k]) : _log_(l, u, u, k, false);
}

pair t_asin (series as, series g, series u, int k, bool trig) {
    CHECK(trig ? u[0] >= -1.0L && u[0] <= 1.0L : 1); CHECK(as != g && as != u && g != u);
    return !k ? (pair){
        .a = as[k] = trig ? asinl(u[k]) : asinhl(u[k]),
        .b =  g[k] = trig ? sqrtl(1.0L - u[k] * u[k]) : sqrtl(1.0L + u[k] * u[k])
    } : (pair){
        .a = as[k] = _log_(as, g, u, k, false),
        .b =  g[k] = _exp_(u, as, k) * (trig ? -1.0L : 1.0L)
    };
}

pair t_acos (series ac, series g, series u, int k, bool trig) {
    CHECK(trig ? u[0] >= -1.0L && u[0] <= 1.0L : u[0] >= 1.0L); CHECK(ac != g && ac != u && g != u);
    return !k ? (pair){
        .a = ac[k] = trig ? acosl(u[k]) : acoshl(u[k]),
        .b =  g[k] = trig ? - sqrtl(1.0L - u[k] * u[k]) : sqrtl(u[k] * u[k] - 1.0L)
    } : (pair){
        .a = ac[k] = _log_(ac, g, u, k, trig),
        .b =  g[k] = _exp_(u, ac, k)
    };
}

pair t_atan (series at, series g, series u, int k, bool trig) {
    CHECK(trig ? 1 : u[0] >= -1.0L && u[0] <= 1.0L); CHECK(at != g && at != u && g != u);
    return !k ? (pair){
        .a = at[k] = trig ? atanl(u[k]) : atanhl(u[k]),
        .b =  g[k] = trig ? 1.0L + u[k] * u[k] : 1.0L - u[k] * u[k]
    } : (pair){
        .a = at[k] = _log_(at, g, u, k, false),
        .b =  g[k] = _exp_(u, u, k) * (trig ? 2.0L : -2.0L)
    };
}

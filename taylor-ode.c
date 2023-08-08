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

static void _out_ (int dp, real x, real y, real z, real t, char *x_tag, char *y_tag, char *z_tag, clock_t since) {
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

static void _diff_ (series3 *j, parameters *p, int n) {
    for (int k = 0; k < n; k++) {
        triplet _ = ode(j->x, j->y, j->z, p, k);
        j->x[k + 1] = _.x / (k + 1);
        j->y[k + 1] = _.y / (k + 1);
        j->z[k + 1] = _.z / (k + 1);
    }
}

static void _next_ (series3 *j, int n, real h) {
    j->x[0] = t_horner(j->x, n, h);
    j->y[0] = t_horner(j->y, n, h);
    j->z[0] = t_horner(j->z, n, h);
}

static char *_tag_ (series jet, real slope, char *min, char *max) {
    return jet[1] * slope < 0.0L ? (jet[2] > 0.0L ? min : max) : "_";
}

void tsm_stdout (int dp, controls *c, series3 *jets, parameters *p, clock_t t0) {
    real slope_x = 0.0L, slope_y = 0.0L, slope_z = 0.0L;
    for (int step = 0; step < c->steps; step++) {
        _diff_(jets, p, c->order);
        _out_(dp, jets->x[0], jets->y[0], jets->z[0], c->step_size * step,
             _tag_(jets->x, slope_x, "x", "X"), _tag_(jets->y, slope_y, "y", "Y"), _tag_(jets->z, slope_z, "z", "Z"), t0);
        slope_x = jets->x[1]; slope_y = jets->y[1]; slope_z = jets->z[1];
        _next_(jets, c->order, c->step_size);
    }
    _out_(dp, jets->x[0], jets->y[0], jets->z[0], c->step_size * c->steps, "_", "_", "_", t0);
}

bool tsm_gen (controls *c, series3 *jets, parameters *p) {
    if (c->looping) goto resume; else c->looping = true;
    for (c->step = 0; c->step < c->steps; c->step++) {
        _diff_(jets, p, c->order);
        _next_(jets, c->order, c->step_size);
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

static real _fwd_ (series g, series u, int k) {
    real _ = 0.0L;
    for (int j = 0; j < k; j++) _ += g[j] * (k - j) * u[k - j];
    return _ / k;
}

static real _rev_ (series u, series g, real f_k, int k, bool neg) {
    real _ = 0.0L;
    for (int j = 1; j < k; j++) _ += g[j] * (k - j) * u[k - j];
    return (f_k + (neg ? _ : -_) / k) / g[0];
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
    return e[k] = !k ? expl(u[k]) : _fwd_(e, u, k);
}

pair t_sin_cos (series s, series c, series u, int k, bool trig) {
    CHECK(s != c && s != u && c != u);
    return !k ? (pair){
        .a = s[k] = trig ? sinl(u[k]) : sinhl(u[k]),
        .b = c[k] = trig ? cosl(u[k]) : coshl(u[k])
    } : (pair){
        .a = s[k] = _fwd_(c, u, k),
        .b = c[k] = _fwd_(s, u, k) * (trig ? -1.0L : 1.0L)
    };
}

pair t_tan_sec2 (series t, series s, series u, int k, bool trig) {
    CHECK(t != s && t != u && s != u);
    return !k ? (pair){
        .a = t[k] = trig ? tanl(u[k]) : tanhl(u[k]),
        .b = s[k] = trig ? 1.0L + t[k] * t[k] : 1.0L - t[k] * t[k]
    } : (pair){
        .a = t[k] = _fwd_(s, u, k),
        .b = s[k] = _fwd_(t, t, k) * (trig ? 2.0L : -2.0L)
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
    return l[k] = !k ? logl(u[k]) : _rev_(l, u, u[k], k, false);
}

pair t_asin (series u, series g, series s, int k, bool trig) {
    CHECK(trig ? s[0] >= -1.0L && s[0] <= 1.0L : 1); CHECK(u != g && u != s && g != s);
    return !k ? (pair){
        .a = u[k] = trig ? asinl(s[k]) : asinhl(s[k]),
        .b = g[k] = trig ?  cosl(u[k]) :  coshl(u[k])
    } : (pair){
        .a = u[k] = _rev_(u, g, s[k], k, false),
        .b = g[k] = _fwd_(s, u, k) * (trig ? -1.0L : 1.0L)
    };
}

pair t_acos (series u, series g, series c, int k, bool trig) {
    CHECK(trig ? c[0] >= -1.0L && c[0] <= 1.0L : c[0] >= 1.0L); CHECK(u != g && u != c && g != c);
    return !k ? (pair){
        .a = u[k] = trig ? acosl(c[k]) : acoshl(c[k]),
        .b = g[k] = trig ? -sinl(u[k]) :  sinhl(u[k])
    } : (pair){
        .a = u[k] = _rev_(u, g, c[k], k, trig),
        .b = g[k] = _fwd_(c, u, k)
    };
}

pair t_atan (series u, series g, series t, int k, bool trig) {
    CHECK(trig ? 1 : t[0] >= -1.0L && t[0] <= 1.0L); CHECK(u != g && u != t && g != t);
    return !k ? (pair){
        .a = u[k] = trig ? atanl(t[k]) : atanhl(t[k]),
        .b = g[k] = trig ? 1.0L + t[k] * t[k] : 1.0L - t[k] * t[k]
    } : (pair){
        .a = u[k] = _rev_(u, g, t[k], k, false),
        .b = g[k] = _fwd_(t, t, k) * (trig ? 2.0L : -2.0L)
    };
}

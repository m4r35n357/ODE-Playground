/*
 * Arbitrary-Order Taylor Series Integrator, together with a well-tested set of recurrence relations
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include <math.h>
#include "taylor-ode.h"

controls *tsm_get_c (int argc, char **argv) {
    PRINT_ARGS(argc, argv);
    controls *_ = malloc(sizeof (controls)); CHECK(_);
    _->dp = (int)strtol(argv[1], NULL, BASE);    CHECK(_->dp >= 0);
    _->order = (int)strtol(argv[2], NULL, BASE); CHECK(_->order >= 2 && _->order <= 64);
    _->h = strtold(argv[3], NULL);               CHECK(_->h > 0.0L);
    _->steps = (int)strtol(argv[4], NULL, BASE); CHECK(_->steps >= 0 && _->steps <= 1000000);
    _->looping = false;
    return _;
}

xyz *tsm_init (char **argv, int o) {
    xyz *_ = malloc(sizeof (xyz)); CHECK(_);
    _->x = tsm_jet(o + 1); _->x[0] = strtold(argv[5], NULL);
    _->y = tsm_jet(o + 1); _->y[0] = strtold(argv[6], NULL);
    _->z = tsm_jet(o + 1); _->z[0] = strtold(argv[7], NULL);
    return _;
}

void tsm_get_p (char **argv, int argc, ...) {
    va_list _;
    va_start(_, argc);
    for (int i = 8; i < argc; i++) *va_arg(_, real *) = strtold(argv[i], NULL);
    va_end(_);
}

series tsm_jet (int n) {
    CHECK(n > 0);
    series _ = malloc((size_t)n * sizeof (real)); CHECK(_);
    for (int i = 0; i < n; i++) _[i] = 0.0L;
    return _;
}

real horner (const series u, int o, real h) {
    real _ = 0.0L;
    for (int i = o; i >= 0; i--) _ = _ * h + u[i];
    CHECK(isfinite(_));
    return _;
}

static void _diff_ (xyz *_, const model *p, int o) {
    for (int k = 0; k < o; k++) {
        triplet v = ode(_->x, _->y, _->z, p, k);
        _->x[k + 1] = v.x / (k + 1);
        _->y[k + 1] = v.y / (k + 1);
        _->z[k + 1] = v.z / (k + 1);
    }
}

static void _next_ (xyz *_, int o, real h) {
    _->x[0] = horner(_->x, o, h);
    _->y[0] = horner(_->y, o, h);
    _->z[0] = horner(_->z, o, h);
}

bool tsm_gen (controls *c, xyz *_, const model *p) {
    if (c->looping) goto resume; else c->looping = true;
    for (c->step = 0; c->step < c->steps; c->step++) {
        _diff_(_, p, c->order);
        _next_(_, c->order, c->h);
        return true;
        resume: ;
    }
    return c->looping = false;
}

static void _out_ (int dp, real x, real y, real z, real t, char x_tag, char y_tag, char z_tag, clock_t since) {
    real cpu = (real)(clock() - since) / CLOCKS_PER_SEC;
    if (dp) {
        printf("%+.*Le %+.*Le %+.*Le %.6Le %c %c %c %.3Lf\n", dp, x, dp, y, dp, z, t, x_tag, y_tag, z_tag, cpu);
    } else {
        printf("%+La %+La %+La %.6Le %c %c %c %.3Lf\n", x, y, z, t, x_tag, y_tag, z_tag, cpu);
    }
}

static char _tp_ (series u, real *v_old, char min) {
    char tag = *v_old * u[1] >= 0.0L ? '_' : (u[2] > 0.0L ? min : (char)toupper(min));
    *v_old = u[1];
    return tag;
}

void tsm (controls *c, xyz *_, const model *p, clock_t t0) {
    real vX = 0.0L, vY = 0.0L, vZ = 0.0L;
    for (int step = 0; step < c->steps; step++) {
        _diff_(_, p, c->order);
        _out_(c->dp, _->x[0], _->y[0], _->z[0], c->h * step, _tp_(_->x, &vX, 'x'), _tp_(_->y, &vY, 'y'), _tp_(_->z, &vZ, 'z'), t0);
        _next_(_, c->order, c->h);
    }
    _out_(c->dp, _->x[0], _->y[0], _->z[0], c->h * c->steps, '_', '_', '_', t0);
}

real t_const (const real value, int k) {
    return k ? 0.0L : value;
}

real t_abs (const series u, int k) {
    if (!k) CHECK(u[0] != 0.0L);
    return u[0] < 0.0L ? -u[k] : u[k];
}

static real _cauchy_ (const series b, const series a, int k, int k0, int k1) {
    real _ = 0.0L;
    for (int j = k0; j <= k1; j++) _ += b[j] * a[k - j];
    return _;
}

real t_mul (const series u, const series v, int k) {
    return _cauchy_(u, v, k, 0, k);
}

real t_div (series q, const series u, const series v, int k) {
    if (k) return q[k] = ((u ? u[k] : 0.0L) - _cauchy_(q, v, k, 0, k - 1)) / v[0];
    CHECK(q != u && q != v);
    CHECK(v[0] != 0.0L);
    return q[k] = (u ? u[k] : 1.0L) / v[0];
}

static real _half_ (const series a, int k, int k0, bool even) {
    return 2.0L * _cauchy_(a, a, k, k0, (k - (even ? 1 : 2)) / 2) + (even ? 0.0L : SQR(a[k / 2]));
}

real t_sqr (const series u, int k) {
    return _half_(u, k, 0, k % 2);
}

real t_sqrt (series r, const series u, int k) {
    if (k) return r[k] = 0.5L * (u[k] - _half_(r, k, 1, k % 2)) / r[0];
    CHECK(r != u);
    CHECK(u[0] > 0.0L);
    return r[k] = sqrtl(u[k]);
}

real t_pwr (series p, const series u, real a, int k) {
    if (k) {
        real _ = 0.0L;
        for (int j = 0; j < k; j++) _ += (a * (k - j) - j) * p[j] * u[k - j];
        return p[k] = _ / (k * u[0]);
    }
    CHECK(p != u);
    CHECK(u[0] > 0.0L);
    return p[k] = powl(u[k], a);
}

static real _chain_ (const series dfdu, const series u, int k, const series fk, int scale) {
    real _ = 0.0L;
    for (int j = fk ? 1 : 0; j < k; j++) _ += dfdu[j] * (k - j) * u[k - j];
    return fk ? (*fk - scale * _ / k) / dfdu[0] : scale * _ / k;  // forward if fk NULL, reverse if non-NULL
}

real t_exp (series e, const series u, int k) {
    if (k) return e[k] = _chain_(e, u, k, NULL, 1);
    CHECK(e != u);
    return e[k] = expl(u[k]);
}

real t_ln (series u, const series e, int k) {
    if (k) return u[k] = _chain_(e, u, k, &e[k], 1);
    CHECK(u != e);
    CHECK(e[0] > 0.0L);
    return u[k] = logl(e[k]);
}

pair t_sin_cos (series s, series c, const series u, int k, bool trig) {
    if (k) return (pair){ s[k] = _chain_(c, u, k, NULL, 1), c[k] = _chain_(s, u, k, NULL, trig ? -1.0L : 1.0L) };
    CHECK(s != c && s != u && c != u);
    return (pair){ s[k] = trig ? sinl(u[k]) : sinhl(u[k]), c[k] = trig ? cosl(u[k]) : coshl(u[k]) };
}

pair t_tan_sec2 (series t, series s, const series u, int k, bool trig) {
    if (k) return (pair){ t[k] = _chain_(s, u, k, NULL, 1), s[k] = _chain_(t, t, k, NULL, trig ? 2.0L : -2.0L) };
    CHECK(t != s && t != u && s != u);
    CHECK(trig ? fabsl(u[0]) < 0.5L * acosl(-1.0L) : true);
    return (pair){ t[k] = trig ? tanl(u[k]) : tanhl(u[k]), s[k] = trig ? 1.0L + SQR(t[k]) : 1.0L - SQR(t[k]) };
}

pair t_asin_cos (series u, series c, const series s, int k, bool trig) {
    if (k) return (pair){ u[k] = _chain_(c, u, k, &s[k], 1), c[k] = _chain_(s, u, k, NULL, trig ? -1.0L : 1.0L) };
    CHECK(u != c && u != s && c != s);
    CHECK(trig ? s[0] > -1.0L && s[0] < 1.0L : true);
    return (pair){ u[k] = trig ? asinl(s[k]) : asinhl(s[k]), c[k] = trig ?  cosl(u[k]) :  coshl(u[k]) };
}

pair t_acos_sin (series u, series s, const series c, int k, bool trig) {
    if (k) return (pair){ u[k] = _chain_(s, u, k, &c[k], trig ? -1.0L : 1.0L), s[k] = _chain_(c, u, k, NULL, 1) };
    CHECK(u != s && u != c && s != c);
    CHECK(trig ? c[0] > -1.0L && c[0] < 1.0L : c[0] > 1.0L);
    return (pair){ u[k] = trig ? acosl(c[k]) : acoshl(c[k]), s[k] = trig ? -sinl(u[k]) :  sinhl(u[k]) };
}

pair t_atan_sec2 (series u, series s, const series t, int k, bool trig) {
    if (k) return (pair){ u[k] = _chain_(s, u, k, &t[k], 1), s[k] = _chain_(t, t, k, NULL, trig ? 2.0L : -2.0L) };
    CHECK(u != s && u != t && s != t);
    CHECK(trig ? true : t[0] > -1.0L && t[0] < 1.0L);
    return (pair){ u[k] = trig ? atanl(t[k]) : atanhl(t[k]), s[k] = trig ? 1.0L + SQR(t[k]) : 1.0L - SQR(t[k]) };
}

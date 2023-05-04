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

const int BASE = 10;

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

static void tsm_step (series3 *j, void *p, int n, real h) {
    for (int k = 0; k < n; k++) {
        components v = ode(j->x, j->y, j->z, p, k);
        j->x[k + 1] = v.x / (k + 1);
        j->y[k + 1] = v.y / (k + 1);
        j->z[k + 1] = v.z / (k + 1);
    }
    j->x[0] = t_horner(j->x, n, h);
    j->y[0] = t_horner(j->y, n, h);
    j->z[0] = t_horner(j->z, n, h);
}

void tsm_stdout (int dp, controls *c, series3 *jets, void *p, clock_t t0) {
    CHECK(dp > 0); CHECK(c); CHECK(jets); CHECK(p); CHECK(t0);
    components s = {0.0L, 0.0L, 0.0L};
    for (int step = 0; step < c->steps; step++) {
        t_out(dp, jets->x[0], jets->y[0], jets->z[0], c->step_size * step,
              tag(jets->x, s.x, "x", "X"), tag(jets->y, s.y, "y", "Y"), tag(jets->z, s.z, "z", "Z"), t0);
        s = (components){jets->x[1], jets->y[1], jets->z[1]};
        tsm_step(jets, p, c->order, c->step_size);
    }
    t_out(dp, jets->x[0], jets->y[0], jets->z[0], c->step_size * c->steps, "_", "_", "_", t0);
}

bool tsm_gen (controls *c, series3 *jets, void *p) {
    CHECK(c); CHECK(jets); CHECK(p);
    static bool looping = false;
    if (looping) goto resume; else looping = true;
    for (c->step = 0; c->step < c->steps; c->step++) {
        tsm_step(jets, p, c->order, c->step_size);
        return true;
        resume: ;
    }
    return looping = false;
}

real t_abs (series u, int k) {
    return u[0] < 0.0L ? -u[k] : u[k];
}

real t_mul (series u, series v, int k) {
    real _m = 0.0L;
    for (int j = 0; j <= k; j++) {
        _m += u[j] * v[k - j];
    }
    return _m;
}

real t_sqr (series u, int k) {
    real _s = 0.0L;
    for (int j = 0; j <= (k - (k % 2 ? 1 : 2)) / 2; j++) {
        _s += u[j] * u[k - j];
    }
    return 2.0L * _s + (k % 2 ? 0.0L : u[k / 2] * u[k / 2]);
}

real t_div (series q, series u, series v, int k) {
    CHECK(v[0] != 0.0L); CHECK(q != u && q != v);
    if (!k) return q[0] = (u ? u[0] : 1.0L) / v[0];
    real _d = 0.0L;
    for (int j = 0; j < k; j++) {
        _d += q[j] * v[k - j];
    }
    return q[k] = ((u ? u[k] : 0.0L) - _d) / v[0];
}

real t_sqrt (series r, series u, int k) {
    CHECK(u[0] > 0.0L); CHECK(r != u);
    if (!k) return r[0] = sqrtl(u[0]);
    real _r = 0.0L;
    for (int j = 1; j <= (k - (k % 2 ? 1 : 2)) / 2; j++) {
        _r += r[j] * r[k - j];
    }
    return r[k] = 0.5L * (u[k] - 2.0L * _r - (k % 2 ? 0.0L : r[k / 2] * r[k / 2])) / r[0];
}

real t_exp (series e, series u, int k) {
    CHECK(e != u);
    if (!k) return e[0] = expl(u[0]);
    real _e = 0.0L;
    for (int j = 0; j < k; j++) {
        _e += e[j] * (k - j) * u[k - j];
    }
    return e[k] = _e / k;
}

pair t_sin_cos (series s, series c, series u, int k, bool trig) {
    CHECK(s != c && s != u && c != u);
    if (!k) return (pair){s[0] = trig ? sinl(u[0]) : sinhl(u[0]), c[0] = trig ? cosl(u[0]) : coshl(u[0])};
    real _s = 0.0L, _c = 0.0L;
    for (int j = 0; j < k; j++) {
        real _ = (k - j) * u[k - j];
        _c += c[j] * _;
        _s += s[j] * _;
    }
    return (pair){s[k] = _c / k, c[k] = (trig ? -_s : _s) / k};
}

pair t_tan_sec2 (series t, series s, series u, int k, bool trig) {
    CHECK(t != s && t != u && s != u);
    if (!k) {
        t[0] = trig ? tanl(u[0]) : tanhl(u[0]);
        return (pair){t[0], s[0] = trig ? 1.0L + t[0] * t[0] : 1.0L - t[0] * t[0]};
    }
    real _t = 0.0L, _s = 0.0L;
    for (int j = 0; j < k; j++) {
        _s += s[j] * (k - j) * u[k - j];
    }
    t[k] = _s / k;
    for (int j = 0; j < k; j++) {
        _t += t[j] * (k - j) * t[k - j];
    }
    return (pair){t[k], s[k] = 2.0L * (trig ? _t : -_t) / k};
}

real t_pwr (series p, series u, real a, int k) {
    CHECK(u[0] > 0.0L); CHECK(p != u);
    if (!k) return p[0] = powl(u[0], a);
    real _p = 0.0L;
    for (int j = 0; j < k; j++) {
        _p += (a * (k - j) - j) * p[j] * u[k - j];
    }
    return p[k] = _p / (k * u[0]);
}

real t_ln (series l, series u, int k) {
    CHECK(u[0] > 0.0L); CHECK(l != u);
    if (!k) return l[0] = logl(u[0]);
    real _l = 0.0L;
    for (int j = 1; j < k; j++) {
        _l += j * l[j] * u[k - j];
    }
    return l[k] = (u[k] - _l / k) / u[0];
}

pair t_asin (series a, series g, series u, int k, bool trig) {
    CHECK(trig ? u[0] >= -1.0L && u[0] <= 1.0L : 1); CHECK(a != g && a != u && g != u);
    if (!k) return (pair){
        a[0] = trig ? asinl(u[0]) : asinhl(u[0]),
        g[0] = trig ? sqrtl(1.0L - u[0] * u[0]) : sqrtl(u[0] * u[0] + 1.0L)
    };
    real _a = 0.0L, _g = 0.0L;
    for (int j = 1; j < k; j++) {
        _a += j * a[j] * g[k - j];
    }
    a[k] = (u[k] - _a / k) / g[0];
    for (int j = 0; j < k; j++) {
        _g += u[j] * (k - j) * a[k - j];
    }
    return (pair){a[k], g[k] = (trig ? -_g : _g) / k};
}

pair t_acos (series a, series g, series u, int k, bool trig) {
    CHECK(trig ? u[0] >= -1.0L && u[0] <= 1.0L : u[0] >= 1.0L); CHECK(a != g && a != u && g != u);
    if (!k) return (pair){
        a[0] = trig ? acosl(u[0]) : acoshl(u[0]),
        g[0] = trig ? - sqrtl(1.0L - u[0] * u[0]) : sqrtl(u[0] * u[0] - 1.0L)
    };
    real _a = 0.0L, _g = 0.0L;
    for (int j = 1; j < k; j++) {
        _a += j * a[j] * g[k - j];
    }
    a[k] = (u[k] - (trig ? -_a : _a) / k) / g[0];
    for (int j = 0; j < k; j++) {
        _g += u[j] * (k - j) * a[k - j];
    }
    return (pair){a[k], g[k] = _g / k};
}

pair t_atan (series a, series g, series u, int k, bool trig) {
    CHECK(trig ? 1 : u[0] >= -1.0L && u[0] <= 1.0L); CHECK(a != g && a != u && g != u);
    if (!k) return (pair){
        a[0] = trig ? atanl(u[0]) : atanhl(u[0]),
        g[0] = trig ? 1.0L + u[0] * u[0] : 1.0L - u[0] * u[0]
    };
    real _a = 0.0L, _g = 0.0L;
    for (int j = 1; j < k; j++) {
        _a += j * a[j] * g[k - j];
    }
    for (int j = 0; j < k; j++) {
        _g += u[j] * (k - j) * u[k - j];
    }
    return (pair){a[k] = (u[k] - _a / k) / g[0], g[k] = 2.0L * (trig ? _g : -_g) / k};
}

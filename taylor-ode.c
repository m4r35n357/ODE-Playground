/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "taylor-ode.h"

series t_jet (int n) {
    series s = calloc((size_t)n, sizeof (real));
    if (!s) {
        fprintf(stderr, "Allocation failure!\n");
        exit(1);
    }
    return s;
}

real t_horner (series s, int n, real h) {
    real sum = 0.0L;
    for (int i = n; i >= 0; i--) {
        sum = sum * h + s[i];
    }
    if (! isfinite(sum)) {
        fprintf(stderr, "Value error!\n");
        exit(2);
    }
    return sum;
}

static char *tag (series jet, real slope, char *min, char *max) {
    return jet[1] * slope < 0.0L ? (jet[2] > 0.0L ? min : max) : "_";
}

void tsm (int dp, int n, real h, int steps, real x0, real y0, real z0, void *p, clock_t t0) {
    series x = t_jet(n + 1); x[0] = x0;
    series y = t_jet(n + 1); y[0] = y0;
    series z = t_jet(n + 1); z[0] = z0;
    components s = (components) {0.0L, 0.0L, 0.0L};
    for (int step = 0; step < steps; step++) {
        t_out(dp, x[0], y[0], z[0], h * step, tag(x, s.x, "x", "X"), tag(y, s.y, "y", "Y"), tag(z, s.z, "z", "Z"), t0);
        s = (components) {x[1], y[1], z[1]};
        for (int k = 0; k < n; k++) {
            components vk = ode(x, y, z, p, k);
            x[k + 1] = vk.x / (k + 1);
            y[k + 1] = vk.y / (k + 1);
            z[k + 1] = vk.z / (k + 1);
        }
        x[0] = t_horner(x, n, h);
        y[0] = t_horner(y, n, h);
        z[0] = t_horner(z, n, h);
    }
    t_out(dp, x[0], y[0], z[0], h * steps, "_", "_", "_", t0);
}

real t_const (real a, int k) {
    return k == 0 ? a : 0.0L;
}

real t_abs (series u, int k) {
    return u[0] < 0.0L ? - u[k] : u[k];
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
    assert(v[0] != 0.0L);
    assert(q != u && q != v);
    if (k == 0) {
        return q[0] = (u ? u[0] : 1.0L) / v[0];
    } else {
        real _d = 0.0L;
        for (int j = 0; j < k; j++) {
            _d += q[j] * v[k - j];
        }
        return q[k] = ((u ? u[k] : 0.0L) - _d) / v[0];
    }
}

real t_sqrt (series r, series u, int k) {
    assert(u[0] > 0.0L);
    assert(r != u);
    if (k == 0) {
        return r[0] = sqrtl(u[0]);
    } else {
        real _r = 0.0L;
        for (int j = 1; j <= (k - (k % 2 ? 1 : 2)) / 2; j++) {
            _r += r[j] * r[k - j];
        }
        return r[k] = 0.5L * (u[k] - 2.0L * _r - (k % 2 ? 0.0L : r[k / 2] * r[k / 2])) / r[0];
    }
}

real t_exp (series e, series u, int k) {
    assert(e != u);
    if (k == 0) {
        return e[0] = expl(u[0]);
    } else {
        real _e = 0.0L;
        for (int j = 0; j < k; j++) {
            _e += e[j] * (k - j) * u[k - j];
        }
        return e[k] = _e / k;
    }
}

pair t_sin_cos (series s, series c, series u, int k, geometry G) {
    assert(s != c && s != u && c != u);
    if (k == 0) {
        return (pair) {s[0] = G == TRIG ? sinl(u[0]) : sinhl(u[0]), c[0] = G == TRIG ? cosl(u[0]) : coshl(u[0])};
    } else {
        real _s = 0.0L, _c = 0.0L;
        for (int j = 0; j < k; j++) {
            real du_dt = (k - j) * u[k - j];
            _c += c[j] * du_dt;
            _s += s[j] * du_dt;
        }
        return (pair) {s[k] = _c / k, c[k] = (G == TRIG ? - _s : _s) / k};
    }
}

pair t_tan_sec2 (series t, series s, series u, int k, geometry G) {
    assert(t != s && t != u && s != u);
    if (k == 0) {
        t[0] = G == TRIG ? tanl(u[0]) : tanhl(u[0]);
        return (pair) {t[0], s[0] = G == TRIG ? 1.0L + t[0] * t[0] : 1.0L - t[0] * t[0]};
    } else {
        real _t = 0.0L, _s = 0.0L;
        for (int j = 0; j < k; j++) {
            _s += s[j] * (k - j) * u[k - j];
        }
        t[k] = _s / k;
        for (int j = 0; j < k; j++) {
            _t += t[j] * (k - j) * t[k - j];
        }
        return (pair) {t[k], s[k] = 2.0L * (G == TRIG ? _t : - _t) / k};
    }
}

real t_pwr (series p, series u, real a, int k) {
    assert(u[0] > 0.0L);
    assert(p != u);
    if (k == 0) {
        return p[0] = powl(u[0], a);
    } else {
        real _p = 0.0L;
        for (int j = 0; j < k; j++) {
            _p += (a * (k - j) - j) * p[j] * u[k - j];
        }
        return p[k] = _p / (k * u[0]);
    }
}

real t_ln (series l, series u, int k) {
    assert(u[0] > 0.0L);
    assert(l != u);
    if (k == 0) {
        return l[0] = logl(u[0]);
    } else {
        real _l = 0.0L;
        for (int j = 1; j < k; j++) {
            _l += j * l[j] * u[k - j];
        }
        return l[k] = (u[k] - _l / k) / u[0];
    }
}

pair t_asin (series a, series g, series u, int k, geometry G) {
    assert(G == TRIG ? u[0] >= -1.0L && u[0] <= 1.0L : 1);
    assert(a != g && a != u && g != u);
    if (k == 0) {
        return (pair) {
            a[0] = G == TRIG ? asinl(u[0]) : asinhl(u[0]),
            g[0] = sqrtl(G == TRIG ? 1.0L - u[0] * u[0] : 1.0L + u[0] * u[0])
        };
    } else {
        real _a = 0.0L, _g = 0.0L;
        for (int j = 1; j < k; j++) {
            _a += j * a[j] * g[k - j];
        }
        a[k] = (u[k] - _a / k) / g[0];
        for (int j = 0; j < k; j++) {
            _g += u[j] * (k - j) * a[k - j];
        }
        return (pair) {a[k], g[k] = (G == TRIG ? - _g : _g) / k};
    }
}

pair t_acos (series a, series g, series u, int k, geometry G) {
    assert(G == TRIG ? u[0] >= -1.0L && u[0] <= 1.0L : u[0] >= 1.0L);
    assert(a != g && a != u && g != u);
    if (k == 0) {
        return (pair) {
            a[0] = G == TRIG ? acosl(u[0]) : acoshl(u[0]),
            g[0] = G == TRIG ? - sqrtl(1.0L - u[0] * u[0]) : sqrtl(u[0] * u[0] - 1.0L)
        };
    } else {
        real _a = 0.0L, _g = 0.0L;
        for (int j = 1; j < k; j++) {
            _a += j * a[j] * g[k - j];
        }
        a[k] = (u[k] - (G == TRIG ? - _a : _a) / k) / g[0];
        for (int j = 0; j < k; j++) {
            _g += u[j] * (k - j) * a[k - j];
        }
        return (pair) {a[k], g[k] = _g / k};
    }
}

pair t_atan (series a, series g, series u, int k, geometry G) {
    assert(G == TRIG ? 1 : u[0] >= -1.0L && u[0] <= 1.0L);
    assert(a != g && a != u && g != u);
    if (k == 0) {
        return (pair) {
            a[0] = G == TRIG ? atanl(u[0]) : atanhl(u[0]),
            g[0] = G == TRIG ? 1.0L + u[0] * u[0] : 1.0L - u[0] * u[0]
        };
    } else {
        real _a = 0.0L, _g = 0.0L;
        for (int j = 0; j < k; j++) {
            _a += j > 0 ? j * a[j] * g[k - j] : 0.0L;
            _g += u[j] * (k - j) * u[k - j];
        }
        return (pair) {a[k] = (u[k] - _a / k) / g[0], g[k] = 2.0L * (G == TRIG ? _g : - _g) / k};
    }
}

/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include "taylor-ode.h"

series t_jet (int n) {
    series s = calloc((size_t)n + 1, sizeof (real));
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
    if (isnan(sum) || isinf(sum)) {
        fprintf(stderr, "Value error!\n");
        exit(2);
    }
    return sum;
}

static char *tag (series jet, real slope, char *min, char *max) {
    return jet[1] * slope < 0.0L ? (jet[2] > 0.0L ? min : max) : "_";
}

void tsm (int dp, int n, real h, int steps, real x0, real y0, real z0, void *p) {
    series x = t_jet(n); x[0] = x0;
    series y = t_jet(n); y[0] = y0;
    series z = t_jet(n); z[0] = z0;
    components s = (components) {0.0L, 0.0L, 0.0L};
    for (int step = 0; step < steps; step++) {
        for (int k = 0; k < n; k++) {
            components vk = ode(x, y, z, p, k);
            x[k + 1] = vk.x / (k + 1);
            y[k + 1] = vk.y / (k + 1);
            z[k + 1] = vk.z / (k + 1);
        }
        t_output(dp, x[0], y[0], z[0], h * step, tag(x, s.x, "x", "X"), tag(y, s.y, "y", "Y"), tag(z, s.z, "z", "Z"));
        s = (components) {x[1], y[1], z[1]};
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

real t_mul (series u, series v, int k) {
    real sum = 0.0L;
    for (int j = 0; j <= k; j++) {
        sum += u[j] * v[k - j];
    }
    return sum;
}

real t_sqr (series u, int k) {
    real sum = 0.0L;
    for (int j = 0; j <= (k - (k % 2 ? 1 : 2)) / 2; j++) {
        sum += u[j] * u[k - j];
    }
    return 2.0L * sum + (k % 2 ? 0.0L : u[k / 2] * u[k / 2]);
}

real t_div (series q, series u, series v, int k) {
    assert(v[0] != 0.0L);
    assert(q != u && q != v);
    if (k == 0) {
        return q[0] = (u ? u[0] : 1.0L) / v[0];
    } else {
        real sum = 0.0L;
        for (int j = 0; j < k; j++) {
            sum += q[j] * v[k - j];
        }
        return q[k] = ((u ? u[k] : 0.0L) - sum) / v[0];
    }
}

real t_sqrt (series r, series u, int k) {
    assert(u[0] > 0.0L);
    assert(r != u);
    if (k == 0) {
        return r[0] = sqrtl(u[0]);
    } else {
        real sum = 0.0L;
        for (int j = 1; j <= (k - (k % 2 ? 1 : 2)) / 2; j++) {
            sum += r[j] * r[k - j];
        }
        return r[k] = 0.5L * (u[k] - 2.0L * sum - (k % 2 ? 0.0L : r[k / 2] * r[k / 2])) / r[0];
    }
}

real t_exp (series e, series u, int k) {
    assert(e != u);
    if (k == 0) {
        return e[0] = expl(u[0]);
    } else {
        real sum = 0.0L;
        for (int j = 0; j < k; j++) {
            sum += e[j] * (k - j) * u[k - j];
        }
        return e[k] = sum / k;
    }
}

pair t_sin_cos (series s, series c, series u, int k, geometry g) {
    assert(s != c && s != u && c != u);
    if (k == 0) {
        return (pair) {s[0] = g == TRIG ? sinl(u[0]) : sinhl(u[0]), c[0] = g == TRIG ? cosl(u[0]) : coshl(u[0])};
    } else {
        real s_sum = 0.0L, c_sum = 0.0L;
        for (int j = 0; j < k; j++) {
            real du_dt = (k - j) * u[k - j];
            c_sum += c[j] * du_dt;
            s_sum += s[j] * du_dt;
        };
        return (pair) {s[k] = c_sum / k, c[k] = (g == TRIG ? - s_sum : s_sum) / k};
    }
}

pair t_tan_sec2 (series t, series s, series u, int k, geometry g) {
    assert(t != s && t != u && s != u);
    if (k == 0) {
        t[0] = g == TRIG ? tanl(u[0]) : tanhl(u[0]);
        return (pair) {t[0], s[0] = g == TRIG ? 1.0L + t[0] * t[0] : 1.0L - t[0] * t[0]};
    } else {
        real t_sum = 0.0L, s_sum = 0.0L;
        for (int j = 0; j < k; j++) {
            s_sum += s[j] * (k - j) * u[k - j];
        };
        t[k] = s_sum / k;
        for (int j = 0; j < k; j++) {
            t_sum += t[j] * (k - j) * t[k - j];
        };
        return (pair) {t[k], s[k] = 2.0L * (g == TRIG ? t_sum : - t_sum) / k};
    }
}

real t_pwr (series p, series u, real a, int k) {
    assert(u[0] > 0.0L);
    assert(p != u);
    if (k == 0) {
        return p[0] = powl(u[0], a);
    } else {
        real sum = 0.0L;
        for (int j = 0; j < k; j++) {
            sum += (a * (k - j) - j) * p[j] * u[k - j];
        }
        return p[k] = sum / (k * u[0]);
    }
}

real t_ln (series l, series u, int k) {
    assert(u[0] > 0.0L);
    assert(l != u);
    if (k == 0) {
        return l[0] = logl(u[0]);
    } else {
        real sum = 0.0L;
        for (int j = 1; j < k; j++) {
            sum += j * l[j] * u[k - j];
        }
        return l[k] = (u[k] - sum / k) / u[0];
    }
}

pair t_asin (series as, series uf, series u, int k, geometry g) {
    assert(g == TRIG ? u[0] >= -1.0L && u[0] <= 1.0L : 1);
    assert(as != uf && as != u && uf != u);
    if (k == 0) {
        return (pair) {
            as[0] = g == TRIG ? asinl(u[0]) : asinhl(u[0]),
            uf[0] = sqrtl(g == TRIG ? 1.0L - u[0] * u[0] : 1.0L + u[0] * u[0])
        };
    } else {
        real as_sum = 0.0L, uf_sum = 0.0L;
        for (int j = 1; j < k; j++) {
            as_sum += j * as[j] * uf[k - j];
        }
        as[k] = (u[k] - as_sum / k) / uf[0];
        for (int j = 0; j < k; j++) {
            uf_sum += u[j] * (k - j) * as[k - j];
        }
        return (pair) {as[k], uf[k] = (g == TRIG ? - uf_sum : uf_sum) / k};
    }
}

pair t_acos (series ac, series uf, series u, int k, geometry g) {
    assert(g == TRIG ? u[0] >= -1.0L && u[0] <= 1.0L : u[0] >= 1.0L);
    assert(ac != uf && ac != u && uf != u);
    if (k == 0) {
        return (pair) {
            ac[0] = g == TRIG ? acosl(u[0]) : acoshl(u[0]),
            uf[0] = g == TRIG ? - sqrtl(1.0L - u[0] * u[0]) : sqrtl(u[0] * u[0] - 1.0L)
        };
    } else {
        real ac_sum = 0.0L, uf_sum = 0.0L;
        for (int j = 1; j < k; j++) {
            ac_sum += j * ac[j] * uf[k - j];
        }
        ac[k] = (u[k] + (g == TRIG ? ac_sum : - ac_sum) / k) / uf[0];
        for (int j = 0; j < k; j++) {
            uf_sum += u[j] * (k - j) * ac[k - j];
        }
        return (pair) {ac[k], uf[k] = uf_sum / k};
    }
}

pair t_atan (series at, series uf, series u, int k, geometry g) {
    assert(g == TRIG ? 1 : u[0] >= -1.0L && u[0] <= 1.0L);
    assert(at != uf && at != u && uf != u);
    if (k == 0) {
        return (pair) {
            at[0] = g == TRIG ? atanl(u[0]) : atanhl(u[0]),
            uf[0] = g == TRIG ? 1.0L + u[0] * u[0] : 1.0L - u[0] * u[0]
        };
    } else {
        real at_sum = 0.0L, uf_sum = 0.0L;
        for (int j = 1; j < k; j++) {
            at_sum += j * at[j] * uf[k - j];
        }
        for (int j = 0; j < k; j++) {
            uf_sum += u[j] * (k - j) * u[k - j];
        }
        return (pair) {
            at[k] = (u[k] - at_sum / k) / uf[0],
            uf[k] = 2.0L * (g == TRIG ? uf_sum : - uf_sum) / k
        };
    }
}

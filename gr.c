/*
 * GR Geodesic equations, solved using Taylor Series integration
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include "gr.h"

series *t_jet4 (long n, vector4 a) {
    series *jet4 = calloc(4, sizeof (series));
    for (int mu = 0; mu < 4; mu++) {
        jet4[mu] = t_jet_c(n + 1, a[mu]);
    }
    return jet4;
}

void mm_mult (matrix4x4 c, matrix4x4 a, matrix4x4 b) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            real sum = 0.0L;
            for (int k = 0; k < 4; k++) {
                sum += a[i][k] * b[k][j];
            }
            c[i][j] = sum;
        }
    }
}

real r_ra2 (real a, real r) {
    return r * r + a * a;
}

real r_delta (real m, real a, real r) {
    return r * r - 2.0L * m * r + a * a;
}

real r_sigma (real a, real r, real theta) {
    return r * r + a * a * cosl(theta) * cosl(theta);
}

dual d_ra2 (real a, dual r) {
    return d_shift(d_sqr(r), a * a);
}

dual d_delta (real m, real a, dual r) {
    return d_add(d_sqr(r), d_shift(d_scale(r, - 2.0L * m), a * a));
}

dual d_sigma (real a, dual r, dual theta) {
    return d_add(d_sqr(r), d_scale(d_sqr(d_cos(theta)), a * a));
}

void gr_control (char **argv, long *dp, long *n, real *h, long *nsteps,
                real *t, real *r, real *theta, real *phi, real *t_dot, real *r_dot, real *theta_dot, real *phi_dot) {
    *dp = strtol(argv[1], NULL, 10);
    *n = strtol(argv[2], NULL, 10);
    assert(*n >= 1 && *n <= 34);
    *h = strtold(argv[3], NULL);
    assert(*h > 0.0L && *h <= 1.0L);
    *nsteps = strtol(argv[4], NULL, 10);
    assert(*nsteps >= 1 && *nsteps <= 100000);
    *t = strtold(argv[5], NULL);
    *r = strtold(argv[6], NULL);
    *theta = strtold(argv[7], NULL);
    *phi = strtold(argv[8], NULL);
    *t_dot = strtold(argv[9], NULL);
    *r_dot = strtold(argv[10], NULL);
    *theta_dot = strtold(argv[11], NULL);
    *phi_dot = strtold(argv[12], NULL);
}

real mod2_v (series4 x, series4 v, parameters p) {
    matrix4x4 g;
    real_metric (g, p.m, p.a, x[1][0], x[2][0]);
    real v_dot_v = 0.0L;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            v_dot_v += g[mu][nu] * v[mu][0] * v[nu][0];
        }
    }
    return v_dot_v;
}

real error (real e) {
    return 10.0L * log10l(fabsl(e) >= 1e-36L ? fabsl(e) : 1e-36L);
}

void gr_output (long dp, series4 x, series4 v, real t, parameters p) {
    char fs[256];
    //sprintf(fs, "%%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.6Le\n",
            //dp, dp, dp, dp, dp, dp, dp, dp);
    //printf(fs, x[0], x[1], x[2], x[3], v[0], v[1], v[2], v[3], t);
    real v_dot_v = mod2_v(x, v, p);
    real ev_dot_v = error(p.v0 - v_dot_v);
    real ra = sqrtl(x[1][0] * x[1][0] + p.a * p.a);
    real sth = sinl(x[2][0]);
    sprintf(fs, "%%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.6Le %%+.3Le %%+.3Le\n", dp, dp, dp);
    printf(fs, ra * sth * cosl(x[3][0]), ra * sth * sinl(x[3][0]), x[1][0] * cosl(x[2][0]), t, v_dot_v, ev_dot_v);
}

void real_metric (matrix4x4 metric, real m, real a, real r, real theta) {
    metric[0][0] = g_t_t(m, a, d_dual(r), d_dual(theta)).val;
    metric[0][1] = metric[1][0] = g_t_r(m, a, d_dual(r), d_dual(theta)).val;
    metric[0][2] = metric[2][0] = g_t_theta(m, a, d_dual(r), d_dual(theta)).val;
    metric[0][3] = metric[3][0] = g_t_phi(m, a, d_dual(r), d_dual(theta)).val;
    metric[1][1] = g_r_r(m, a, d_dual(r), d_dual(theta)).val;
    metric[1][2] = metric[2][1] = g_r_theta(m, a, d_dual(r), d_dual(theta)).val;
    metric[1][3] = metric[3][1] = g_r_phi(m, a, d_dual(r), d_dual(theta)).val;
    metric[2][2] = g_theta_theta(m, a, d_dual(r), d_dual(theta)).val;
    metric[2][3] = metric[3][2] = g_theta_phi(m, a, d_dual(r), d_dual(theta)).val;
    metric[3][3] = g_phi_phi(m, a, d_dual(r), d_dual(theta)).val;
}

void real_inverse (matrix4x4 inverse, real m, real a, real r, real theta) {
    inverse[0][0] = i_t_t(m, a, r, theta);
    inverse[0][1] = inverse[1][0] = i_t_r(m, a, r, theta);
    inverse[0][2] = inverse[2][0] = i_t_theta(m, a, r, theta);
    inverse[0][3] = inverse[3][0] = i_t_phi(m, a, r, theta);
    inverse[1][1] = i_r_r(m, a, r, theta);
    inverse[1][2] = inverse[2][1] = i_r_theta(m, a, r, theta);
    inverse[1][3] = inverse[3][1] = i_r_phi(m, a, r, theta);
    inverse[2][2] = i_theta_theta(m, a, r, theta);
    inverse[2][3] = inverse[3][2] = i_theta_phi(m, a, r, theta);
    inverse[3][3] = i_phi_phi(m, a, r, theta);
}

void christoffel (matrix4x4x4 symbols, matrix4x4 inverse, real m, real a, real r, real theta) {
    matrix4x4 dr;
    dual r_dual = d_var(r);
    dual theta_dual = d_dual(theta);
    dr[0][0] = g_t_t(m, a, r_dual, theta_dual).dot;
    dr[0][1] = dr[1][0] = g_t_r(m, a, r_dual, theta_dual).dot;
    dr[0][2] = dr[2][0] = g_t_theta(m, a, r_dual, theta_dual).dot;
    dr[0][3] = dr[3][0] = g_t_phi(m, a, r_dual, theta_dual).dot;
    dr[1][1] = g_r_r(m, a, r_dual, theta_dual).dot;
    dr[1][2] = dr[2][1] = g_r_theta(m, a, r_dual, theta_dual).dot;
    dr[1][3] = dr[3][1] = g_r_phi(m, a, r_dual, theta_dual).dot;
    dr[2][2] = g_theta_theta(m, a, r_dual, theta_dual).dot;
    dr[2][3] = dr[3][2] = g_theta_phi(m, a, r_dual, theta_dual).dot;
    dr[3][3] = g_phi_phi(m, a, r_dual, theta_dual).dot;
    matrix4x4 dtheta;
    r_dual = d_dual(r);
    theta_dual = d_var(theta);
    dtheta[0][0] = g_t_t(m, a, r_dual, theta_dual).dot;
    dtheta[0][1] = dtheta[1][0] = g_t_r(m, a, r_dual, theta_dual).dot;
    dtheta[0][2] = dtheta[2][0] = g_t_theta(m, a, r_dual, theta_dual).dot;
    dtheta[0][3] = dtheta[3][0] = g_t_phi(m, a, r_dual, theta_dual).dot;
    dtheta[1][1] = g_r_r(m, a, r_dual, theta_dual).dot;
    dtheta[1][2] = dtheta[2][1] = g_r_theta(m, a, r_dual, theta_dual).dot;
    dtheta[1][3] = dtheta[3][1] = g_r_phi(m, a, r_dual, theta_dual).dot;
    dtheta[2][2] = g_theta_theta(m, a, r_dual, theta_dual).dot;
    dtheta[2][3] = dtheta[3][2] = g_theta_phi(m, a, r_dual, theta_dual).dot;
    dtheta[3][3] = g_phi_phi(m, a, r_dual, theta_dual).dot;
    matrix4x4x4 d_g;
    for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
            d_g[j][k][0] = 0.0L;
            d_g[j][k][1] = dr[j][k];
            d_g[j][k][2] = dtheta[j][k];
            d_g[j][k][3] = 0.0L;
        }
    }
    for (int i = 0; i < 4; i++) {
        for (int k = 0; k < 4; k++) {
            for (int l = 0; l < 4; l++) {
                real sum = 0.0L;
                for (int m = 0; m < 4; m++) {
                    sum += inverse[i][m] * (d_g[m][k][l] + d_g[m][l][k] - d_g[k][l][m]);
                }
                symbols[i][k][l] = 0.5L * sum;
            }
        }
    }
}

void geodesic (matrix4x4x4 symbols, series4 x_dot, series4 v_dot, series4 v, int k) {
    for (int mu = 0; mu < 4; mu++) {
        real sum = 0.0L;
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                sum -= symbols[mu][a][b] * t_prod(v[a], v[b], k);
            }
        }
        v_dot[mu][k] = sum;
        x_dot[mu][k] = v[mu][k];
    }
}

void tsm4 (int argc, char **argv) {
    long n, steps, dp;
    real h;
    vector4 x, v, xdot, vdot;
    gr_control(argv, &dp, &n, &h, &steps, &x[0], &x[1], &x[2], &x[3], &v[0], &v[1], &v[2], &v[3]);
    series4 x4 = t_jet4(n, x), v4 = t_jet4(n, v), xdot4 = t_jet4(n, xdot), vdot4 = t_jet4(n, vdot);
    parameters p = (parameters) {
        .m = strtold(argv[13], NULL),
        .a = strtold(argv[14], NULL)
    };
    p.v0 = mod2_v(x4, v4, p);
    gr_output(dp, x4, v4, 0.0, p);
    for (long step = 1; step <= steps; step++) {
        matrix4x4 inverse;
        real_inverse(inverse, p.m, p.a, x4[1][0], x4[2][0]);
        matrix4x4x4 symbols;
        christoffel(symbols, inverse, p.m, p.a, x4[1][0], x4[2][0]);
        for (int k = 0; k < n; k++) {
            geodesic(symbols, xdot4, vdot4, v4, k);
            for (int mu = 0; mu < 4; mu++) {
                x4[mu][k + 1] = xdot4[mu][k] / (k + 1);
                v4[mu][k + 1] = vdot4[mu][k] / (k + 1);
            }
        }
        for (int mu = 0; mu < 4; mu++) {
            *x4[mu] = t_horner(x4[mu], n, h);
            *v4[mu] = t_horner(v4[mu], n, h);
        }
        gr_output(dp, x4, v4, h * step, p);
    }
}
